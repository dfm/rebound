/**
 * @file 	integrator.c
 * @brief 	IAS15 integrator.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail	This file implements the IAS15 integration scheme.  
 * IAS stands for Integrator with Adaptive Step-size control, 15th 
 * order. This scheme is a fifteenth order integrator well suited for 
 * high accuracy orbit integration with non-conservative forces.
 * For more details see Rein & Spiegel 2014. Also see Everhart, 1985,
 * ASSL Vol. 115, IAU Colloq. 83, Dynamics of Comets, Their Origin 
 * and Evolution, 185 for the original implementation by Everhart.
 * Part of this code is based a function from the ORSE package.
 * See orsa.sourceforge.net for more details on their implementation.
 *
 * 
 * @section 	LICENSE
 * Copyright (c) 2011-2012 Hanno Rein, Dave Spiegel.
 * Copyright (c) 2002-2004 Pasquale Tricarico.
 *
 * This file is part of rebound.
 *
 * rebound is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * rebound is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with rebound.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "particle.h"
#include "main.h"
#include "gravity.h"
#include "boundaries.h"
#include "problem.h"
#include "output.h"
#include "tools.h"

#ifdef TREE
#error IAS15 integrator not working with TREE module.
#endif
#ifdef MPI
#error IAS15 integrator not working with MPI.
#endif

int 	integrator_force_is_velocitydependent	= 1;	// Turn this off to safe some time if the force is not velocity dependent.
double 	integrator_epsilon 			= 5e-9;	// Precision parameter 
							// If it is zero, then a constant timestep is used. 
int	integrator_epsilon_global		= 1;	// if 1: estimate the fractional error by max(acceleration_error)/max(acceleration), where max is take over all particles.
							// if 0: estimate the fractional error by max(acceleration_error/acceleration).
double 	integrator_min_dt 			= 0;	// Minimum timestep used as a floor when adaptive timestepping is enabled.
double	integrator_error			= 0;	// Error estimate in last timestep (used for debugging only)
unsigned int integrator_iterations_max		= 12;	// Maximum number of iterations in predictor/corrector loop
unsigned long integrator_iterations_max_exceeded= 0;	// Count how many times the iteration did not converge
const double safety_factor 			= 0.25;	// Maximum increase/deacrease of consecutve timesteps.


const double h[9]	= { 
	0.0, 
	0.04463395528996985073312102185830776152357713301,
	0.1443662570421455714852185202282149697135218395,
	0.2868247571444305189486862397490926585531162091,
	0.4548133151965733509677277700467869448670199439,
	0.6280678354167276975691460395173707553318999366,
	0.7856915206043692416424587324183296773255702825,
	0.9086763921002060439962585419254586468929595243,
	0.9822200848526365481867948989623209387335116018
	}; // Gauss Lobatto spacings

double r[36],c[28],d[28],s[10]; // These constants will be set dynamically.

int N3allocated 		= 0; 	// Size of allocated arrays.
int integrator_ias15_init_done 	= 0;	// Calculate coefficients once.

double* at   = NULL;	// Temporary buffer for acceleration
double* x0  = NULL;	// Temporary buffer for position (used for initial values at h=0) 
double* v0  = NULL;	//                      velocity
double* a0  = NULL;	//                      acceleration
double* csx  = NULL;	//                      compensated summation
double* csv  = NULL;	//                      compensated summation

double* g[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double* b[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double* e[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;
double* br[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;	// Used for resetting after timestep rejection
double* er[8] = {NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL} ;

void copybuffers(double* _a[8], double* _b[8], int N3);
void predict_next_step(double ratio, int N3, double* _e[8], double* _b[8]);
double dt_last_success;

void integrator_part1(){
	// Do nothing here. This is only used in a leapfrog-like DKD integrator.
}

// This function updates the acceleration on all particles. 
// It uses the current position and velocity data in the (struct particle*) particles structure.
// Note: this does currently not work with MPI or any TREE module.
void integrator_update_acceleration(){
	PROFILING_STOP(PROFILING_CAT_INTEGRATOR)
	PROFILING_START()
	gravity_calculate_acceleration();
	if (problem_additional_forces) problem_additional_forces();
	PROFILING_STOP(PROFILING_CAT_GRAVITY)
	PROFILING_START()
}

int integrator_ias15_step(); // Does the actual timestep.

void integrator_part2(){
	if (!integrator_ias15_init_done){ 	// Generate coefficients.
		int l=0;
		for (int j=1;j<9;++j) {
			for(int k=0;k<j;++k) {
				r[l] = 1.0 / (h[j] - h[k]);
				++l;
			}
		}
		c[0] = -h[1];
		d[0] =  h[1];
		l=0;
		for (int j=2;j<8;++j) {
			++l;
			c[l] = -h[j] * c[l-j+1];
			d[l] =  h[1] * d[l-j+1];
			for(int k=2;k<j;++k) {
				++l;
				c[l] = c[l-j] - h[j] * c[l-j+1];
				d[l] = d[l-j] + h[k] * d[l-j+1];
			}
			++l;
			c[l] = c[l-j] - h[j];
			d[l] = d[l-j] + h[j]; 
		}
		integrator_ias15_init_done = 1;
	}
	// Try until a step was successful.
	while(!integrator_ias15_step());
}
 
int integrator_ias15_step() {
	const int N3 = 3*N;
	if (N3 > N3allocated) {
		for (int l=0;l<8;++l) {
			g[l] = realloc(g[l],sizeof(double)*N3);
			b[l] = realloc(b[l],sizeof(double)*N3);
			e[l] = realloc(e[l],sizeof(double)*N3);
			br[l] = realloc(br[l],sizeof(double)*N3);
			er[l] = realloc(er[l],sizeof(double)*N3);
			for (int k=0;k<N3;k++){
				b[l][k] = 0;
				e[l][k] = 0;
				br[l][k] = 0;
				er[l][k] = 0;
			}
		}
		at = realloc(at,sizeof(double)*N3);
		x0 = realloc(x0,sizeof(double)*N3);
		v0 = realloc(v0,sizeof(double)*N3);
		a0 = realloc(a0,sizeof(double)*N3);
		csx= realloc(csx,sizeof(double)*N3);
		csv= realloc(csv,sizeof(double)*N3);
		for (int i=0;i<N3;i++){
			// Kill compensated summation coefficients
			csx[i] = 0;
			csv[i] = 0;
		}
		N3allocated = N3;
	}
	
	// integrator_update_acceleration(); // Not needed. Forces are already calculated in main routine.

	for(int k=0;k<N;k++) {
		x0[3*k]   = particles[k].x;
		x0[3*k+1] = particles[k].y;
		x0[3*k+2] = particles[k].z;
		v0[3*k]   = particles[k].vx;
		v0[3*k+1] = particles[k].vy;
		v0[3*k+2] = particles[k].vz;
		a0[3*k]   = particles[k].ax;
		a0[3*k+1] = particles[k].ay;  
		a0[3*k+2] = particles[k].az;
	}

	for(int k=0;k<N3;k++) {
		g[0][k] = b[7][k]*d[21] + b[6][k]*d[15] + b[5][k]*d[10] + b[4][k]*d[6] + b[3][k]*d[3]  + b[2][k]*d[1]  + b[1][k]*d[0]  + b[0][k];
		g[1][k] = b[7][k]*d[22] + b[6][k]*d[16] + b[5][k]*d[11] + b[4][k]*d[7] + b[3][k]*d[4]  + b[2][k]*d[2]  + b[1][k];
		g[2][k] = b[7][k]*d[23] + b[6][k]*d[17] + b[5][k]*d[12] + b[4][k]*d[8] + b[3][k]*d[5]  + b[2][k];
		g[3][k] = b[7][k]*d[24] + b[6][k]*d[18] + b[5][k]*d[13] + b[4][k]*d[9] + b[3][k];
		g[4][k] = b[7][k]*d[25] + b[6][k]*d[19] + b[5][k]*d[14] + b[4][k];
		g[5][k] = b[7][k]*d[26] + b[6][k]*d[20] + b[5][k];
		g[6][k] = b[7][k]*d[27] + b[6][k];
		g[7][k] = b[7][k];
	}

	double predictor_corrector_error = 1;
	double predictor_corrector_error_last = 2;
	int iterations = 0;	
	// Predictor corrector loop
	// Stops if accuracy better than 1e-16 or if accuracy starts to oscillate
	while(predictor_corrector_error>1e-16 ){	
		predictor_corrector_error_last = predictor_corrector_error;
		if (iterations>=integrator_iterations_max){
			integrator_iterations_max_exceeded++;
			const int integrator_iterations_warning = 10;
			if (integrator_iterations_max_exceeded==integrator_iterations_warning ){
				fprintf(stderr,"\n\033[1mWarning!\033[0m At least %d predictor corrector loops in integrator_ias15.c did not converge. This is typically an indication of the timestep being too large.\n",integrator_iterations_warning);
			}
			break;								// Quit predictor corrector loop
		}
		predictor_corrector_error = 0;
		iterations++;

		for(int n=1;n<9;n++) {							// Loop over interval using Gauss-Radau spacings

			s[0] = dt * h[n];
			s[1] = s[0] * s[0] / 2.;
			s[2] = s[1] * h[n] / 3.;
			s[3] = s[2] * h[n] / 2.;
			s[4] = 3. * s[3] * h[n] / 5.;
			s[5] = 2. * s[4] * h[n] / 3.;
			s[6] = 5. * s[5] * h[n] / 7.;
			s[7] = 3. * s[6] * h[n] / 4.;
			s[8] = 7. * s[7] * h[n] / 9.;
			s[9] = 8. * s[8] * h[n] / 10.;

			// Prepare particles arrays for force calculation
			for(int i=0;i<N;i++) {						// Predict positions at interval n using b values
				const int k0 = 3*i+0;
				const int k1 = 3*i+1;
				const int k2 = 3*i+2;

				double xk0  = csx[k0] + (s[9]*b[7][k0] + s[8]*b[6][k0] + s[7]*b[5][k0] + s[6]*b[4][k0] + s[5]*b[3][k0] + s[4]*b[2][k0] + s[3]*b[1][k0] + s[2]*b[0][k0] + s[1]*a0[k0] + s[0]*v0[k0] );
				particles[i].x = xk0 + x0[k0];           
				double xk1  = csx[k1] + (s[9]*b[7][k1] + s[8]*b[6][k1] + s[7]*b[5][k1] + s[6]*b[4][k1] + s[5]*b[3][k1] + s[4]*b[2][k1] + s[3]*b[1][k1] + s[2]*b[0][k1] + s[1]*a0[k1] + s[0]*v0[k1] );
				particles[i].y = xk1 + x0[k1];          
				double xk2  = csx[k2] + (s[9]*b[7][k2] + s[8]*b[6][k2] + s[7]*b[5][k2] + s[6]*b[4][k2] + s[5]*b[3][k2] + s[4]*b[2][k2] + s[3]*b[1][k2] + s[2]*b[0][k2] + s[1]*a0[k2] + s[0]*v0[k2] );
				particles[i].z = xk2 + x0[k2];
			}
			
			if (integrator_force_is_velocitydependent){
				s[0] = dt * h[n];
				s[1] =      s[0] * h[n] / 2.;
				s[2] = 2. * s[1] * h[n] / 3.;
				s[3] = 3. * s[2] * h[n] / 4.;
				s[4] = 4. * s[3] * h[n] / 5.;
				s[5] = 5. * s[4] * h[n] / 6.;
				s[6] = 6. * s[5] * h[n] / 7.;
				s[7] = 7. * s[6] * h[n] / 8.;
				s[8] = 8. * s[7] * h[n] / 9.;

				for(int i=0;i<N;i++) {					// Predict velocities at interval n using b values
					const int k0 = 3*i+0;
					const int k1 = 3*i+1;
					const int k2 = 3*i+2;

					double vk0 =  csv[k0] + s[8]*b[7][k0] + s[7]*b[6][k0] + s[6]*b[5][k0] + s[5]*b[4][k0] + s[4]*b[3][k0] + s[3]*b[2][k0] + s[2]*b[1][k0] + s[1]*b[0][k0] + s[0]*a0[k0];
					particles[i].vx = vk0 + v0[k0];
					double vk1 =  csv[k1] + s[8]*b[7][k1] + s[7]*b[6][k1] + s[6]*b[5][k1] + s[5]*b[4][k1] + s[4]*b[3][k1] + s[3]*b[2][k1] + s[2]*b[1][k1] + s[1]*b[0][k1] + s[0]*a0[k1];
					particles[i].vy = vk1 + v0[k1];
					double vk2 =  csv[k2] + s[8]*b[7][k2] + s[7]*b[6][k2] + s[6]*b[5][k2] + s[5]*b[4][k2] + s[4]*b[3][k2] + s[3]*b[2][k2] + s[2]*b[1][k2] + s[1]*b[0][k2] + s[0]*a0[k2];
					particles[i].vz = vk2 + v0[k2];
				}
			}


			integrator_update_acceleration();				// Calculate forces at interval n

			for(int k=0;k<N;++k) {
				at[3*k]   = particles[k].ax;
				at[3*k+1] = particles[k].ay;  
				at[3*k+2] = particles[k].az;
			}
			switch (n) {							// Improve b and g values
				case 1: 
					for(int k=0;k<N3;++k) {
						double tmp = g[0][k];
						g[0][k]  = (at[k] - a0[k]) * r[0];
						b[0][k] += g[0][k] - tmp;
					} break;
				case 2: 
					for(int k=0;k<N3;++k) {
						double tmp = g[1][k];
						double gk = at[k] - a0[k];
						g[1][k] = (gk*r[1] - g[0][k])*r[2];
						tmp = g[1][k] - tmp;
						b[0][k] += tmp * c[0];
						b[1][k] += tmp;
					} break;
				case 3: 
					for(int k=0;k<N3;++k) {
						double tmp = g[2][k];
						double gk = at[k] - a0[k];
						g[2][k] = ((gk*r[3] - g[0][k])*r[4] - g[1][k])*r[5];
						tmp = g[2][k] - tmp;
						b[0][k] += tmp * c[1];
						b[1][k] += tmp * c[2];
						b[2][k] += tmp;
					} break;
				case 4:
					for(int k=0;k<N3;++k) {
						double tmp = g[3][k];
						double gk = at[k] - a0[k];
						g[3][k] = (((gk*r[6] - g[0][k])*r[7] - g[1][k])*r[8] - g[2][k])*r[9];
						tmp = g[3][k] - tmp;
						b[0][k] += tmp * c[3];
						b[1][k] += tmp * c[4];
						b[2][k] += tmp * c[5];
						b[3][k] += tmp;
					} break;
				case 5:
					for(int k=0;k<N3;++k) {
						double tmp = g[4][k];
						double gk = at[k] - a0[k];
						g[4][k] = ((((gk*r[10] - g[0][k])*r[11] - g[1][k])*r[12] - g[2][k])*r[13] - g[3][k])*r[14];
						tmp = g[4][k] - tmp;
						b[0][k] += tmp * c[6];
						b[1][k] += tmp * c[7];
						b[2][k] += tmp * c[8];
						b[3][k] += tmp * c[9];
						b[4][k] += tmp;
					} break;
				case 6:
					for(int k=0;k<N3;++k) {
						double tmp = g[5][k];
						double gk = at[k] - a0[k];
						g[5][k] = (((((gk*r[15] - g[0][k])*r[16] - g[1][k])*r[17] - g[2][k])*r[18] - g[3][k])*r[19] - g[4][k])*r[20];
						tmp = g[5][k] - tmp;
						b[0][k] += tmp * c[10];
						b[1][k] += tmp * c[11];
						b[2][k] += tmp * c[12];
						b[3][k] += tmp * c[13];
						b[4][k] += tmp * c[14];
						b[5][k] += tmp;
					} break;
				case 7:
				{
					for(int k=0;k<N3;++k) {
						double tmp = g[6][k];
						double gk = at[k] - a0[k];
						g[6][k] = ((((((gk*r[21] - g[0][k])*r[22] - g[1][k])*r[23] - g[2][k])*r[24] - g[3][k])*r[25] - g[4][k])*r[26] - g[5][k])*r[27];
						tmp = g[6][k] - tmp;	
						b[0][k] += tmp * c[15];
						b[1][k] += tmp * c[16];
						b[2][k] += tmp * c[17];
						b[3][k] += tmp * c[18];
						b[4][k] += tmp * c[19];
						b[5][k] += tmp * c[20];
						b[6][k] += tmp;
					} break;
				case 8:
				{
					double maxak = 0.0;
					double maxb6ktmp = 0.0;
					for(int k=0;k<N3;++k) {
						double tmp = g[7][k];
						double gk = at[k] - a0[k];
						g[7][k] = (((((((gk*r[28] - g[0][k])*r[29] - g[1][k])*r[30] - g[2][k])*r[31] - g[3][k])*r[32] - g[4][k])*r[33] - g[5][k])*r[34] - g[6][k])*r[35];
						tmp = g[7][k] - tmp;	
						b[0][k] += tmp * c[21];
						b[1][k] += tmp * c[22];
						b[2][k] += tmp * c[23];
						b[3][k] += tmp * c[24];
						b[4][k] += tmp * c[25];
						b[5][k] += tmp * c[26];
						b[6][k] += tmp * c[27];
						b[7][k] += tmp;
						
						// Monitor change in b[6][k] relative to at[k]. The predictor corrector scheme is converged if it is close to 0.
						if (integrator_epsilon_global){
							const double ak  = fabs(at[k]);
							if (isnormal(ak) && ak>maxak){
								maxak = ak;
							}
							const double b6ktmp = fabs(tmp);  // change of b6ktmp coefficient
							if (isnormal(b6ktmp) && b6ktmp>maxb6ktmp){
								maxb6ktmp = b6ktmp;
							}
						}else{
							const double ak  = at[k];
							const double b6ktmp = tmp; 
							const double errork = fabs(b6ktmp/ak);
							if (isnormal(errork) && errork>predictor_corrector_error){
								predictor_corrector_error = errork;
							}
						}
					} 
					if (integrator_epsilon_global){
						predictor_corrector_error = maxb6ktmp/maxak;
					}
					
					}
					break;
				}
			}
		}
	}
	// Find new timestep
	const double dt_done = dt;
	
	if (integrator_epsilon>0){
		// Estimate error (given by last term in series expansion) 
		// There are two options:
		// integrator_epsilon_global==1  (default)
		//   First, we determine the maximum acceleration and the maximum of the last term in the series. 
		//   Then, the two are divided.
		// integrator_epsilon_global==0
		//   Here, the fractional error is calculated for each particle individually and we use the maximum of the fractional error.
		//   This might fail in cases where a particle does not experience any (physical) acceleration besides roundoff errors. 
		integrator_error = 0.0;
		if (integrator_epsilon_global){
			double maxak = 0.0;
			double maxb6k = 0.0;
			for(int k=0;k<N3;k++) {  // Looping over all particles and all 3 components of the acceleration. 
				const double ak  = fabs(at[k]);
				if (isnormal(ak) && ak>maxak){
					maxak = ak;
				}
				const double b6k = fabs(b[7][k]); 
				if (isnormal(b6k) && b6k>maxb6k){
					maxb6k = b6k;
				}
			}
			integrator_error = maxb6k/maxak;
		}else{
			for(int k=0;k<N3;k++) {
				const double ak  = at[k];
				const double b6k = b[7][k]; 
				const double errork = fabs(b6k/ak);
				if (isnormal(errork) && errork>integrator_error){
					integrator_error = errork;
				}
			}
		}

		double dt_new;
		if  (isnormal(integrator_error)){ 	
			// if error estimate is available increase by more educated guess
		 	dt_new = pow(integrator_epsilon/integrator_error,1./8.)*dt_done;
		}else{					// In the rare case that the error estimate doesn't give a finite number (e.g. when all forces accidentally cancel up to machine precission).
		 	dt_new = dt_done/safety_factor; // by default, increase timestep a little
		}
		
		if (dt_new<integrator_min_dt) dt_new = integrator_min_dt;
		
		if (fabs(dt_new/dt_done) < safety_factor) {	// New timestep is significantly smaller.
			// Reset particles
			for(int k=0;k<N;++k) {
				particles[k].x = x0[3*k+0];	// Set inital position
				particles[k].y = x0[3*k+1];
				particles[k].z = x0[3*k+2];

				particles[k].vx = v0[3*k+0];	// Set inital velocity
				particles[k].vy = v0[3*k+1];
				particles[k].vz = v0[3*k+2];
			}
			dt = dt_new;
			double ratio = dt/dt_last_success;
			predict_next_step(ratio, N3, er, br);
			
			return 0; // Step rejected. Do again. 
		}		
		if (fabs(dt_new/dt_done) > 1.0) {	// New timestep is larger.
			if (dt_new/dt_done > 1./safety_factor) dt_new = dt_done /safety_factor;	// Don't increase the timestep by too much compared to the last one.
		}
		dt = dt_new;
	}

	// Find new position and velocity values at end of the sequence
	const double dt_done2 = dt_done * dt_done;
	for(int k=0;k<N3;++k) {
		{
			double a = x0[k];
			csx[k]  +=  (b[7][k]/90. + b[6][k]/72. + b[5][k]/56. + b[4][k]/42. + b[3][k]/30. + b[2][k]/20. + b[1][k]/12. + b[0][k]/6. + a0[k]/2.) 
					* dt_done2 + v0[k] * dt_done;
			x0[k]    = a + csx[k];
			csx[k]  += a - x0[k]; 
		}
		{
			double a = v0[k]; 
			csv[k]  += (b[7][k]/9. + b[6][k]/8. + b[5][k]/7. + b[4][k]/6. + b[3][k]/5. + b[2][k]/4. + b[1][k]/3. + b[0][k]/2. + a0[k])
					* dt_done;
			v0[k]    = a + csv[k];
			csv[k]  += a - v0[k];
		}
	}

	t += dt_done;
	// Swap particle buffers

	for(int k=0;k<N;++k) {
		particles[k].x = x0[3*k+0];	// Set final position
		particles[k].y = x0[3*k+1];
		particles[k].z = x0[3*k+2];

		particles[k].vx = v0[3*k+0];	// Set final velocity
		particles[k].vy = v0[3*k+1];
		particles[k].vz = v0[3*k+2];
	}
	dt_last_success = dt_done;
	copybuffers(e,er,N3);		
	copybuffers(b,br,N3);		
	double ratio = dt/dt_done;
	predict_next_step(ratio, N3, e, b);
	return 1; // Success.
}

void predict_next_step(double ratio, int N3, double* _e[8], double* _b[8]){
	// Predict new B values to use at the start of the next sequence. The predicted
	// values from the last call are saved as E. The correction, BD, between the
	// actual and predicted values of B is applied in advance as a correction.
	//
	const double q1 = ratio;
	const double q2 = q1 * q1;
	const double q3 = q1 * q2;
	const double q4 = q2 * q2;
	const double q5 = q2 * q3;
	const double q6 = q3 * q3;
	const double q7 = q3 * q4;

	for(int k=0;k<N3;++k) {
		double be0 = _b[0][k] - _e[0][k];
		double be1 = _b[1][k] - _e[1][k];
		double be2 = _b[2][k] - _e[2][k];
		double be3 = _b[3][k] - _e[3][k];
		double be4 = _b[4][k] - _e[4][k];
		double be5 = _b[5][k] - _e[5][k];
		double be6 = _b[6][k] - _e[6][k];
		double be7 = _b[7][k] - _e[7][k];


		e[0][k] = q1*(_b[6][k]* 7.0 + _b[5][k]* 6.0 + _b[4][k]* 5.0 + _b[3][k]* 4.0 + _b[2][k]* 3.0 + _b[1][k]*2.0 + _b[0][k]);
		e[1][k] = q2*(_b[6][k]*21.0 + _b[5][k]*15.0 + _b[4][k]*10.0 + _b[3][k]* 6.0 + _b[2][k]* 3.0 + _b[1][k]);
		e[2][k] = q3*(_b[6][k]*35.0 + _b[5][k]*20.0 + _b[4][k]*10.0 + _b[3][k]* 4.0 + _b[2][k]);
		e[3][k] = q4*(_b[6][k]*35.0 + _b[5][k]*15.0 + _b[4][k]* 5.0 + _b[3][k]);
		e[4][k] = q5*(_b[6][k]*21.0 + _b[5][k]* 6.0 + _b[4][k]);
		e[5][k] = q6*(_b[6][k]* 7.0 + _b[5][k]);
		e[6][k] = q7* _b[6][k];
		e[7][k] = 0; // ignoring this for now TODO
		

		b[0][k] = e[0][k] + be0;
		b[1][k] = e[1][k] + be1;
		b[2][k] = e[2][k] + be2;
		b[3][k] = e[3][k] + be3;
		b[4][k] = e[4][k] + be4;
		b[5][k] = e[5][k] + be5;
		b[6][k] = e[6][k] + be6;
		b[7][k] = e[7][k] + be7;
	}
}

void copybuffers(double* _a[8], double* _b[8], int N3){
	for (int i=0;i<N3;i++){	
		_b[0][i] = _a[0][i];
		_b[1][i] = _a[1][i];
		_b[2][i] = _a[2][i];
		_b[3][i] = _a[3][i];
		_b[4][i] = _a[4][i];
		_b[5][i] = _a[5][i];
		_b[6][i] = _a[6][i];
		_b[7][i] = _a[7][i];
	}
}
