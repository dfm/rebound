/**
 * @file 	problem.c
 * @brief 	Example problem: circular orbit.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * @detail 	This example uses the Wisdom Holman integrator
 * to integrate the outer planets of the solar system. The initial 
 * conditions are taken from Applegate et al 1986. Pluto is a test
 * particle.
 * 
 * @section 	LICENSE
 * Copyright (c) 2011 Hanno Rein, Shangfei Liu
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
#include <gmp.h>
#include <time.h>
#include <sys/time.h>
#include "main.h"
#include "output.h"
#include "input.h"
#include "integrator.h"
#include "tools.h"
#include "particle.h"
#include "boundaries.h"

double ss_pos[6][3] = 
{
	{-4.06428567034226e-3,	-6.08813756435987e-3,	-1.66162304225834e-6	}, // Sun
	{+3.40546614227466e+0,	+3.62978190075864e+0,	+3.42386261766577e-2	}, // Jupiter
	{+6.60801554403466e+0,	+6.38084674585064e+0,	-1.36145963724542e-1	}, // Saturn
	{+1.11636331405597e+1,	+1.60373479057256e+1,	+3.61783279369958e-1	}, // Uranus
	{-3.01777243405203e+1,	+1.91155314998064e+0,	-1.53887595621042e-1	}, // Neptune
	{-2.13858977531573e+1,	+3.20719104739886e+1,	+2.49245689556096e+0	}  // Pluto
};
double ss_vel[6][3] = 
{
	{+6.69048890636161e-6,	-6.33922479583593e-6,	-3.13202145590767e-9	}, // Sun
	{-5.59797969310664e-3,	+5.51815399480116e-3,	-2.66711392865591e-6	}, // Jupiter
	{-4.17354020307064e-3,	+3.99723751748116e-3,	+1.67206320571441e-5	}, // Saturn
	{-3.25884806151064e-3,	+2.06438412905916e-3,	-2.17699042180559e-5	}, // Uranus
	{-2.17471785045538e-4,	-3.11361111025884e-3,	+3.58344705491441e-5	}, // Neptune
	{-1.76936577252484e-3,	-2.06720938381724e-3,	+6.58091931493844e-4	}  // Pluto
};

double ss_mass[6] =
{
	1.00000597682, 	// Sun + inner planets
	1./1047.355,	// Jupiter
	1./3501.6,	// Saturn
	1./22869.,	// Uranus
	1./19314.,	// Neptune
	0.		// Pluto
};

const double k	 	= 0.01720209895;	// Gaussian constant 
#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

double energy();
const double jupiter_orbital_period = 4332.589;
double timing_start;

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= input_get_double(argc,argv,"dt",0.1);			// days
	N_active	= 5;
	tmax		= jupiter_orbital_period*input_get_double(argc,argv,"orbits",1e6);	// number of orbits
	G		= k*k;
#ifdef INTEGRATOR_IAS15
	integrator_epsilon = input_get_double(argc,argv,"integrator_epsilon",1e-3);
#endif // INTEGRATOR_IAS15

#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(200); 			// Init box with width 200 astronomical units

	// Initial conditions
	for (int i=0;i<6;i++){
		struct particle p;
		p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
		p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
		p.ax = 0; 			p.ay = 0; 			p.az = 0;
		p.m  = ss_mass[i];
		particles_add(p); 
	}
#ifdef INTEGRATOR_WH
	// Move to heliocentric frame (required by WHM)
	for (int i=1;i<N;i++){
		particles[i].x -= particles[0].x;	particles[i].y -= particles[0].y;	particles[i].z -= particles[0].z;
		particles[i].vx -= particles[0].vx;	particles[i].vy -= particles[0].vy;	particles[i].vz -= particles[0].vz;
	}
	particles[0].x = 0;	particles[0].y = 0;	particles[0].z = 0;
	particles[0].vx= 0;	particles[0].vy= 0;	particles[0].vz= 0;
#else
	// Move to barycentric frame
	tools_move_to_center_of_momentum();
#endif // INTEGRATOR_WH
	mpf_set_default_prec(512);
	energy(); // calculate and store initial energy
	
	struct timeval tim;
	gettimeofday(&tim, NULL);
	timing_start = tim.tv_sec+(tim.tv_usec/1000000.0);
}

void problem_inloop(){
}

int energy_init_set = 0;
mpf_t energy_init;
double energy(){
	mpf_t energy_kinetic;
	mpf_init(energy_kinetic);
	mpf_t energy_potential;
	mpf_init(energy_potential);
	mpf_t mass;
	mpf_init(mass);
	mpf_t vx;
	mpf_init(vx);
	mpf_t vy;
	mpf_init(vy);
	mpf_t vz;
	mpf_init(vz);
	mpf_t temp;
	mpf_init(temp);
	for (int i=0;i<N;i++){
		mpf_set_d(temp,particles[i].vx*particles[i].m);
		mpf_add(vx, vx, temp);
		mpf_set_d(temp,particles[i].vy*particles[i].m);
		mpf_add(vy, vy, temp);
		mpf_set_d(temp,particles[i].vz*particles[i].m);
		mpf_add(vz, vz, temp);

		mpf_set_d(temp,particles[i].m);
		mpf_add(mass, mass, temp);
	}
	mpf_div(vx,vx,mass);
	mpf_div(vy,vy,mass);
	mpf_div(vz,vz,mass);

	for (int i=0;i<N;i++){
		for (int j=0;j<i;j++){
			double dx = particles[i].x - particles[j].x;
			double dy = particles[i].y - particles[j].y;
			double dz = particles[i].z - particles[j].z;
			double r = sqrt(dx*dx + dy*dy + dz*dz + softening*softening);
			mpf_set_d(temp,-G*particles[i].m*particles[j].m/r);
			mpf_add(energy_potential, energy_potential,temp);
			
		}
	
		double dvx = particles[i].vx-mpf_get_d(vx);
		double dvy = particles[i].vy-mpf_get_d(vy);
		double dvz = particles[i].vz-mpf_get_d(vz);
		mpf_set_d(temp,1./2.*particles[i].m * (dvx*dvx + dvy*dvy + dvz*dvz));
		mpf_add(energy_kinetic, energy_kinetic,temp);
		
	}
	mpf_add(energy_kinetic,energy_kinetic,energy_potential);
	if (!energy_init_set){
		mpf_init(energy_init);
		mpf_set(energy_init,energy_kinetic);
		energy_init_set = 1;
	}
	mpf_sub(energy_kinetic,energy_kinetic,energy_init);
	mpf_div(energy_kinetic,energy_kinetic,energy_init);
	mpf_abs(energy_kinetic,energy_kinetic);
	
	double return_value = mpf_get_d(energy_kinetic);
	mpf_clear(energy_kinetic);
	mpf_clear(temp);
	mpf_clear(energy_potential);
	mpf_clear(mass);
	mpf_clear(vx);
	mpf_clear(vy);
	mpf_clear(vz);
	return return_value;
}

void problem_finish(){
}

double output_interval 	= 1.01;
double output_next	= 100;
void problem_output(){
	if (output_next<t){
		output_next *= output_interval;
		char filename[4096];
		sprintf(filename,"energy_%s.txt",input_arguments);
		FILE* of = fopen(filename,"a+"); 
		struct timeval tim;
		gettimeofday(&tim, NULL);
		double timing =  (tim.tv_sec+(tim.tv_usec/1000000.0)) - timing_start;
		fprintf(of,"%e\t%e\t%e\n",t/jupiter_orbital_period,energy(),timing);
		fclose(of);
	}

}
