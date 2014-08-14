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

double energy();
const double jupiter_orbital_period = 2.*M_PI;
double timing_start;

void harmonic(){
	particles[0].ax += -particles[0].x/particles[0].m;
}

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= input_get_double(argc,argv,"dt",0.1);			// days
	tmax		= jupiter_orbital_period*input_get_double(argc,argv,"orbits",1e6);	// number of orbits
	problem_additional_forces = harmonic();
	G		= k*k;
#ifdef INTEGRATOR_IAS15
	integrator_epsilon = input_get_double(argc,argv,"integrator_epsilon",1e-3);
#endif // INTEGRATOR_IAS15

	init_boxwidth(200); 			// Init box with width 200 astronomical units

	// Initial conditions
	for (int i=0;i<1;i++){
		struct particle p;
		p.x  = 0; 		p.y  = 0;	 	p.z  = 0;
		p.vx = 1; 		p.vy = 0;	 	p.vz = 0;
		p.ax = 0; 		p.ay = 0; 		p.az = 0;
		p.m  = 1;
		particles_add(p); 
	}
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
		double dx = particles[i].x;
		mpf_set_d(temp,0.5*particles[i].m*dx*dx);
		mpf_add(energy_potential, energy_potential,temp);
			
	
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
