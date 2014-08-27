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
#include "main.h"
#include "output.h"
#include "input.h"
#include "integrator.h"
#include "tools.h"
#include "particle.h"
#include "boundaries.h"

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

double energy();
double* angular_momentum();
double energy_init;
double* aminit;
void input_binary_special(char* filename);
int init_N;
double timescale = 1;
int outputenergy;
double* jacobi();
double* jacobi_init;

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= input_get_double(argc,argv,"dt",1.);			// days
	outputenergy	= input_get_int(argc,argv,"outputenergy",0);
	const double k 	= 0.01720209895;	// Gaussian constant 
	G		= k*k;

	integrator_epsilon = input_get_double(argc,argv,"integrator_epsilon",integrator_epsilon);
	integrator_force_is_velocitydependent = 0;
	
#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(150); 			// Init box with width 200 astronomical units

	input_binary_special("particles.bin");
	dt *= timescale;
	init_N = N;

#ifndef INTEGRATOR_WH
	// Move to barycentric frame
	tools_move_to_center_of_momentum();
#endif // INTEGRATOR_WH
	mpf_set_default_prec(512);
	energy_init = energy();
	jacobi_init = jacobi();
	aminit = angular_momentum();
}

void problem_inloop(){
}

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
	for (int i=0;i<N;i++){
		mpf_t temp;
		mpf_init(temp);
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
			mpf_t temp;
			mpf_init(temp);
			mpf_set_d(temp,-G*particles[i].m*particles[j].m/r);
			mpf_add(energy_potential, energy_potential,temp);
			
		}
	
		double dvx = particles[i].vx-mpf_get_d(vx);
		double dvy = particles[i].vy-mpf_get_d(vy);
		double dvz = particles[i].vz-mpf_get_d(vz);
		mpf_t temp;
		mpf_init(temp);
		
		mpf_set_d(temp,1./2.*particles[i].m * (dvx*dvx + dvy*dvy + dvz*dvz));
		mpf_add(energy_kinetic, energy_kinetic,temp);
		
	}
	mpf_add(energy_kinetic,energy_kinetic,energy_potential);

	return mpf_get_d(energy_kinetic);
}

double* jacobi(){
	tools_move_to_center_of_momentum();

	double* j = malloc((N-2)*sizeof(double));
	double dx = particles[1].x - particles[0].x;
	double dy = particles[1].y - particles[0].y;
	double dz = particles[1].z - particles[0].z;
	double r = sqrt(dx*dx + dy*dy + dz*dz );
	
	double dvx = particles[1].vx - particles[0].vx;
	double dvy = particles[1].vy - particles[0].vy;
	double dvz = particles[1].vz - particles[0].vz;
	double v = sqrt(dvx*dvx + dvy*dvy + dvz*dvz );

	double mu = G*(particles[1].m+particles[0].m);
	double a = 1./(2./r-v*v/mu);
	double n = sqrt(mu/(a*a*a));

	for(int i=2;i<N;i++){
		double jac = 0;
		double r1,r2;
		{
			double dx = particles[i].x - particles[0].x;
			double dy = particles[i].y - particles[0].y;
			double dz = particles[i].z - particles[0].z;
			r1 = sqrt(dx*dx + dy*dy + dz*dz );
		}
		{
			double dx = particles[i].x - particles[1].x;
			double dy = particles[i].y - particles[1].y;
			double dz = particles[i].z - particles[1].z;
			r2 = sqrt(dx*dx + dy*dy + dz*dz );
		}
		jac += 2.*G*particles[0].m/r1;
		jac += 2.*G*particles[1].m/r2;
		jac += 2.*n*(particles[i].x*particles[i].vy-particles[i].y*particles[i].vx);
		jac -= particles[i].vx * particles[i].vx;
		jac -= particles[i].vy * particles[i].vy;
		jac -= particles[i].vz * particles[i].vz;
		j[i-2] = jac;
	}
	return j;
}

void problem_output(){
	if (outputenergy){
		if(output_check(tmax/10000.)){
			FILE* of = fopen("energy_timeseries.txt","a+"); 
			double rel_energy = fabs((energy()-energy_init)/energy_init);
			fprintf(of,"%e\t%e\n",t,rel_energy);
			fclose(of);
		}
	}
	{
	FILE* of = fopen("jacobi.txt","a+"); 
	double* jacobi_final = jacobi();
	for(int i=2;i<N;i++){
		double rel_jacobi = fabs((jacobi_final[i-2]-jacobi_init[i-2])/jacobi_init[i-2]);
		fprintf(of,"%.20e %.20e\n",t, rel_jacobi);
	}
	fclose(of);
	}
}

void problem_finish(){
	{
		FILE* of = fopen("energy.txt","w"); 
		double rel_energy = fabs((energy()-energy_init)/energy_init);
		double* am = angular_momentum();
		if (init_N!=N){
			rel_energy = -1;
		}
		fprintf(of,"%e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t%.20e\t",rel_energy,am[0],am[1],am[2],aminit[0],aminit[1],aminit[2]);
		fclose(of);
	}
	{
	FILE* of = fopen("jacobi.txt","w"); 
	double* jacobi_final = jacobi();
	for(int i=2;i<N;i++){
		double rel_jacobi = fabs((jacobi_final[i-2]-jacobi_init[i-2])/jacobi_init[i-2]);
		fprintf(of,"%.20e %.20e\n",t, rel_jacobi);
	}
	fclose(of);
	}
}

void input_binary_special(char* filename){
	FILE* inf = fopen(filename,"rb"); 
	long objects = 0;
	int _N;
	objects += fread(&_N,sizeof(int),1,inf);
	objects += fread(&tmax,sizeof(double),1,inf);
	objects += fread(&timescale,sizeof(double),1,inf);
	for (int i=0;i<_N;i++){
		struct particle p;
		objects += fread(&p,sizeof(struct particle),1,inf);
		particles_add(p);
	}
	fclose(inf);
}

double* angular_momentum(){
	mpf_t amx;
	mpf_init(amx);
	mpf_set_d(amx,0.);
	mpf_t amy;
	mpf_init(amy);
	mpf_set_d(amy,0.);
	mpf_t amz;
	mpf_init(amz);
	mpf_set_d(amz,0.);
	mpf_t mass;
	mpf_init(mass);
	mpf_t vx;
	mpf_init(vx);
	mpf_t vy;
	mpf_init(vy);
	mpf_t vz;
	mpf_init(vz);
	mpf_t x;
	mpf_init(x);
	mpf_t y;
	mpf_init(y);
	mpf_t z;
	mpf_init(z);
	for (int i=0;i<N;i++){
		mpf_t temp;
		mpf_init(temp);
		mpf_set_d(temp,particles[i].vx*particles[i].m);
		mpf_add(vx, vx, temp);
		mpf_set_d(temp,particles[i].vy*particles[i].m);
		mpf_add(vy, vy, temp);
		mpf_set_d(temp,particles[i].vz*particles[i].m);
		mpf_add(vz, vz, temp);
		
		mpf_set_d(temp,particles[i].x*particles[i].m);
		mpf_add(x, x, temp);
		mpf_set_d(temp,particles[i].y*particles[i].m);
		mpf_add(y, y, temp);
		mpf_set_d(temp,particles[i].z*particles[i].m);
		mpf_add(z, z, temp);


		mpf_set_d(temp,particles[i].m);
		mpf_add(mass, mass, temp);
	}
	mpf_div(vx,vx,mass);
	mpf_div(vy,vy,mass);
	mpf_div(vz,vz,mass);
	
	mpf_div(x,x,mass);
	mpf_div(y,y,mass);
	mpf_div(z,z,mass);

	for (int i=0;i<N;i++){
	
		double dvx = particles[i].vx-mpf_get_d(vx);
		double dvy = particles[i].vy-mpf_get_d(vy);
		double dvz = particles[i].vz-mpf_get_d(vz);
		
		double dx = particles[i].x-mpf_get_d(x);
		double dy = particles[i].y-mpf_get_d(y);
		double dz = particles[i].z-mpf_get_d(z);
		mpf_t temp;
		mpf_init(temp);
		
		mpf_set_d(temp,particles[i].m * (dy*dvz-dz*dvy));
		mpf_add(amx, amx,temp);
		mpf_set_d(temp,particles[i].m * (dz*dvx-dx*dvz));
		mpf_add(amy, amy,temp);
		mpf_set_d(temp,particles[i].m * (dx*dvy-dy*dvx));
		mpf_add(amz, amz,temp);
		
	}
	double* am = malloc(sizeof(double)*3);
	am[0] = mpf_get_d(amx);
	am[1] = mpf_get_d(amy);
	am[2] = mpf_get_d(amz);

	return am;
}
