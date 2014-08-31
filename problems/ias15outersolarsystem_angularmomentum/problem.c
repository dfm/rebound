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
double energy_init;
double* angular_momentum();
double* aminit;

void problem_init(int argc, char* argv[]){
	// Setup constants
	dt 		= input_get_double(argc,argv,"dt",40);			// days
	N_active	= 5;
	tmax		= 4300e2;		// 1 Myr
	G		= k*k;
#ifdef INTEGRATOR_IAS15
	integrator_epsilon = input_get_double(argc,argv,"integrator_epsilon",0);
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
	aminit = angular_momentum();
	mpf_set_default_prec(512);
	energy_init = energy();
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

void problem_output(){
}

void problem_finish(){
	FILE* of = fopen("energy.txt","a+"); 
	double rel_energy = fabs((energy()-energy_init)/energy_init);
	double* amfini = angular_momentum();
	double _ai = sqrt(aminit[0]*aminit[0] + aminit[1]*aminit[1] + aminit[2]*aminit[2] );
	double _af = sqrt(amfini[0]*amfini[0] + amfini[1]*amfini[1] + amfini[2]*amfini[2] );
	fprintf(of,"%e\t%e\t%e\n",dt,rel_energy,fabs((_ai-_af)/_ai));
	fclose(of);

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
