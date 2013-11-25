/**
 * @file 	problem.c
 * @brief 	Example problem: Close Encounter.
 * @author 	Hanno Rein <hanno@hanno-rein.de>
 * 
 * @section 	LICENSE
 * Copyright (c) 2013 Hanno Rein, Dave Spiegel
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
#include "main.h"
#include "tools.h"
#include "output.h"
#include "particle.h"
#include "integrator.h"

#ifdef OPENGL
extern int display_wire;
#endif // OPENGL

void problem_init(int argc, char* argv[]){
	dt 		= 0.1*2.*M_PI;			// initial timestep = 0.1 years
	integrator_epsilon = 1e-3;			// Accuracy
	tmax		= 5e2*2.*M_PI;			// 500 yrs

#ifdef OPENGL
	display_wire	= 1;			// Show orbits.
#endif // OPENGL
	init_boxwidth(10); 			// Init box with width 10 astronomical units

	struct particle star;
	star.m = 1;
	star.x = 0; 	star.y = 0; 	star.z = 0;
	star.vx = 0; 	star.vy = 0; 	star.vz = 0;
	particles_add(star);
	
	// Add planets between 1 and 2 AU
	int N_planets = 7;
	for (int i=0;i<N_planets;i++){
		double a = 1.+(double)i/(double)(N_planets-1);
		double v = sqrt(1./a); 	// Circular orbit
		struct particle planet;
		planet.m = 1e-3; 	// 1e-3 Solar Masses 
		planet.x = a; 	planet.y = 0; 	planet.z = 0;
		planet.vx = 0; 	planet.vy = v; 	planet.vz = 0;
		particles_add(planet); 
	}
	tools_move_to_center_of_momentum();
}

void problem_inloop(){
}

void problem_output(){
	if (output_check(10.*2.*M_PI)){  // Every 10 years
		output_timing();
	}
}

void problem_finish(){
}
