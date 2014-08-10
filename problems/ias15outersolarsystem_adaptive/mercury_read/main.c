#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include "../../../src/particle.h"

const double k	 	= 0.01720209895;	// Gaussian constant 

double energy();
double* angular_momentum();
double energy_init;
double G;
long N=7;
struct particle* particles;


void readin(char* filename){
	// find central mass 
	system("grep mass ../param.in | awk '{print $5}' > tmp.txt");
	FILE* parf =  fopen("tmp.txt", "r");	
	fscanf(parf,"%lf",&(particles[0].m));
	fclose(parf);

	
	// Initial conditions
	particles[0].x = 0; particles[0].y = 0; particles[0].z = 0;
	particles[0].vx = 0; particles[0].vy = 0; particles[0].vz = 0;

	FILE* inpf = fopen(filename, "r");	
	for (int i=0;i<6;i++){
		char buffer[100];
		fgets(buffer, 100, inpf);
	}
	N=1;
	while(1){
		char buffer[1000];
		fgets(buffer, 1000, inpf);
		sscanf(&(buffer[8]),"%lf",&(particles[N].m));
		int d =fscanf(inpf,"%lf",&(particles[N].x));
		if(d!=1) break;

		fscanf(inpf,"%lf",&(particles[N].y));
		fscanf(inpf,"%lf",&(particles[N].z));
		fscanf(inpf,"%lf",&(particles[N].vx));
		fscanf(inpf,"%lf",&(particles[N].vy));
		fscanf(inpf,"%lf",&(particles[N].vz));
		fscanf(inpf,"%lf",&(particles[N].ax));
		fscanf(inpf,"%lf",&(particles[N].ay));
		fscanf(inpf,"%lf",&(particles[N].az));
		fgets(buffer, 1000, inpf);
		N++;
	}
	for (int i=0;i<N;i++){
	//	printf("%i particles mass = %e\n",i, particles[i].m);
	}
}

int main (int argc, char* argv[]){
	G		= k*k;
	particles = calloc(sizeof(struct particle),10000);

	mpf_set_default_prec(512);
	readin(argv[1]);
	int N_1 = N;
	energy_init = energy();
	double* aminit;
	aminit = angular_momentum();
	readin(argv[2]);
	int N_2 = N;
	if (N_1!=N_2){
		printf("%e\n",-1.);
		exit(0);
	}
		

	double rel_energy = fabs((energy()-energy_init)/energy_init);
	double* am;
	am = angular_momentum();
	printf("%e\t%e\t%e\t%e\t%e\t%e\t%e\n",rel_energy,am[0],am[1],am[2],aminit[0],aminit[1],aminit[2]);
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
			double r = sqrt(dx*dx + dy*dy + dz*dz);
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
