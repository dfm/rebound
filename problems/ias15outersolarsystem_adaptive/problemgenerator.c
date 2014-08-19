#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <getopt.h>
#include "../../src/particle.h"

int N = 0;
struct particle* particles;
const double G	= 0.01720209895*0.01720209895;
double tmax;
double timescale = 1;

int input_get_int(int argc, char** argv, const char* argument, int _default);
void output_binary(char* filename);

int main(int argc, char* argv[]){
	particles = calloc(sizeof(struct particle),1000); // more then needed
	int testcase = input_get_int(argc,argv,"testcase",10000);
	
	switch (testcase){
		case 0: // Kozai  89.8
		{
			double semia = 1.0;

			struct particle star; 
			star.m  = 1;
			star.x  = 0; star.y  = 0; star.z  = 0; 
			star.vx = 0; star.vy = 0; star.vz = 0;
			particles[N++] = star;

			
			// The planet 
			struct particle planet; 
			planet.m  = 0.01*star.m;
			planet.x  = semia; 
			planet.y  = 0; 
			planet.z  = 0; 
			planet.vx = 0; 
			planet.vy = sqrt(G*star.m/semia); 
			planet.vz = 0;
			particles[N++] = planet;
			
			// The perturber
			struct particle perturber; 
			perturber.m  = 0.1*star.m;
			perturber.x  = semia*10.; 
			perturber.y  = 0; 
			perturber.z  = 0; 
			double inc_perturber = 88; // 80.
			perturber.vx = 0; 
			perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			particles[N++] = perturber;
			tmax	= 1e4*365.*sqrt(semia*semia*semia/star.m);
			break;
		}
		case 1: // Kozai  80
		{
			double semia = 1.0;

			struct particle star; 
			star.m  = 1;
			star.x  = 0; star.y  = 0; star.z  = 0; 
			star.vx = 0; star.vy = 0; star.vz = 0;
			particles[N++] = star;

			
			// The planet 
			struct particle planet; 
			planet.m  = 0.01*star.m;
			planet.x  = semia; 
			planet.y  = 0; 
			planet.z  = 0; 
			planet.vx = 0; 
			planet.vy = sqrt(G*star.m/semia); 
			planet.vz = 0;
			particles[N++] = planet;
			
			// The perturber
			struct particle perturber; 
			perturber.m  = 0.1*star.m;
			perturber.x  = semia*10.; 
			perturber.y  = 0; 
			perturber.z  = 0; 
			double inc_perturber =  80.;
			perturber.vx = 0; 
			perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			particles[N++] = perturber;
			tmax	= 1e4*365.*sqrt(semia*semia*semia/star.m);
			break;
		}
		case 2: // Kozai  89.8, size = 0.01
		{
			double semia = 0.01;

			struct particle star; 
			star.m  = 1;
			star.x  = 0; star.y  = 0; star.z  = 0; 
			star.vx = 0; star.vy = 0; star.vz = 0;
			particles[N++] = star;

			
			// The planet 
			struct particle planet; 
			planet.m  = 0.01*star.m;
			planet.x  = semia; 
			planet.y  = 0; 
			planet.z  = 0; 
			planet.vx = 0; 
			planet.vy = sqrt(G*star.m/semia); 
			planet.vz = 0;
			particles[N++] = planet;
			
			// The perturber
			struct particle perturber; 
			perturber.m  = 0.1*star.m;
			perturber.x  = semia*10.; 
			perturber.y  = 0; 
			perturber.z  = 0; 
			double inc_perturber = 89.9; // 80.
			perturber.vx = 0; 
			perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			particles[N++] = perturber;
			tmax	= 1e4*365.*sqrt(semia*semia*semia/star.m);
			timescale  = sqrt(semia*semia*semia/star.m);
			break;
		}
		case 3: // Kozai  80 size = 0.01
		{
			double semia = 0.01;

			struct particle star; 
			star.m  = 1;
			star.x  = 0; star.y  = 0; star.z  = 0; 
			star.vx = 0; star.vy = 0; star.vz = 0;
			particles[N++] = star;

			
			// The planet 
			struct particle planet; 
			planet.m  = 0.01*star.m;
			planet.x  = semia; 
			planet.y  = 0; 
			planet.z  = 0; 
			planet.vx = 0; 
			planet.vy = sqrt(G*star.m/semia); 
			planet.vz = 0;
			particles[N++] = planet;
			
			// The perturber
			struct particle perturber; 
			perturber.m  = 0.1*star.m;
			perturber.x  = semia*10.; 
			perturber.y  = 0; 
			perturber.z  = 0; 
			double inc_perturber =  80.;
			perturber.vx = 0; 
			perturber.vy = cos(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			perturber.vz = sin(inc_perturber/180.*M_PI)*sqrt(G*(star.m+perturber.m)/perturber.x); 
			particles[N++] = perturber;
			tmax	= 1e4*365.*sqrt(semia*semia*semia/star.m);
			timescale  = sqrt(semia*semia*semia/star.m);
			break;
		}
		case 4: // outer solar system
		{
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
			// Initial conditions
			for (int i=0;i<6;i++){
				struct particle p;
				p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
				p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
				p.ax = 0; 			p.ay = 0; 			p.az = 0;
				p.m  = ss_mass[i];
				particles[N++] = p; 
			}
			for (int i=1;i<N;i++){
				particles[i].x -= particles[0].x;	particles[i].y -= particles[0].y;	particles[i].z -= particles[0].z;
				particles[i].vx -= particles[0].vx;	particles[i].vy -= particles[0].vy;	particles[i].vz -= particles[0].vz;
			}
			particles[0].x = 0;	particles[0].y = 0;	particles[0].z = 0;
			particles[0].vx= 0;	particles[0].vy= 0;	particles[0].vz= 0;
			tmax		= 365e4;		// 10000 yr
			break;
		}
		case 5: // outer solar system, scaled 0.01
		{
			double scale = 0.01;
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
			// Initial conditions
			for (int i=0;i<6;i++){
				struct particle p;
				p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
				p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
				p.ax = 0; 			p.ay = 0; 			p.az = 0;
				p.m  = ss_mass[i];
				particles[N++] = p; 
			}
			for (int i=1;i<N;i++){
				particles[i].x -= particles[0].x;	particles[i].y -= particles[0].y;	particles[i].z -= particles[0].z;
				particles[i].vx -= particles[0].vx;	particles[i].vy -= particles[0].vy;	particles[i].vz -= particles[0].vz;
				particles[i].x *=scale;
				particles[i].y *=scale;
				particles[i].z *=scale;
				particles[i].vx *=sqrt(1./scale);
				particles[i].vy *=sqrt(1./scale);
				particles[i].vz *=sqrt(1./scale);
			}
			particles[0].x = 0;	particles[0].y = 0;	particles[0].z = 0;
			particles[0].vx= 0;	particles[0].vy= 0;	particles[0].vz= 0;
			tmax		= 365e4*sqrt(scale*scale*scale);		// 10000 yr
			timescale  = sqrt(scale*scale*scale);
			break;
		}
		case 6: // Eccentric orbit 0.99
		{
			double semia = 1.0;

			struct particle star; 
			star.m  = 1;
			star.x  = 0; star.y  = 0; star.z  = 0; 
			star.vx = 0; star.vy = 0; star.vz = 0;
			particles[N++] = star;

			
			// The planet 
			struct particle planet; 
			double e = 0.99;
			planet.m  = 0.01*star.m;
			planet.x  = semia*(1.+e); 
			planet.y  = 0; 
			planet.z  = 0; 
			planet.vx = 0; 
			planet.vy = sqrt((1.-e)/(1.+e)*G*star.m/semia); 
			planet.vz = 0;
			particles[N++] = planet;
			
			tmax	= 1e3*365.*sqrt(semia*semia*semia/star.m);
			break;
		}
		case 7: // Eccentric orbit 0.99999
		{
			double semia = 1.0;

			struct particle star; 
			star.m  = 1;
			star.x  = 0; star.y  = 0; star.z  = 0; 
			star.vx = 0; star.vy = 0; star.vz = 0;
			particles[N++] = star;

			
			// The planet 
			struct particle planet; 
			double e = 0.99999;
			planet.m  = 0.01*star.m;
			planet.x  = semia*(1.+e); 
			planet.y  = 0; 
			planet.z  = 0; 
			planet.vx = 0; 
			planet.vy = sqrt((1.-e)/(1.+e)*G*star.m/semia); 
			planet.vz = 0;
			particles[N++] = planet;
			
			tmax	= 1e3*365.*sqrt(semia*semia*semia/star.m);
			break;
		}
		case 8: // outer solar system, scaled 0.01, long term
		{
			double scale = 0.01;
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
			// Initial conditions
			for (int i=0;i<6;i++){
				struct particle p;
				p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
				p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
				p.ax = 0; 			p.ay = 0; 			p.az = 0;
				p.m  = ss_mass[i];
				particles[N++] = p; 
			}
			for (int i=1;i<N;i++){
				particles[i].x -= particles[0].x;	particles[i].y -= particles[0].y;	particles[i].z -= particles[0].z;
				particles[i].vx -= particles[0].vx;	particles[i].vy -= particles[0].vy;	particles[i].vz -= particles[0].vz;
				particles[i].x *=scale;
				particles[i].y *=scale;
				particles[i].z *=scale;
				particles[i].vx *=sqrt(1./scale);
				particles[i].vy *=sqrt(1./scale);
				particles[i].vz *=sqrt(1./scale);
			}
			particles[0].x = 0;	particles[0].y = 0;	particles[0].z = 0;
			particles[0].vx= 0;	particles[0].vy= 0;	particles[0].vz= 0;
			tmax		= 365e7*sqrt(scale*scale*scale);		// 10.000.000 yr
			timescale  = sqrt(scale*scale*scale);
			break;
		}
		case 9: // kozai test getting stuck
		{
			double ss_pos[3][3] = 
			{
				{-7.078841e-01,   -2.615421e-03,   -5.543430e-01},
				{-7.080543e-01,   -2.698939e-03,   -5.543991e-01},
				{ 7.149646e+00,    2.642411e-02,    5.598870e+00 },
			};
			double ss_vel[3][3] = 
			{
				{ -4.045634e-03,   1.072532e-02,    -1.313862e-02},
				{ 4.367986e-01,    -1.072698e+00,   1.273079e+00},
				{ -3.223527e-03,   1.654099e-05,    4.078342e-03},
			};

			double ss_mass[3] =
			{
				1.,0.01,0.1
			};
			// Initial conditions
			for (int i=0;i<3;i++){
				struct particle p;
				p.x  = ss_pos[i][0]; 		p.y  = ss_pos[i][1];	 	p.z  = ss_pos[i][2];
				p.vx = ss_vel[i][0]; 		p.vy = ss_vel[i][1];	 	p.vz = ss_vel[i][2];
				p.ax = 0; 			p.ay = 0; 			p.az = 0;
				p.m  = ss_mass[i];
				particles[N++] = p; 
			}
			tmax		= 365000;		// 10.000.000 yr
			timescale  = 0.0000001;
			break;
		}
		default:
			printf("test case not found\n");
			exit(-1);
			break;
	}


	output_binary("particles.bin");

	FILE* of;
	
	// create files for mercury
	of  = fopen("big.in","w"); 
	fprintf(of,")O+_06 Big-body initial data  (WARNING: Do not delete this line!!)\n");
	fprintf(of,") Lines beginning with `)' are ignored.\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of," style (Cartesian, Asteroidal, Cometary) = Cartesian\n");
	fprintf(of," epoch (in days) = 0.0\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	for (int i=1;i<N;i++){
		fprintf(of,"NAME%i m=%.16e \n",i,particles[i].m);
		fprintf(of," %.16e %.16e %.16e\n",particles[i].x,particles[i].y,particles[i].z);
		fprintf(of," %.16e %.16e %.16e\n",particles[i].vx,particles[i].vy,particles[i].vz);
		fprintf(of," 0. 0. 0.\n");
	}
	fclose(of);
	
	of = fopen("timescale.txt","w"); 
	fprintf(of,"%.16f\n",timescale);
	fclose(of);

	
	of = fopen("param.in","w"); 
	fprintf(of,")O+_06 Integration parameters  (WARNING: Do not delete this line!!)\n");
	fprintf(of,") Lines beginning with `)' are ignored.\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of,") Important integration parameters:\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of," algorithm (MVS, BS, BS2, RADAU, HYBRID etc) = SCHEME\n");
	fprintf(of," start time (days)=  0d0\n");
	fprintf(of," stop time (days) = %.16e\n",tmax);
	fprintf(of," output interval (days) = %.16e\n",tmax);
	fprintf(of," timestep (days) = STEPSIZE\n");
	fprintf(of," accuracy parameter = EPSILON\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of,") Integration options:\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of," stop integration after a close encounter = no\n");
	fprintf(of," allow collisions to occur = no\n");
	fprintf(of," include collisional fragmentation = no\n");
	fprintf(of," express time in days or years = days\n");
	fprintf(of," express time relative to integration start time = no\n");
	fprintf(of," output precision = high\n");
	fprintf(of," < not used at present >\n");
	fprintf(of," include relativity in integration= no\n");
	fprintf(of," include user-defined force = no\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of,") These parameters do not need to be adjusted often:\n");
	fprintf(of,")---------------------------------------------------------------------\n");
	fprintf(of," ejection distance (AU)= 1000\n");
	fprintf(of," radius of central body (AU) = 0\n");
	fprintf(of," central mass (solar) = %.16e\n",particles[0].m);
	fprintf(of," central J2 = 0\n");
	fprintf(of," central J4 = 0\n");
	fprintf(of," central J6 = 0\n");
	fprintf(of," < not used at present >\n");
	fprintf(of," < not used at present >\n");
	fprintf(of," Hybrid integrator changeover (Hill radii) = 3.\n");
	fprintf(of," number of timesteps between data dumps = 50000000\n");
	fprintf(of," number of timesteps between periodic effects = 10000000\n");
	fclose(of);

}


char* input_get_argument(int argc, char** argv, const char* argument){
	opterr = 0;
	optind = 1;
  	while (1) {
      		struct option long_options[] = {
	  		{NULL, required_argument, 0, 'a'},
			{0,0,0,0}
		};

		long_options[0].name = argument;

      		/* getopt_long stores the option index here.   */
      		int option_index = 0;
		//				short options. format abc:d::
      		int c = getopt_long (argc, argv, "", long_options, &option_index);

      		/* Detect the end of the options.   */
      		if (c == -1) break;

      		switch (c)
		{
			case 'a':
				return optarg;
				break;
			default:
				break;
		}
  	}
	return NULL;
}
int input_get_int(int argc, char** argv, const char* argument, int _default){
	char* value = input_get_argument(argc,argv,argument);
	if (value){
		return atoi(value);
	}
	return _default;
}

void output_binary(char* filename){
	FILE* of = fopen(filename,"wb"); 
	fwrite(&N,sizeof(int),1,of);
	fwrite(&tmax,sizeof(double),1,of);
	fwrite(&timescale,sizeof(double),1,of);
	for (int i=0;i<N;i++){
		struct particle p = particles[i];
		fwrite(&(p),sizeof(struct particle),1,of);
	}
	fclose(of);
}
