from __future__ import division

cimport cython

from libc.math cimport sqrt, M_PI


cdef extern from "boundaries.h":
    cdef void boundaries_check()

cdef extern from "gravity.h":
    cdef void gravity_calculate_acceleration()

cdef extern from "problem.h":
    cdef void problem_inloop()

cdef extern from "tools.h":
    cdef void tools_move_to_center_of_momentum()

cdef extern from "integrator.h":
    cdef void integrator_part1()
    cdef void integrator_part2()

cdef extern from "particle.h":
    struct particle:
        double x
        double y
        double z
        double vx
        double vy
        double vz
        double ax
        double ay
        double az
        double m
    cdef void particles_add(particle pt)

cdef extern from "main.h":
    cdef double G
    cdef double boxsize
    cdef double dt
    cdef double tmax
    cdef int run_main(int argc, char* argv[])
    cdef void init_box()


def run_sim():
    cdef int argc = 0
    cdef char** argv
    run_main(argc, argv)


cdef public void problem_init(int argc, char** argv):
    global G, boxsize, dt, tmax

    G = 1
    cdef double e_testparticle = 1.-1e-7
    cdef double mass_scale = 1.
    cdef double size_scale = 1

    boxsize = 25.*size_scale
    init_box()

    cdef particle star
    star.m  = mass_scale
    star.x  = 0
    star.y  = 0
    star.z  = 0
    star.vx = 0
    star.vy = 0
    star.vz = 0
    particles_add(star)

    cdef  particle planet
    planet.m  = 0;
    planet.x  = size_scale*(1.-e_testparticle)
    planet.y  = 0
    planet.z  = 0
    planet.vx = 0
    planet.vy = sqrt((1.+e_testparticle)/(1.-e_testparticle)*mass_scale/size_scale)
    planet.vz = 0
    particles_add(planet)

    tools_move_to_center_of_momentum()

    dt = 1e-13*sqrt(size_scale*size_scale*size_scale/mass_scale)
    tmax = 1e2*2.*M_PI*sqrt(size_scale*size_scale*size_scale/mass_scale)

cdef public void problem_inloop():
    pass

cdef public void problem_output():
    pass

cdef public void problem_finish():
    pass

# if(output_check(tmax/10000.)){		// outputs to the screen
#     output_timing();
#     }
# }

# void problem_output(){
#         // Output the time and the current timestep. Plot it to see how IAS15 automatically reduces the timestep at pericentre.
#         FILE* of = fopen("timestep.txt","a");
#         fprintf(of,"%e\t%e\t\n",t/tmax,dt/tmax);
#         fclose(of);
#         }

# void problem_finish(){
#        }
