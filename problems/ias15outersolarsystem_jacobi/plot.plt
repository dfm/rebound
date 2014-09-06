#!/bin/gnuplot
set terminal pdf monochrome size 2.8in,2.8in
set xrange [-6*5.2:6*5.2]
set yrange [-6*5.2:6*5.2]
set xlabel "x [AU]"
set ylabel "y [AU]"
set size square

set output "jacobi_1.pdf" 
plot "orbits_init.txt" u ($1*5.2):($2*5.2) w l lc rgb "#555555" notit, \
"orbits_jupiter.txt" u ($1*5.2):($2*5.2) w l lw 5 notit 

set output "jacobi_2.pdf" 
plot "orbits_final.txt" u ($1*5.2):($2*5.2) w l lc rgb "#555555" notit, \
"orbits_jupiter.txt" u ($1*5.2):($2*5.2) w l lw 5 notit
