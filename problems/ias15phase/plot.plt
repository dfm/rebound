#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 3in,2in
set xlabel "timestep [days]"
set ylabel "phase error after 100 orbits [rad]"
set logscale xy
set autoscale fix
set key top left
set format y '%.0e'

set st d lp

plot \
"energy_wh.txt" u 3:4  lc rgb "red"  t "    WH (REBOUND)",  \
"energy_mvs.txt" u 1:(abs($2-$3)) lc rgb "dark-red"  t "MVS (MERCURY)" , \
"energy_ias15.txt" u 3:4 t "IAS15"  lc rgb "dark-green", \
