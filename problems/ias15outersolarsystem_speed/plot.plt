#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 3in,2in
set xlabel "time to complete 100 orbits [s]"
set ylabel "relative energy error"
set logscale xy
set autoscale fix
set key top right
set yrange [1e-15:1]
set xrange [0.005:0.1]

set st d lp

plot \
"energy_wh.txt" t "WH (REBOUND)",  \
"energy_ias15.txt" t "IAS15"  lc rgb "dark-green", \
