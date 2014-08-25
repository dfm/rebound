#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 3in,2in
set xlabel "timestep [days]"
set ylabel "relative energy after 100 orbits"
set logscale xy
set autoscale fix
set yrange [1e-16:0.1]
set xrange[:1.2e4]
set key right bottom

set st d lp

plot \
1e-3*(x/1000.)**2 lt 4 t "dt^{2}", \
1e-12*(x/1000.)**15 lt 5  t "dt^{15}", \
"energy_wh.txt" lt 1 t "WH (REBOUND)",  \
"energy_mvs.txt"  lt 2  lc rgb "blue"  t "MVS (MERCURY)" , \
"energy_ias15.txt" lt 3 lc rgb "dark-green" t "IAS15"  , \
