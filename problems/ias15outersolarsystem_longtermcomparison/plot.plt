#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf size 9in,2in
set logscale xyx2y2
set autoscale fix
set yrange [5e-17:1]
set xrange [1e-3:1e4]
set xlabel "wall time [s]"
set ylabel "relative energy error"
set multiplot layout 1,3

set format y "%g" 
set format y "10^{%S}"
set ytics 1000


#do for [i=0:5] {
do for [i in "1 3 5"] {

set title "".(10**i)." orbits"


plot \
"<cat  testcase_".i."/energy_ias15_canonical.txt"  	u 1:($2+1e-16) pt 7 ps 2 lt 1 lc rgb "dark-green" t "IAS15" w p, \
"<cat  testcase_".i."/energy_wh.txt"  			u 1:($2+1e-16) ps 1 lt 1 lc rgb "dark-red" t "WH (REBOUND)" w lp,  \
"<cat  testcase_".i."/energy_mvs.txt"  			u 1:($2+1e-16) ps 1 lt 2 lc rgb "dark-blue" t "MVS (MERCURY)" w lp, \



}
