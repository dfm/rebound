#!/bin/gnuplot
set output "plot.png" 
set terminal png enhanced
set logscale xyx2y2
set autoscale fix
set yrange [5e-17:]
set xlabel "wall time [s]"
set ylabel "relative energy error"


plot \
for [i=0:7] "<cat  testcase_".i."/energy_ias15_canonical.txt"  	u 1:($2+1e-16) notit ps 2 lt (i+1) w p, \
for [i=0:7] "<cat  testcase_".i."/energy_wh.txt"  		u 1:($2+1e-16) notit ps 1  lt (i+1) w l, \
for [i=0:7] "<cat  testcase_".i."/energy_mvs.txt"  		u 1:($2+1e-16) notit ps 1  lt (i+1) w lp, \
for [i=0:7] "<cat  testcase_".i."/energy_wh.txt"       		u (1.):(1.e-20) t "10^{".i."} orbits" ps 1  lt (i+1) w lp


