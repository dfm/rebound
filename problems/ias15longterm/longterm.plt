#!/bin/gnuplot
set output "longterm.pdf"
set term pdf enhanced color dashed  size 5in,2in
set key bottom right
set logscale xy
set autoscale fix
set pointsize 0.4
set yrange [1e-16:1e-3]
set st d p

set label "dt = 1 day" at 5.8e6, 9.9e-9
set label "dt = 0.1 day" at 484286, 3.8e-7
set label "dt = 10 day" at 4.6e7, 4.1e-11
set label "dt = 100 day" at 1.1e8, 1.8e-9
set xrange [1:1e9]
set xtics 100
set xlabel "orbits"
set ylabel "relative energy error"
set format y "10^{%T}"
set format x "10^{%T}"
plot \
"<cat mercury*/energy.txt"  u ($2/4332.0):(abs($5)) lc rgb "dark-gray" pt 7 notit, 1/0 w p t "Mercury MVS" lc rgb "dark-gray" ps 1.5  pt 7 , \
"<cat ias15/energy_orbits_1.000e+10__integrator_epsilon_5.000e-09.txt" w p lt 1 pt 7 notit , 1/0 w p t "IAS15" lt 1 pt 7 ps 1.5, \
"<cat ias15_cs/energy_orbits_1.000e+10__integrator_epsilon_5.000e-09.txt" w p lt 2 pt 7 notit , 1/0 w p t "IAS15 cs" lt 2 pt 7 ps 1.5, \
"<cat ias15_cs_14/energy_orbits_1.000e+10__integrator_epsilon_5.000e-09.txt" w p lt 4 pt 7 notit , 1/0 w p t "IAS15 cs 14" lt 4 pt 7 ps 1.5, \
5e-16*sqrt(x) ls 3 lw 3 t "t^{0.5}", \





