#!/bin/gnuplot
set output "kozaionly.pdf"
set term pdf enhanced   size 7in,2.0in
set logscale xyy2
set autoscale fix
set pointsize 1
set yrange [1e-16:1]
set xrange [0.9e-2:2e1]
set ytics 1000
set format y '%.0e'

set multiplot layout 1,2

set key center left
set xlabel "cpu time to complete one Kozai cycle [s]"
set ylabel "relative energy error"
#set format y "{10^{%T}}"
#set format y2 "{10^{%T}}"
#set format x "{10^{%T}}"
plot \
"<cat testcase_0/energy_ias15*" w p lc rgb "dark-green" pt 7 notit , 1/0 w p t "IAS15" lc rgb "dark-green" pt 7 ps 1.5, \
"testcase_0/energy_ias15_canonical.txt" w p lc rgb "dark-green" pt 6 ps 3 notit , 1/0 w p t "IAS15 (canonical)" lc rgb "dark-green" pt 6 ps 1.5, \
"testcase_0/energy_bs.txt"  lc rgb "purple" pt 7 notit, 1/0 w p  lc rgb "dark-gray" ps 1  pt 7 notit , 1/0 w p lc rgb "purple" pt 7 ps 1.5 t "BS (MERCURY)", \
"testcase_0/energy_wh.txt" w p lc rgb "red" pt 1 notit , 1/0 w p t "WH (REBOUND)" lc rgb "red" pt 1 ps 1.5, \
"testcase_0/energy_mvs.txt" w p lc rgb "dark-red" pt 2 notit , 1/0 w p t "MVS (MERCURY)" lc rgb "dark-red" pt 2 ps 1.5, \



set key top left
set ylabel "relative angular momentum error"
plot \
	"<cat testcase_0/energy_ias15*"		u 1:(abs((sqrt($3*$3+$4*$4+$5*$5)-sqrt($6*$6+$7*$7+$8*$8))/sqrt($6*$6+$7*$7+$8*$8))+1e-16) w p lc rgb "dark-green" pt 7 notit , 1/0 w p t "IAS15" lc rgb "dark-green" pt 7 ps 1.5, \
	"testcase_0/energy_ias15_canonical.txt" u 1:(abs((sqrt($3*$3+$4*$4+$5*$5)-sqrt($6*$6+$7*$7+$8*$8))/sqrt($6*$6+$7*$7+$8*$8))+1e-16) w p lc rgb "dark-green" pt 6 ps 3 notit , 1/0 w p t "IAS15 (canonical)" lc rgb "dark-green" pt 6 ps 1.5, \
	"testcase_0/energy_bs.txt" 		u 1:(abs((sqrt($3*$3+$4*$4+$5*$5)-sqrt($6*$6+$7*$7+$8*$8))/sqrt($6*$6+$7*$7+$8*$8))+1e-16) w p lc rgb "purple" pt 7 notit, 1/0  lc rgb "dark-gray" ps 1  pt 7 notit , 1/0 w p lc rgb "purple" pt 7 ps 1.5 t "BS (MERCURY)", \
	"testcase_0/energy_wh.txt" 		u 1:(abs((sqrt($3*$3+$4*$4+$5*$5)-sqrt($6*$6+$7*$7+$8*$8))/sqrt($6*$6+$7*$7+$8*$8))+1e-16) w p lc rgb "red" pt 1 notit , 1/0 w p t "WH (REBOUND)" lc rgb "red" pt 1 ps 1.5, \
	"testcase_0/energy_mvs.txt" 		u 1:(abs((sqrt($3*$3+$4*$4+$5*$5)-sqrt($6*$6+$7*$7+$8*$8))/sqrt($6*$6+$7*$7+$8*$8))+1e-16) w p lc rgb "dark-red" pt 2 notit , 1/0 w p t "MVS (MERCURY)" lc rgb "dark-red" pt 2 ps 1.5, \
