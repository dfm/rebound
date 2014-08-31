#!/bin/gnuplot
set output "longterm.pdf"
set term pdf enhanced color dashed  size 7in,3in
set key top left
set logscale xy
set autoscale fix
set pointsize 0.4
set yrange [1e-16:1e-3]
set st d p

set label "dt = 0.1 day" at 1242860, 10.8e-7
set label "dt = 1 day" at 1.2e7, 24.9e-9
set label "dt = 10 day" at 14.6e7, 1.1e-11
set label "dt = 100 day" at 1.1e8, 3.8e-9
set xrange [1:1e9]
set xtics 100
set xlabel "orbits"
set ylabel "relative energy error"
set format y "10^{%T}"
set format x "10^{%T}"
plot \
"<awk '{if (NR%5==0) print $0}' mercury_bs_*/energy.txt"  u ($2/4332.0):(abs($5)) lc rgb "purple" pt 7 notit, 1/0 w p lc rgb "purple" pt 7 ps 1.5 t "BS (MERCURY)", \
"<awk '{if (NR%2==0) print $0}' mercury_radau_1e-{04,06,08,16}/energy.txt"  u ($2/4332.0):(abs($5)) lc rgb "orange" pt 7 notit, 1/0 w p lc rgb "orange" pt 7 ps 1.5 t "    RADAU (MERCURY)", \
"<awk '{if (NR%2==0) print $0}' mercury_{0.1,1,10,100}/energy.txt"  u ($2/4332.0):(abs($5)) lc rgb "dark-red" pt 7 notit, 1/0 w p lc rgb "dark-red" pt 7 ps 1.5 t "   MVS (MERCURY)", \
"<awk '{if (NR%5==0) print $0}'  ias15_cs_rms/energy_orbits_1.000*" w p ps 0.5 lc rgb "light-green" pt 7 notit , 1/0 w p t "IAS15" lc rgb "light-green" pt 7 ps 1.5, \
"< sort -k 1 -g ias15_cs_rms/energy_orbits_1.000* | awk -f rms.awk" w p lc rgb "dark-green" pt 7 notit , 1/0 w p t "IAS15 RMS" lc rgb "dark-green" pt 7 ps 1.5, \
0.6e-16*sqrt(x) ls 0 lw 5 t "t^{0.5}", \





