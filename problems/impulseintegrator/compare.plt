#!/usr/bin/gnuplot
set output "compare.png"
#set output "compare.pdf"
set key top left
set terminal png enhanced 
set xlabel "timestep"
set title "v=100"
set ylabel "relative energy error"
#set terminal pdf enhanced dashed
set logscale xy
set yrange []
plot \
"error_radau15.txt" 		u 1:2 w lp lc 1 lt 1 t "radau15",\
"error_leapfrog.txt" 		u 1:2 w lp lc 2 lt 1 t "leapfrog", \
"error_euler.txt" 		u 1:2 w lp lc 7 lt 1 t "euler", \
"error_modifiedeuler.txt" 	u 1:2 w lp lc 8 lt 1 t "modified euler", \
"error_impulse.txt" 		u 1:2 w lp lc 3 lt 1 t "impulse (1st order)", \
"error_impulse2ndorder.txt" 	u 1:2 w lp lc 4 lt 1 t "impulse (2nd order)", \
1*x 				lc 5 lt 3 t "x", \
100*x**2 			lc 6 lt 3 t "x^2"

#plot "error.txt" u 1:2 w lp, "" u 1:3 w lp, "" u 1:4 w lp, "" u 1:5 w lp, x t "~x", x**2 t "~x^2"

