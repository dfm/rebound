#!/bin/gnuplot
set output "timeseries.pdf" 
set terminal pdf color enhanced size 6in,5in
set logscale yy2
set autoscale fix
set multiplot layout 3,2
unset key
set ytics 1000
set st d p
set yrange [1e-16:1]

bottommargin = 0.08
keymargin = 0.1
topmargin = 0.05
ny = 3
px = 0
unset xtics

unset key;

do for [i=0:5]{

	if (i==2) {
		set key at screen 0.5,0.05 center horizontal 
	}else{
		unset key;
	}

	if (i%2==0){
		set rmargin 0
		set lmargin 12
		set format y "%g" 
		set format y "10^{%S}"
		if (i/2==(ny-1)/2){
			set ylabel "relative energy error"
		}else{
			unset ylabel
		}
		py = 0
	}else{
		set lmargin 0
		set rmargin 12
		unset ylabel
		set format y  ""
		py = 0.5
	}
	if (i/2==ny-1){
		px=keymargin
	}else{
		px=(1.-bottommargin-topmargin-keymargin)/ny*(ny-i/2-1)+bottommargin+keymargin

	}

	set origin py,px

	if (i/2==ny-1){
		set size 0.5,(1.-bottommargin-topmargin-keymargin)/ny+bottommargin
		set xlabel "time"
		set bmargin 5
		set format x "%g" 
	}else{
		set size 0.5,(1.-bottommargin-topmargin-keymargin)/ny
		set format x ""
		set bmargin 0
		set tmargin 0
	}

	titfile=system("sed '".(i+1)."q;d' titles.txt");
	set label 1 titfile at graph 0.91,0.1 right


	plot \
	"testcase_".i."/energy_timeseries_ias15.txt" t "IAS15,  {/Symbol e}=0.01", \
	"testcase_".i."/energy_timeseries_wh.txt" t "IAS15,  {/Symbol e}=0.01", \
	1e-14*x, \
	1e-14*sqrt(x)
}
