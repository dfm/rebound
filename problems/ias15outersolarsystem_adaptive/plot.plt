#!/bin/gnuplot
set output "plot.pdf" 
set terminal pdf color enhanced size 6in,7in
set logscale xyx2y2
set autoscale fix
set yrange [1e-16:1]
set xrange [0.07:15]
set multiplot 
unset key
set ytics 1000
set st d p

bottommargin = 0.08
keymargin = 0.08
topmargin = 0.05
ny = 4
px = 0

do for [i=0:7]{
	if (i==8) {
		set xrange [7:1200]
	}
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
			set ylabel " "
		}
		py = 0
	}else{
		set lmargin 0
		set rmargin 12
		unset ylabel
		set format y  ""
		py = 0.5
	}
	if (i/2==ny-1 || (i-1)/2==ny-1){
		px=keymargin
	}else{
		px=(1.-bottommargin-topmargin-keymargin)/ny*(ny-i/2-1)+bottommargin+keymargin

	}

	set origin py,px

	if (i/2==ny-1 || (i-1)/2==ny-1){
		set size 0.5,(1.-bottommargin-topmargin-keymargin)/ny+bottommargin
		set xlabel "time to complete run [s]"
		set bmargin 6
		set format x "%g" 
	}else{
		set size 0.5,(1.-bottommargin-topmargin-keymargin)/ny
		set format x ""
		set bmargin 0
		set tmargin 0
	}

	titfile=system("sed '".(i+1)."q;d' titles.txt");
	set label 1 titfile at graph 0.01,0.1 left

	if (i==7){
		set rmargin 14.95
		set bmargin 7.57
	}

	plot \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_ias15.txt" 		lt 1 t "IAS15", \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_ias15_canonical.txt" u (0.0001):(1.) t "IAS15, default" ps 2 lt 6, \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_ra15.txt" 		lt 3 t "RA15 (REBOUND)", \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_wh.txt" 		lt 4 t "      WH (REBOUND)",  \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_bs2.txt" 		lt 5 t "BS2 (MERCURY)",  \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_radau.txt" 		lt 2 t "  RADAU (MERCURY)" ,  \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_mvs.txt" 		lt 7 t "      MWS (MERCURY)",  \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_ias15_canonical.txt" notit ps 4 lt 6, \
	"<awk -v type=regular -f regular.awk testcase_".i."/energy_ias15_canonical.txt" notit lt 1, \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_ias15.txt" 		u 1:2:(0):(0.5) notit with vectors head  lt 1, \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_ra15.txt" 		u 1:2:(0):(0.5) notit with vectors head  lt 3, \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_wh.txt" 		u 1:2:(0):(0.5) notit with vectors head  lt 4,  \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_bs2.txt" 		u 1:2:(0):(0.5) notit with vectors head  lt 5,  \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_radau.txt" 		u 1:2:(0):(0.5) notit with vectors head  lt 2,  \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_mvs.txt" 		u 1:2:(0):(0.5) notit with vectors head  lt 7,  \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_ias15.txt" 		notit lt 1 ps 0.5 , \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_ra15.txt" 		notit lt 3 ps 0.5 , \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_wh.txt" 		notit lt 4 ps 0.5 ,  \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_bs2.txt" 		notit lt 5 ps 0.5 ,  \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_radau.txt" 		notit lt 2 ps 0.5 ,  \
	"<awk -v type=limit   -f regular.awk testcase_".i."/energy_mvs.txt" 		notit lt 7 ps 0.5 ,  \
}
