#!/bin/bash
function runepsilon {
	echo "Running mercury epsilon $1"
	rm -f ../energy_$1.txt
	points=20
	min=-14
	max=-5
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=16; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=16; e($exp*l(10))"  | bc -l )
		rm -f *.tmp
		rm -f *.dmp
		rm -f *.out
		sed "s/EPSILON/$e/g" param.template > param.template2
		sed "s/SCHEME/$1/g" param.template2 > param.template3
		sed "s/STEPSIZE/1/g" param.template3 > param.in
		utime="$( TIMEFORMAT='%R';time ( ./mercury6 ) 2>&1 1>/dev/null )"
		energy="$(../mercury_read/mercury_energy big.in big.dmp)"
		echo "$utime $energy $e" >> ../energy_$1.txt 
		echo "$utime $energy $e"  
		rm -f param.template?

	done
}
function rundt {
	echo "Running mercury dt $1 /  $2 /  $3"
	rm -f ../energy_$1.txt
	points=20
	min=$2
	max=$3
	for i in $(seq 0 $points)
	do 
		timescale=$(cat ../timescale.txt)
		exp=$(echo "scale=16; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=16; e($exp*l(10))*$timescale"  | bc -l )
		rm -f *.tmp
		rm -f *.dmp
		rm -f *.out
		sed "s/EPSILON/1/g" param.template > param.template2
		sed "s/SCHEME/$1/g" param.template2 > param.template3
		sed "s/STEPSIZE/$e/g" param.template3 > param.in
		utime="$( TIMEFORMAT='%R';time ( ./mercury6 ) 2>&1 1>/dev/null )"
		energy="$(../mercury_read/mercury_energy big.in big.dmp)"
		echo "$utime $energy $e" >> ../energy_$1.txt 
		echo "$utime $energy $e"  
		rm -f param.template?

	done
}
function runepsilonnbody {
	echo "Running REBOUND epsilon $1"
	rm -f energy_$1.txt
	points=10
	min=$2
	max=$3
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=16; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=16; e($exp*l(10))"  | bc -l )
		utime="$( TIMEFORMAT='%R';time ( ./nbody --integrator_epsilon=$e 2>&1 ) 2>&1 1>/dev/null )"
		if [ ! -f energy.txt ]; then
			energy="-1"
		else
			energy="$(cat energy.txt)"
		fi
		echo "Did not finish in time."
		echo "$4 1. $e" >> energy_$1.txt 
		rm -f energy.txt

	done
}
function runnbodycanonical {
	echo "Running REBOUND canonical $1"
	dt="10."
	utime="$( TIMEFORMAT='%R';time ( ./nbody --dt=$dt  --outputenergy=$2 2>&1 ) 2>&1 1>/dev/null )"
	if [ ! -f energy.txt ]; then
		energy="-1"
	else
		energy="$(cat energy.txt)"
	fi
	if [ $2 -eq 1 ]; then 
		mv energy_timeseries.txt energy_timeseries_$1.txt
	else
		echo "$utime $energy " >> energy_$1_canonical.txt 
		echo "$utime $energy "  
	fi
	rm -f energy.txt

}

function rundtnbody {
	echo "Running REBOUND dt $1"
	rm -f energy_$1.txt
	points=20
	min=$2
	max=$3
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=16; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=16; e($exp*l(10))"  | bc -l )
		utime="$( TIMEFORMAT='%R';time ( ./nbody --dt=$e 2>&1 ) 2>&1 1>/dev/null )"
		if [ ! -f energy.txt ]; then
			energy="-1"
		else
			energy="$(cat energy.txt)"
		fi
		echo "$utime $energy $e" >> energy_$1.txt 
		echo "$utime $energy $e"  
		rm -f energy.txt

	done
}

make problemgenerator
rm -rf energy_*.txt


for t in $(seq 5 5)
do
	echo "###################################"
	echo "Running test case $t"


	./problemgenerator --testcase=$t

	make -s ias15
	#runepsilonnbody ias15 -11 -5 $runtime
        runnbodycanonical ias15 0 
	#runnbodycanonical ias15 1 $runtime
     
	#make -s ra15
	#runepsilonnbody ra15 -14 -6 $runtime

	make -s wh
	mindt=$(echo "scale=16; l(e($t*l(10))/1000000.)/l(10) " |bc -l)
	rundtnbody wh $mindt 4 

	pushd mercury
	rm -f *.tmp
	rm -f *.dmp
	rm -f *.out
	rm -f output.txt
	#runepsilon bs $runtime  
	#runepsilon bs2 $runtime
	#runepsilon radau $runtime
	#runepsilon hybrid 
	rundt mvs $mindt 4
	popd

	rm -rf testcase_$t
	mkdir testcase_$t
	mv energy*.txt testcase_$t/

done
exit

