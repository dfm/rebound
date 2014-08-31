#!/bin/bash
function rundt {
	points=100
	min=1
	max=4
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=10; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=10; e($exp*l(10))"  | bc -l )
		./nbody --dt=$e 
	done
}
function rundtmercury {
	echo "Running mercury dt $1"
	points=20
	min=1
	max=4
	for i in $(seq 0 $points)
	do 
		exp=$(echo "scale=16; ($max-($min))/$points*$i+($min) " |bc)
		e=$(echo "scale=16; e($exp*l(10))"  | bc -l )
		rm -f *.tmp
		rm -f *.dmp
		rm -f *.out
		cp big.template big.in
		sed "s/EPSILON/1/g" param.template > param.template2
		sed "s/SCHEME/$1/g" param.template2 > param.template3
		sed "s/STEPSIZE/$e/g" param.template3 > param.in
		utime="$( TIMEFORMAT='%R';time ( doalarm $2 ./mercury6 ) 2>&1 1>data.tmp )"
		energy="$(cat data.tmp | tail -n 2 | head -n 1 | awk '{print $5}')"
		am="$(cat data.tmp | tail -n 2 | head -n 1 | awk '{print $7}')"
		echo "$e $energy $am" >> energy.txt 
		echo "$e $energy $am"

		rm -f param.template?

	done
}

rm energy_*.txt
make wh
rundt
mv energy.txt energy_wh.txt

make ias15
rundt
mv energy.txt energy_ias15.txt

pushd mercury 
rundtmercury mvs 10
mv energy.txt ../energy_mvs.txt
popd
