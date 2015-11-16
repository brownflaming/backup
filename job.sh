#!/bin/bash

cd bin

bendersFlag=1
imprvdBendersFlag=0

for T in $(seq 5 5)
do
	echo "Time horizon: $T; Benders' cut: $bendersFlag; Improved Benders' cut: $imprvdBendersFlag" >> table.txt
	cp data_$T/* data/
	for n in 5
	do
		cd data
		echo $n > numFWsample.dat
		cd ..
		./gep $bendersFlag $imprvdBendersFlag
	done
done
