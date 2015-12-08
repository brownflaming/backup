#!/bin/bash

cd bin

bendersFlag=1
imprvdBendersFlag=1

for T in $(seq 5 10)
do
	echo "Time horizon: $T; Benders' cut: $bendersFlag; Improved Benders' cut: $imprvdBendersFlag" >> large_instance_2.txt
	cp data_large/data_${T}_large/* data/
	./gep $bendersFlag $imprvdBendersFlag 10000
done



#	for n in 5
#	do
#		cd data
#		echo $n > numFWsample.dat
#		cd ..
#		./gep $bendersFlag $imprvdBendersFlag
