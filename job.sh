#!/bin/bash

cd bin

bendersFlag=1
imprvdBendersFlag=1

for K in 3 5 10 15 20 30 40 50
do
	cp test_data/5_${K}/* data/
	./gep $bendersFlag $imprvdBendersFlag
done


#for T in $(seq 5 9)
#do
#	echo "Time horizon: $T; Benders' cut: $bendersFlag; Improved Benders' cut: $imprvdBendersFlag" >> large_instance_2.txt
#	cp data_large/data_${T}_large/* data/
#	./gep $bendersFlag $imprvdBendersFlag
#done



#	for n in 5
#	do
#		cd data
#		echo $n > numFWsample.dat
#		cd ..
#		./gep $bendersFlag $imprvdBendersFlag
