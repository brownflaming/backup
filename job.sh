#!/bin/bash

cd bin

for T in 3 4 5
do
	for K in 10 15 20
	do
		cp data/${T}_${K}/* data/
		./pfopt 0 1 0 1
		clear
	done
done

# cd ../tree

# for T in 3 4
# do
# 	for K in 10 15 25
# 	do
# 		cp data/${T}_${K}/* data/
# 		./tree
# 		clear
# 	done
# done


# cd bin
# for K in 10 15 25
# do
# 	cp data/3_${K}/* data/
# 	./pfopt 0 1 0 1
# 	clear
# done

# cd ../tree
# for K in 10 15 25
# do
# 	cp data/3_${K}/* data/
# 	./tree
# 	clear
# done



# cd bin

# cp inst_data/data_5/* data/
# ./gep 1 0 0 0
# ./gep 0 1 0 0
# ./gep 0 0 1 0
# ./gep 0 0 0 1
# ./gep 1 0 1 0
# ./gep 1 0 0 1
# ./gep 0 1 1 0
# ./gep 0 1 0 1
# ./gep 0 1 1 1

#########################################
# for T in $(seq 5 9)
# do
# 	cp data_large/data_${T}_large/* data/
# 	./gep 1 0 1 0
# 	./gep 0 1 1 0
# 	./gep 0 1 1 1
# done

########################################
# cp inst_data/data_10/* data/

# strengthened benders + lagrangian
# for K in 10 20 50
# do
# 	cd data
# 	echo $K > numFWsample.dat
# 	cd ..
# 	for i in $(seq 1 3)
# 	do
# 		./gep 0 1 1 0
# 	done
# done

# strengthened benders + integer
# for K in 1 2 3 5 10 20 50
# do
# 	cd data
# 	echo $K > numFWsample.dat
# 	cd ..
# 	for i in $(seq 1 3)
# 	do
# 		./gep 0 1 0 1
# 	done
# done

# # strengthened benders + lagrangian + integer
# for K in 1 2 3 5 10 20 50
# do
# 	cd data
# 	echo $K > numFWsample.dat
# 	cd ..
# 	for i in $(seq 1 2)
# 	do
# 		./gep 0 1 1 1
# 	done
# done

# # lagrangian
# for K in 1 2 3 5 10 20 50
# do
# 	cd data
# 	echo $K > numFWsample.dat
# 	cd ..
# 	for i in $(seq 1 1)
# 	do
# 		./gep 0 0 1 0
# 	done
# done

# # benders + lagrangian
# for K in 1 2 3 5 10 20 50
# do
# 	cd data
# 	echo $K > numFWsample.dat
# 	cd ..
# 	for i in $(seq 1 1)
# 	do
# 		./gep 1 0 1 0
# 	done
# done

# # benders + integer
# for K in 1 2 3 5 10 20 50
# do
# 	cd data
# 	echo $K > numFWsample.dat
# 	cd ..
# 	for i in $(seq 1 1)
# 	do
# 		./gep 1 0 0 1
# 	done
# done

# integer
# for K in 1 2 3 5 10 20 50
# do
#	cd data
#	echo $K > numFWsample.dat
#	cd ..
#	./gep 0 0 0 1
#done
#########################################

#bendersFlag=1
#imprvdBendersFlag=0
#lagrangianFlag=1
#integerFlag=0

#for T in $(seq 3 10)
#do
#	cp inst_data/data_${T}/* data/
#	for K in 1 2 3 5 10 20 50
#	do
#		cd data
#		echo $K > numFWsample.dat
#		cd ..
#		./gep $bendersFlag $imprvdBendersFlag $lagrangianFlag $integerFlag
#	done
#done

#for T in $(seq 10 10)
#do
#	cp data_large/data_${T}_large/* data/
#	./gep $bendersFlag $imprvdBendersFlag $lagrangianFlag $integerFlag
#done



#	for n in 5
#	do
#		cd data
#		echo $n > numFWsample.dat
#		cd ..
#		./gep $bendersFlag $imprvdBendersFlag
