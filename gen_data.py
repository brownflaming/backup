#!/usr/bin/env python

import numpy as np
import math as m
import os
import pandas as pd
import pdb
from scipy.sparse import find, issparse, isspmatrix
from scipy.stats import beta


def write(filename, object):
	myFile = open(filename, "w")
	myFile.write(object)
	myFile.close()


if __name__ == "__main__":

	numStage = 3
	scenPerStage = 3
	numScen = np.ones(numStage-1, np.int32) * scenPerStage
	numScen = np.insert(numScen, 0, 1)
	numFWsample = 1

	ODI = 12 # AH HA BH HB CH HC ABH BHA AHCCHA BHC CHB
	FC = 6 # B1 B2 E1 E2 E3 E4
	LEG = 6 # AH HA BH HB CH HC
	CABIN = 2
	SEAT = np.array([24, 216])
	C0 = np.tile(SEAT, 6)
	c1 = np.array([500, 340, 200, 160, 130, 100])
	c2 = np.array([800, 540, 320, 260, 210, 160])
	clnRate = np.array([0.1, 0.1, 0.05, 0.05, 0.0, 0.0]) # for each FC
	rI = np.diag(np.tile(clnRate, ODI))  # diag(r,...r)
	rIinv = np.diag(np.tile((1-clnRate)**(-1), ODI)) # diag((1-r)^{-1}, ..., (1-r)^{-1})

	df = pd.read_csv("resource.csv", sep=',', header=None)
	df.fillna(0, inplace=True)
	R = df.values
	print R.shape

	nrow = ODI * FC * 4 + CABIN * LEG
	x_ub = SEAT * 1.1 ## upper bound on x variables, this is problem dependent
	numBin = np.floor(np.log2(np.ceil(x_ub))).astype(int) + 1  # number of binary needed for integer part
	dimB = ODI * FC
	# dimX_old = dimB * 2 # (B, C, P)
	dimX = (numBin[0] * 2 + numBin[1] * 4) * ODI  # (10+32)*12
	T = np.zeros((dimB, dimX))   # binary expansion matrix hat{x} = T*x
	print T.shape
	# pdb.set_trace()
	for i in range(ODI):
		for j in range(FC):
			row = i*FC + j
			if j < 2:
				for k in range(numBin[0]):
					col = i*dimX/ODI + j*numBin[0] + k
					T[row][col] = 2 ** k
			else:
				for k in range(numBin[1]):
					col = i*dimX/ODI + 2*numBin[0] + (j-2)*numBin[1] + k
					T[row][col] = 2 ** k
	# pdb.set_trace()

	initState = np.zeros(dimX)
	thetaLB = -(800.0 * 24 + 320 * 216) * 6 * np.ones(numStage)
	# temp = np.tile([-0.02*epsilon, 0.02*epsilon], (numStock, 1))
	# constrSlack = np.row_stack((temp, np.zeros((nrow - numStock, 2))))
	# constrSlack = np.transpose(constrSlack)
	# print constrSlack

	# dataDir = "bin/data/" + str(numStage) + "_" + str(scenPerStage) + "/"
	dataDir = "bin/data/"
	if not os.path.exists(dataDir):
		os.makedirs(dataDir)

	write(dataDir + "numStage.dat", str(numStage))
	write(dataDir + "numScen.dat", str(numScen.tolist()))
	write(dataDir + "numFWsample.dat", str(numFWsample))
	write(dataDir + "initState.dat", str(initState.tolist()))
	write(dataDir + "thetaLB.dat", str(thetaLB.tolist()))
	# write(dataDir + "constrSlack.dat", str(constrSlack.tolist()))
	

	#######################################################################
	# coefficient for x
	xObj = np.zeros((numStage, dimX)) # B in binary
	write(dataDir + "x.dat", str(xObj.tolist()))
	print "xObj: ", xObj.shape
	# print xObj

	# coefficient for y1
	y1Obj = np.zeros((numStage, dimB * 3)) # [Bt+1, Bt, bt+1, Ct+1, Ct, ct+1]
	cc = np.concatenate((np.tile(c1, ODI/2), np.tile(c2, ODI/2)))
	for t in range(1,numStage):
		y1Obj[t][dimB * 2: dimB * 3] = - cc
		# y1Obj[t][dimB * 5: dimB * 6] = cc
	write(dataDir + "y1.dat", str(y1Obj.tolist()))
	print "y1Obj: ", y1Obj.shape

	# coefficient for y2 [ ]
	y2Obj = np.zeros((numStage, 0))
	write(dataDir + "y2.dat", str(y2Obj.tolist()))
	print "y2Obj: ", y2Obj.shape


	## matrix of x variables
	A = []
	matX = np.row_stack((T, np.zeros((nrow - dimB, dimX))))
	print "matX", matX.shape
	# pdb.set_trace()
	matX = np.array(find(matX)).tolist()
	for t in range(0, numStage):
		A.append(matX)
	write(dataDir + "A.dat", str(A))
	print "A: ", len(A)
	for i in A: print len(i[0]),
	print " "

	## matrix of z variables
	B = []
	matZ = np.row_stack(( np.zeros((dimB, dimX)), T, np.zeros((nrow - dimB * 2, dimX))))
	print "matZ", matZ.shape
	matZ = np.array(find(matZ)).tolist()
	for t in range(numStage):
		B.append(matZ)
	write(dataDir + "B.dat", str(B))
	print "B: ", len(B)
	for i in B: print len(i[0]),
	print " "

	## matrix of y1 (local integer) variables [Bt+1, Bt, bt+1, Ct+1, Ct, ct+1]
	m1 = np.column_stack(( np.zeros((dimB, dimB * 0)), np.identity(dimB), np.zeros((dimB, dimB * 2)) ))
	m2 = np.column_stack(( np.zeros((dimB, dimB * 1)), np.identity(dimB), np.zeros((dimB, dimB * 1)) ))
	m3 = np.column_stack(( np.zeros((dimB, dimB * 2)), np.identity(dimB), np.zeros((dimB, dimB * 0)) ))
	# m4 = np.column_stack(( np.zeros((dimB, dimB * 3)), np.identity(dimB), np.zeros((dimB, dimB * 2)) ))
	# m5 = np.column_stack(( np.zeros((dimB, dimB * 4)), np.identity(dimB), np.zeros((dimB, dimB * 1)) ))
	# m6 = np.column_stack(( np.zeros((dimB, dimB * 5)), np.identity(dimB), np.zeros((dimB, dimB * 0)) ))
	m7 = np.column_stack(( R, np.zeros((R.shape[0], dimB * 2)) ))

	W1 = []
	# matY1 = np.row_stack((-m1, -m2, m1-m2-m3, m4-m5-m6, np.dot(rI, m1)-m4, -np.dot(rI, m1)+m4,\
						  # np.dot(rI, m2)-m5, -np.dot(rI, m2)+m5, m7, m3 ))
	matY1 = np.row_stack(( -m1, -m2, m1-m2-m3, m7, m3 ))
	print "matY1", matY1.shape, T.shape
	matY1 = np.array(find(matY1)).tolist()
	for t in xrange(numStage):
		W1.append(matY1)
	write(dataDir + "W1.dat", str(W1))
	print "W1: ", len(W1)
	for i in W1: print len(i[0]),
	print " "

	## matrix of y2 (local continuous) variables
	W2 = []
	matY2 = np.array(find(np.zeros((nrow, 0))))
	for t in xrange(numStage):
		W2.append(matY2.tolist())
	write(dataDir + "W2.dat", str(W2))
	print "W2: ", len(W2)
	for i in W2: print len(i[0]),
	print " "

	## slack matrix S
	S = np.identity(nrow)
	for i in range(dimB * 3):
		S[i][i] = 0
	S = np.array(find(S)).tolist()
	write(dataDir + "S.dat", str(S))

	## rhs of constraints
	b = np.zeros((numStage, nrow))
	rhs = np.concatenate(( np.zeros(dimB*3), C0, np.zeros(dimB) ))
	for t in xrange(numStage):
		b[t] = rhs
	write(dataDir + "rhs.dat", str(b.tolist()))
	print "b: ", b.shape
	# pdb.set_trace()
	# print rhs

	## uncertainty data index
	uncertainData = np.array([0, 0, 0, 0, 0, 0, 0, 1])
	write(dataDir + "uncertainData.dat", str(uncertainData.tolist()))

	# ####################################################################
	shape1 = np.array([3,3,10,40/3.0,22,30])
	shape2 = np.array([2,2,5,20/3.0,11,15])
	shape3 = np.array([2,2,7.5,10,16.5,22.5])
	shape4 = np.array([3,3,15,20,33,45])

	gammaShape = np.row_stack((np.tile(shape1, (4,1)), np.tile(shape2,(2,1)),\
							np.tile(shape3,(2,1)), np.tile(shape4, (4,1))))
	gammaScale = np.tile(np.array([1.5,1.5,1.2,1.2,1.0,1.0]), (ODI, 1))

	betaShape = np.array([[12,8,6,4,3,2],[1.5,2.0,2.0,3.0,4.0,4.0]])
	# dcp = (182 - np.array([182, 126, 84, 56, 35, 21, 14, 10, 7, 5, 3, 2, 1, 0]))/182.0
	# dcp = (182 - np.array([182, 126, 84, 56, 35, 21, 14, 10, 7, 5, 2, 0]))/182.0
	# dcp = (182 - np.array([182, 35, 7, 0]))/182.0
	dcp = np.arange(0,183,182/(numStage-1))/182.0
	prop = np.zeros((FC, numStage))
	for i in range(FC):
		prop[i] = beta.cdf(dcp, betaShape[0][i], betaShape[1][i])
		for j in range(len(prop[i])-1, 0, -1):
			prop[i,j] = prop[i,j] - prop[i,j-1]

	col = list('ij') + range(numStage)
	df = pd.DataFrame(columns = col)
	M = 200
	for i in range(ODI):
		for j in range(FC):
			arrival = np.random.gamma(gammaShape[i,j], gammaScale[i,j], M) / (1.0 - clnRate[j])
			path = np.round(np.outer(arrival, prop[j]))
			path = np.column_stack((np.tile([i,j], (M,1)), path))
			dfnew = pd.DataFrame(data=path, columns=col)
			df = df.append(dfnew, ignore_index=True)
	# print df.shape

	scenTree = np.zeros((numStage, scenPerStage, ODI * FC))
	for i in range(ODI):
		for j in range(FC):
			dfnew = df.loc[(df['i'] == i) & (df['j'] == j), :]
			for t in range(1,numStage):
				scenTree[t,:,i * FC + j] = dfnew.loc[:,t].sample(numScen[t]).values
	
	# print scenTree
	# pdb.set_trace()

	## bRhsScenarios (T x scenPerStage x numRows)
	bScen = np.zeros((numStage, scenPerStage, nrow))
	for t in range(1, numStage):
		for k in xrange(scenPerStage):
			new_rhs = rhs
			new_rhs[-dimB : ] = scenTree[t,k]
			bScen[t,k] = new_rhs
	write(dataDir + "bScen.dat", str(bScen.tolist()))
	print "bScen: ", bScen.shape

	newScenTree = scenTree.reshape((numStage, scenPerStage, ODI, FC))
	# print newScenTree
	# dataDir = "tree/data/" + str(numStage) + "_" + str(scenPerStage) + "/"
	dataDir = "tree/data/"
	if not os.path.exists(dataDir):
		os.makedirs(dataDir)
	write(dataDir + "demand.dat", str(newScenTree.tolist()))
	write(dataDir + "stage.dat", str(numStage))
	write(dataDir + "branch.dat", str(scenPerStage))

	####################################################################
	# extensive formulation data
	# dataDir = "tree/data/" + str(numStage) + "_" + str(scenPerStage) + "/"
	# dataDir = "tree/data/"
	# if not os.path.exists(dataDir):
		# os.makedirs(dataDir)

	# write(dataDir + "numStage.dat", str(numStage))
	# write(dataDir + "numStock.dat", str(numStock))
	# write(dataDir + "numChi.dat", str(scenPerStage))
	# write(dataDir + "assetLimit.dat", str(assetLimit.tolist()))
	# write(dataDir + "buyTran.dat", str(buy.tolist()))
	# write(dataDir + "sellTran.dat", str(sell.tolist()))
	# write(dataDir + "initState.dat", str(np.dot(T, initState).tolist()))
	# write(dataDir + "scenario.dat", str(scenTree.tolist()))