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

	numStage = 14
	scenPerStage = 3
	numScen = np.ones(numStage, np.int32) * scenPerStage
	numScen = np.insert(numScen, 0, 1)
	numFWsample = 2

	ODI = 12 # AH HA BH HB CH HC ABH BHA AHCCHA BHC CHB
	FC = 6 # B1 B2 E1 E2 E3 E4
	LEG = 6 # AH HA BH HB CH HC
	CABIN = 2
	SEAT = np.array([24, 216])
	C0 = np.tile(SEAT, 6)
	c1 = np.array([500,340,200,160,130,100])
	c2 = np.array([800,540,320,260,210,160])
	clnRate = np.array([0.1,0.1,0.05,0.05,0.0,0.0]) # for each FC
	rI = np.diag(np.tile(clnRate, ODI))  # diag(r,...r)
	rIinv = np.diag(np.tile((1-clnRate)**(-1), ODI)) # diag((1-r)^{-1}, ..., (1-r)^{-1})

	df = pd.read_csv("matrix.csv", sep=',', header=None)
	df.fillna(0, inplace=True)
	R = df.values
	print R.shape

	nrow = ODI * FC * 10 + CABIN * LEG
	x_ub = SEAT * 2 ## upper bound on x variables, this is problem dependent
	numBin = np.floor(np.log2(np.ceil(x_ub))).astype(int) + 1  # number of binary needed for integer part
	dimB = ODI * FC
	dimX_old = ODI * FC * 3 # (B, P)
	dimX = (numBin[0] * 2 + numBin[1] * 4) * ODI * 3  # (12+36)*12*3
	T = np.zeros((dimX_old, dimX))   # binary expansion matrix hat{x} = T*x
	print T.shape
	# pdb.set_trace()
	for var in range(3):
		for i in range(ODI):
			for j in range(FC):
				row = var*ODI*FC + i*FC + j
				if j < 2:
					for k in range(numBin[0]):
						col = var*dimX/3 + i*dimX/ODI/3 + j*numBin[0] + k
						T[row][col] = 2 ** k
				else:
					for k in range(numBin[1]):
						col = var*dimX/3 + i*dimX/ODI/3 + 2*numBin[0] + (j-2)*numBin[1] + k
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
	xObj = np.zeros((numStage, dimX)) # [B C P] in binary
	write(dataDir + "x.dat", str(xObj.tolist()))
	print "xObj: ", xObj.shape
	# print xObj

	# coefficient for y1
	y1Obj = np.zeros((numStage, dimB * 3)) # [b c u]
	cc = np.concatenate((np.tile(c1, 6), np.tile(c2, 6)))
	for t in range(1, numStage):
		y1Obj[t][0:len(cc)] = - cc
		y1Obj[t][len(cc):2*len(cc)] = cc
	write(dataDir + "y1.dat", str(y1Obj.tolist()))
	print "y1Obj: ", y1Obj.shape

	# coefficient for y2 [Up, Ud]
	y2Obj = np.zeros((numStage, dimB*4))
	penalty = 1e10 * np.ones(dimB * 2)
	for t in range(1,numStage):
		y2Obj[t][dimB*2:dimB*4] = penalty 
	write(dataDir + "y2.dat", str(y2Obj.tolist()))
	print "y2Obj: ", y2Obj.shape

	## matrix of x variables
	m0 = np.zeros((dimB, dimB * 3))
	m1 = np.column_stack((np.identity(dimB), np.zeros((dimB, 2 * dimB))))
	m2 = np.column_stack((np.zeros((dimB, dimB)), np.identity(dimB), np.zeros((dimB, dimB))))
	m3 = np.column_stack((np.zeros((dimB, 2 * dimB)), np.identity(dimB)))
	m4 = np.column_stack((rI, np.zeros((dimB, 2 * dimB))))
	m5 = np.column_stack((np.zeros((CABIN*LEG,dimB*2)), R))
	# print m0.shape, m1.shape, m2.shape, m3.shape, m4.shape, m5.shape
	matX = np.row_stack((np.tile(m0, (10,1)), m5))
	matX = np.array(find(np.dot(matX, T)))
	A = [matX.tolist()]
	matX = np.row_stack((m1, m2, m0, m1, -m1, -m4 + m2, m4 - m2, m0, m0, m0, m5))
	print matX.shape, T.shape
	# pdb.set_trace()
	matX = np.array(find(np.dot(matX, T))).tolist()
	for t in range(1, numStage):
		A.append(matX)
	write(dataDir + "A.dat", str(A))
	print "A: ", len(A)
	for i in A: print len(i[0]),
	print " "

	## matrix of z variables
	m6 = np.column_stack((np.zeros((dimB, 2 * dimB)), rIinv))
	matZ = np.array(find(np.zeros((nrow, dimX))))
	B = [matZ.tolist()]
	matZ = np.row_stack((-m1, -m2, m0, -m6, m6, np.tile(m0,(5,1)), np.zeros((CABIN*LEG, dimB*3))))
	matZ = np.array(find(np.dot(matZ, T))).tolist()
	for t in range(1, numStage):
		B.append(matZ)
	write(dataDir + "B.dat", str(B))
	print "B: ", len(B)
	for i in B: print len(i[0]),
	print " "


	## matrix of y1 (local integer) variables
	m7 = np.column_stack((np.zeros((dimB, 2 * dimB)), sum(x_ub) * np.identity(dimB)))
	matY1 = np.array(find(np.zeros((nrow, dimB * 3))))
	W1 = [matY1.tolist()]
	matY1 = np.row_stack((-m1, -m2, m1, np.tile(m0, (4,1)), m7, -m7, m3, np.zeros((CABIN*LEG, dimB*3))))
	matY1 = np.array(find(matY1)).tolist()
	for t in xrange(1, numStage):
		W1.append(matY1)
	write(dataDir + "W1.dat", str(W1))
	print "W1: ", len(W1)
	for i in W1: print len(i[0]),
	print " "

	## matrix of y2 (local continuous) variables
	s0 = np.zeros((dimB, dimB * 4))
	s1 = np.column_stack((np.identity(dimB), np.zeros((dimB, dimB*3))))
	s2 = np.column_stack((np.zeros((dimB, dimB)), np.identity(dimB), np.zeros((dimB,dimB*2)) ))
	s3 = np.column_stack((np.zeros((dimB, dimB*2)), np.identity(dimB), np.zeros((dimB,dimB)) ))
	s4 = np.column_stack((np.zeros((dimB, dimB*3)), np.identity(dimB)))
	matY2 = np.array(find(np.zeros((nrow, dimB * 4))))
	W2 = [matY2.tolist()]
	matY2 = np.row_stack((s0, s4, s2, s1-s3, -s1+s3, np.tile(s0,(2,1)), s1, s2, s0, np.zeros((CABIN*LEG, dimB*4))))
	matY2 = np.array(find(matY2)).tolist()
	for t in xrange(1, numStage):
		W2.append(matY2)
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
	rhs0 = np.concatenate((np.zeros(nrow-CABIN*LEG), C0 * 2))
	rhs1 = np.concatenate((np.zeros(dimB*4),np.ones(dimB), np.zeros(dimB),\
						np.ones(dimB), sum(x_ub) * np.ones(dimB), np.zeros(dimB),\
						np.ones(dimB), C0 ))
	rhs2 = np.concatenate((np.zeros(dimB*4),np.ones(dimB), np.zeros(dimB),\
						np.ones(dimB), sum(x_ub) * np.ones(dimB), np.zeros(dimB),\
						np.ones(dimB), C0 * 2 ))
	for t in xrange(numStage):
		if t == 0:
			b[t] = rhs0
		elif (t == numStage - 2):
			b[t] = rhs1
		else:
			b[t] = rhs2
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
	dcp = (182 - np.array([182, 126, 84, 56, 35, 21, 14, 10, 7, 5, 3, 2, 1, 0]))/182.0
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
			if (t == numStage - 2):
				new_rhs = rhs1
				new_rhs[dimB*2 : dimB*3] = scenTree[t,k]
				bScen[t,k] = new_rhs
			else:
				new_rhs = rhs2
				new_rhs[dimB*2 : dimB*3] = scenTree[t,k]
				bScen[t,k] = new_rhs
	write(dataDir + "bScen.dat", str(bScen.tolist()))
	print "bScen: ", bScen.shape

	

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