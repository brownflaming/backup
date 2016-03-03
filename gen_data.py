#!/usr/bin/env python

import numpy as np
import math as m
import os
import pandas as pd
import pdb


def write(filename, object):
	myFile = open(filename, "w")
	myFile.write(object)
	myFile.close()


if __name__ == "__main__":

	numStage = 5
	scenPerStage = 20
	numScen = np.ones(numStage, np.int32) * scenPerStage
	numScen = np.insert(numScen, 0, 1)
	numFWsample = 2
	numStock = 80
	nrow = 2* numStock + 2
	buy = 0.00 * np.ones(numStock) # transaction fee
	sell = 0.00 * np.ones(numStock) # transaction fee

	x_ub = 100 ## upper bound on x variables, this is problem dependent
	epsilon = 1e-2 ## accuracy of binary approximation
	k1 = int(m.log(np.ceil(x_ub)-1, 2)) + 1  # number of binary needed for integer part
	k2 = int(m.log(1.0/epsilon, 2)) + 1 # number of binary needed for decimal part
	k = k1 + k2
	print k1, k2
	
	assetLimit = int(0.5 * numStock) # transaction limit

	dimX = k * numStock
	T = np.zeros((numStock, dimX))   # binary expansion matrix hat{x} = T*x
	for i in range(numStock):
		for j in range(k1, -k2, -1):
			T[i][k1 - j + k * i] = 2**(j-1)
	# print T

	initState = np.zeros(k * numStock)
	initState[2] = 1
	#initState[k1-1] = 1
	thetaLB = -500.0 * np.ones(numStage)
	temp = np.tile([-0.2*epsilon, 0.2*epsilon], (numStock, 1))
	constrSlack = np.row_stack((temp, np.zeros((nrow-numStock, 2))))
	constrSlack = np.transpose(constrSlack)
	# print constrSlack

	dataDir = "bin/data/"
	if not os.path.exists(dataDir):
		os.makedirs(dataDir)

	write(dataDir + "numStage.dat", str(numStage))
	write(dataDir + "numScen.dat", str(numScen.tolist()))
	write(dataDir + "numFWsample.dat", str(numFWsample))
	write(dataDir + "initState.dat", str(initState.tolist()))
	write(dataDir + "thetaLB.dat", str(thetaLB.tolist()))
	write(dataDir + "constrSlack.dat", str(constrSlack.tolist()))

	#######################################################################
	# coefficient for x
	xObj = np.zeros((numStage - 1, dimX))
	xObj = np.row_stack((xObj, np.dot(-np.ones(numStock), T)))
	write(dataDir + "x.dat", str(xObj.tolist()))
	print "xObj: ", xObj.shape
	# print xObj
	# pdb.set_trace()

	# coefficient for y1
	y1Obj = np.zeros((numStage, numStock))
	write(dataDir + "y1.dat", str(y1Obj.tolist()))
	print "y1Obj", y1Obj.shape

	# coefficient for y2
	y2Obj = np.zeros((numStage, numStock * 2))
	write(dataDir + "y2.dat", str(y2Obj.tolist()))
	print "y2Obj", y2Obj.shape

	## matrix of x variables
	A = np.zeros((numStage, nrow, dimX))
	matX = np.row_stack((np.identity(numStock), np.identity(numStock), np.zeros((2, numStock))))
	matX = np.dot(matX, T)
	for t in xrange(numStage):
		A[t] = matX
	write(dataDir + "A.dat", str(A.tolist()))
	print "A: ", A.shape

	## matrix of z variables
	B = np.zeros((numStage, nrow, dimX))
	matZ = np.row_stack((-np.identity(numStock), np.zeros((numStock+2, numStock))))
	matZ = np.dot(matZ, T)
	for t in xrange(numStage):
		B[t] = matZ
	write(dataDir + "B.dat", str(B.tolist()))
	print "B: ", B.shape

	## matrix of y1 (local integer) variables
	W1 = np.zeros((numStage, nrow, numStock))
	matY1 = np.row_stack(( np.zeros((numStock, numStock)), - x_ub*np.identity(numStock),\
	 np.ones(numStock), np.zeros(numStock) ))
	for t in xrange(numStage):
		W1[t] = matY1
	write(dataDir + "W1.dat", str(W1.tolist()))
	print "W1: ", W1.shape

	## matrix of y2 (local continuous) variables
	W2 = np.zeros((numStage, nrow, numStock * 2))
	matY2 = np.column_stack((-np.identity(numStock), np.identity(numStock)))
	matY2 = np.row_stack((matY2, np.zeros((numStock + 1, numStock * 2))))
	temp = np.concatenate((1 + buy, sell - 1))
	matY2 = np.row_stack((matY2, temp))
	for t in xrange(numStage):
		W2[t] = matY2
	write(dataDir + "W2.dat", str(W2.tolist()))
	print "W2: ", W2.shape

	## slack matrix S
	S = np.identity(nrow)
	for i in range(numStock): S[i][i] = 0
	S[-1][-1] = 0
	write(dataDir + "S.dat", str(S.tolist()))
	print "S: ", S.shape

	## rhs of constraints
	b = np.zeros((numStage, nrow))
	rhs = np.concatenate((np.zeros(2*numStock), [assetLimit, 0]))
	for t in xrange(numStage): b[t] = rhs
	write(dataDir + "rhs.dat", str(b.tolist()))
	print "b: ", b.shape
	# print rhs

	## uncertainty data index
	uncertainData = np.array([1, 0, 0, 0, 1, 0, 0, 0])
	write(dataDir + "uncertainData.dat", str(uncertainData.tolist()))

	####################################################################
	universe = pd.read_csv("return.csv")
	idxLength = len(universe.index)
	scenTree = np.ones((numStage+1, scenPerStage, numStock))
	# stock_ind = list(np.floor(np.random.sample(numStock) * 99).astype(int))
	stock_ind = range(3, numStock+3)
	# print list(stock_ind.astype(int))
	# print universe.ix[1,list(stock_ind.astype(int))]

	for t in range(1, numStage+1):
		# sample_ind = np.array(range(scenPerStage)) + (t-1) * scenPerStage
		# sample_ind = np.array(range(scenPerStage)) + (t+5) * scenPerStage
		sample_ind = np.ceil(np.random.sample(numScen[t]) * idxLength)
		# print np.array(universe.ix[sample_ind,stock_ind])
		# print sample_ind, stock_ind
		scenTree[t][:,0:numStock] = np.array(universe.ix[sample_ind,stock_ind])

	## xObjScenarios (T x scenPerStage x numStock)
	xScen = np.zeros((numStage, scenPerStage, dimX))
	temp = np.dot(np.sum(scenTree[-1], axis=0)/scenTree[-1].shape[0], T)
	for k in xrange(scenPerStage): xScen[-1][k] = - temp
	write(dataDir + "xScen.dat", str(xScen.tolist()))
	print "xScen: ", xScen.shape
	# print xScen[1]

	## BmatrixScenarios (T x scenPerStage, 2 * numStock, numStock)
	BScen = np.zeros((numStage, scenPerStage, nrow, dimX))
	for t in xrange(numStage):
		for k in xrange(scenPerStage):
			matZ = np.row_stack((-np.diag(scenTree[t][k]), np.zeros((numStock+2, numStock))))
			BScen[t][k] = np.dot(matZ, T)
	write(dataDir + "BScen.dat", str(BScen.tolist()))
	print "BScen: ", BScen.shape
	# print BScen
	print scenTree

	####################################################################
	# extensive formulation data
	dataDir = "tree/data/"
	if not os.path.exists(dataDir):
		os.makedirs(dataDir)

	write(dataDir + "numStage.dat", str(numStage))
	write(dataDir + "numStock.dat", str(numStock))
	write(dataDir + "numChi.dat", str(scenPerStage))
	write(dataDir + "assetLimit.dat", str(assetLimit))
	write(dataDir + "buyTran.dat", str(buy))
	write(dataDir + "sellTran.dat", str(sell))
	write(dataDir + "initState.dat", str(np.dot(T, initState).tolist()))
	write(dataDir + "scenario.dat", str(scenTree.tolist()))



