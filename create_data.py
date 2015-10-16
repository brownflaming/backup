#!/usr/bin/env python
import numpy as np

if __name__ == "__main__":

    HORIZON = 3
    GENERATOR_LIST = {0: 'BaseLoad', 1: 'CC', 2: 'CT', 3: 'Nuclear', 4: 'Wind', 5: 'IGCC'}
    MAX_OUTPUT = np.array([1130.0, 390.0, 380.0, 1180.0, 175.0, 560.0])
    MAX_UNIT = np.array([4, 10, 10, 1, 45, 4], np.int32)
    SUBPERIOD = 3
    SUBPERIOD_HOUR = [271.0, 6556.0, 1933.0]

    CONSTRUCTION_COST = np.array([1.446, 0.795, 0.575, 1.613, 1.650, 1.671])
    MAX_CAPACITY = np.array([1200.0, 400.0, 400.0, 1200.0, 500.0, 600.0])

    FUEL_PRICE = np.array([3.37, 9.11 * 1e-6 / 1.028, 9.11 * 1e-6 / 1.028, 0.93e-3, 0, 3.37])
    RATIO = [8.844 / 0.4, 7.196 / 0.56, 10.842 / 0.4, 10.400 / 0.45, 0.0, 8.613 / 0.48]
    FUEL_PRICE_GROWTH = 0.02
    OPERATING_COST = np.array([4.7, 2.11, 3.66, 0.51, 5.00, 2.98]) * 1e-6
    OPER_COST_GROWTH = 0.03
    PENALTY_COST = 1e-1
    LAMBDA = np.array([1.38, 1.04, 0.80])
    HOURS_PER_YEAR = 8760.0
    rate = 0.08     # interest rate

    nType = len(GENERATOR_LIST)     # number of technology
    nUnit = sum(MAX_UNIT)           # sum of maximum construction

    numFWsample = 20

    numScen = np.ones(HORIZON - 1, np.int32) * 5
    numScen = np.insert(numScen, 0, 1)

    initState = np.zeros(2 * (nUnit + nType))
    for i in xrange(nType):
        index = nUnit + nType + sum(MAX_UNIT[0: i]) + i
        initState[index] = 1

    valueLB = np.zeros(HORIZON)

    D0 = 0.57
    scenarios = [[0] * numScen[i] for i in xrange(HORIZON)]
    scenarios[0][0] = [0] * 2 * SUBPERIOD
    for t in xrange(1, HORIZON):
        sample = D0 * (1.01 ** t - 1.005 ** t) * np.random.sample((5,)) + D0 * 1.005 ** t
        scenarios[t] = list(sample)
        for j in xrange(numScen[t]):
            temp = LAMBDA * max(scenarios[t][j] - D0, 0.0) * 1e9 / HOURS_PER_YEAR
            temp = np.concatenate((temp, -temp))
            scenarios[t][j] = list(temp)

    ''' Parameters used by SDDP '''
    # number of stages
    myFile = open("bin/data/numStage.dat", "w")
    myFile.write(str(HORIZON))
    myFile.close()

    # initial state variable
    myFile = open("bin/data/initState.dat", "w")
    myFile.write("[")
    for i in xrange(initState.shape[0]):
        myFile.write(str(initState[i]))
        if i != initState.shape[0] - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    # number of forward samples
    myFile = open("bin/data/numFWsample.dat", "w")
    myFile.write(str(numFWsample))
    myFile.close()

    # number of scenarios at each stage
    myFile = open("bin/data/numScen.dat", "w")
    myFile.write("[")
    for i in xrange(numScen.shape[0]):
        myFile.write(str(numScen[i]))
        if i != numScen.shape[0] - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    # lower bounds for value functions
    myFile = open("bin/data/valueLB.dat", "w")
    myFile.write("[")
    for i in xrange(valueLB.shape[0]):
        myFile.write(str(valueLB[i]))
        if i != valueLB.shape[0] - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    # scenarios
    myFile = open("bin/data/scenarios.dat", "w")
    myFile.write(str(scenarios))
    myFile.close()

    ''' Construct data for objective function c x + b y '''
    # construct c
    totalCost = np.multiply(CONSTRUCTION_COST, MAX_CAPACITY)
    xCoef = np.zeros((HORIZON, 2 * (nUnit + nType)))
    for t in np.arange(HORIZON):
        temp = np.arange(MAX_UNIT[0] + 1) * totalCost[0] / (1 + rate) ** t
        for i in np.arange(1, nType):
            temp = np.concatenate((temp, np.arange(MAX_UNIT[i] + 1) * totalCost[i] / (1 + rate) ** t))
        temp = np.concatenate((temp, np.zeros(nUnit + nType)))
        xCoef[t] = temp

    print "c.shape", xCoef.shape

    myFile = open("bin/data/xCoef.dat", "w")
    myFile.write("[")
    for row in xrange(xCoef.shape[0]):
        myFile.write("[")
        for column in xrange(xCoef.shape[1]):
            myFile.write(str(xCoef[row, column]))
            if column != xCoef.shape[1] - 1:
                myFile.write(",")
        myFile.write("]")
        if row != xCoef.shape[0] - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    # construct b
    yCoef = np.zeros((HORIZON, SUBPERIOD * (nType + 1)))
    for t in np.arange(HORIZON):
        for k in np.arange(SUBPERIOD):
            temp = np.multiply(FUEL_PRICE, RATIO) * (1 + FUEL_PRICE_GROWTH) ** t + \
                   OPERATING_COST * (1 + OPER_COST_GROWTH) ** t
            temp = np.append(temp, PENALTY_COST)
            temp = SUBPERIOD_HOUR[k] * temp / (1 + rate) ** t
            yCoef[t, k * (nType + 1): (k + 1) * (nType + 1)] = temp

    print "b.shape", yCoef.shape

    myFile = open("bin/data/yCoef.dat", "w")
    myFile.write("[")
    for row in xrange(yCoef.shape[0]):
        myFile.write("[")
        for column in xrange(yCoef.shape[1]):
            myFile.write(str(yCoef[row, column]))
            if column != yCoef.shape[1] - 1:
                myFile.write(",")
        myFile.write("]")
        if row != yCoef.shape[0] - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    ''' Construct data for constraints A x + W y + B z >= rhs '''

    # construct matrix A (x variable)

    P = np.zeros((nType, nUnit + nType))
    for i in xrange(nType):
        P[i, sum(MAX_UNIT[0:i]) + i: sum(MAX_UNIT[0:i + 1]) + i + 1] = range(0, MAX_UNIT[i] + 1)

    Q = np.zeros((nType, nUnit + nType))
    for i in xrange(nType):
        Q[i, sum(MAX_UNIT[0:i]) + i: sum(MAX_UNIT[0:i + 1]) + i + 1] = np.ones(MAX_UNIT[i] + 1)

    matrix0 = np.zeros((nType, nUnit + nType))

    pattern1 = np.hstack((matrix0, P))
    pattern2 = np.hstack((-P, P))
    pattern3 = np.hstack((Q, matrix0))
    pattern4 = np.hstack((matrix0, Q))
    pattern5 = np.zeros((2 * SUBPERIOD, 2 * (nUnit + nType)))

    A = np.vstack((pattern1, pattern1, pattern1,
                   pattern2, -pattern2, -pattern1,
                   pattern3, -pattern3, pattern4, -pattern4, pattern5))

    print "A.shape:", A.shape

    myFile = open("bin/data/Amatrix.dat", "w")
    myFile.write("[")
    for t in np.arange(HORIZON):
        myFile.write("[")
        for row in xrange(A.shape[0]):
            myFile.write("[")
            for column in xrange(A.shape[1]):
                myFile.write(str(A[row, column]))
                if column != A.shape[1] - 1:
                    myFile.write(",")
            myFile.write("]")
            if row != A.shape[0] - 1:
                myFile.write(",")
        myFile.write("]")
        if t != HORIZON - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    # construct matrix W (y variable)
    temp = -1.0 / MAX_OUTPUT
    R = np.diagflat([temp, temp, temp])
    for k in xrange(SUBPERIOD):
        R = np.insert(R, nType * (k + 1), 0.0, axis=1)

    row1 = np.hstack((np.ones(nType + 1),
                      np.zeros((SUBPERIOD - 1) * (nType + 1))))
    row2 = np.hstack((np.zeros(nType + 1),
                      np.ones(nType + 1),
                      np.zeros(nType + 1)))
    row3 = np.hstack((np.zeros((SUBPERIOD - 1) * (nType + 1)),
                      np.ones(nType + 1)))

    W = np.vstack((R, np.zeros((7 * nType, SUBPERIOD * (nType + 1))),
                   row1, row2, row3, -row1, -row2, -row3))

    print "W.shape:", W.shape

    myFile = open("bin/data/Wmatrix.dat", "w")
    myFile.write("[")
    for t in np.arange(HORIZON):
        myFile.write("[")
        for row in xrange(W.shape[0]):
            myFile.write("[")
            for column in xrange(W.shape[1]):
                myFile.write(str(W[row, column]))
                if column != W.shape[1] - 1:
                    myFile.write(",")
            myFile.write("]")
            if row != W.shape[0] - 1:
                myFile.write(",")
        myFile.write("]")
        if t != HORIZON - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    # construct B matrix (z variables)

    B = np.vstack((np.zeros((3 * nType, nUnit + nType)), -P, P,
                   np.zeros((5 * nType + 2 * SUBPERIOD, nUnit + nType))))

    B = np.hstack((np.zeros((B.shape[0], B.shape[1])), B))

    print "B.shape:", B.shape

    myFile = open("bin/data/Bmatrix.dat", "w")
    myFile.write("[")
    for t in np.arange(HORIZON):
        myFile.write("[")
        for row in xrange(B.shape[0]):
            myFile.write("[")
            for column in xrange(B.shape[1]):
                myFile.write(str(B[row, column]))
                if column != B.shape[1] - 1:
                    myFile.write(",")
            myFile.write("]")
            if row != B.shape[0] - 1:
                myFile.write(",")
        myFile.write("]")
        if t != HORIZON - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    # construct rhs
    rhs = np.hstack((np.zeros((5 * nType)), -MAX_UNIT,
                     np.ones(nType), -np.ones(nType),
                     np.ones(nType), -np.ones(nType),
                     np.zeros(2 * SUBPERIOD)))

    print "rhs.shape:", rhs.shape

    myFile = open("bin/data/bRhs.dat", "w")
    myFile.write("[")
    for t in np.arange(HORIZON):
        myFile.write("[")
        for row in xrange(rhs.shape[0]):
            myFile.write(str(rhs[row]))
            if row != rhs.shape[0] - 1:
                myFile.write(",")
        myFile.write("]")
        if t != HORIZON - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()

    # uncertain indices
    uncertainIndex = np.arange(rhs.shape[0] - 6, rhs.shape[0])

    myFile = open("bin/data/uncertainIndex.dat", "w")
    myFile.write("[")
    for i in xrange(uncertainIndex.shape[0]):
        myFile.write(str(uncertainIndex[i]))
        if i != uncertainIndex.shape[0] - 1:
            myFile.write(",")
    myFile.write("]")
    myFile.close()
