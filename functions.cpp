#include <ilcplex/ilocplex.h>
#include <cmath>
#include <unordered_set>
#include <string>
#include <sstream>
#include <vector>
#include <cstddef>  // std::size_t
#include <numeric>
#include <random>

#include "global.h"
#include "functions.h"
#include "mt64.h"

using namespace std;


void sparse2full3d(IloNumArray3 & fullMat, IloInt dim [2], const char * fileName)
{
	IloEnv env = fullMat.getEnv();
	int stage = fullMat.getSize();
	IloNumArray3 sparseMat = IloNumArray3(env, stage);
	readArray<IloNumArray3> (sparseMat, fileName);
	for ( int t = 0; t < stage; ++t )
	{
		IloNumArray2 M = IloNumArray2(env, dim[0]);
		for ( int i = 0; i < dim[0]; ++i )
			M[i] = IloNumArray(env, dim[1]);
		int k = sparseMat[t][0].getSize(); // number of nonzero entries
		for ( int i = 0; i < k; ++i )
			M[int(sparseMat[t][0][i])][int(sparseMat[t][1][i])] = sparseMat[t][2][i];
		fullMat[t] = M;
	}
	sparseMat.end();
}

void sparse2full2d(IloNumArray2 & fullMat, IloInt dim [2], const char * fileName)
{
	IloEnv env = fullMat.getEnv();
	IloNumArray2 sparseMat = IloNumArray2(env, 3);
	readArray<IloNumArray2> (sparseMat, fileName);
	IloNumArray2 M = IloNumArray2(env, dim[0]);
	for ( int i = 0; i < dim[0]; ++i )
		M[i] = IloNumArray(env, dim[1]);
	int k = sparseMat[0].getSize(); // number of nonzero entries
	for ( int i = 0; i < k; ++i )
		M[int(sparseMat[0][i])][int(sparseMat[1][i])] = sparseMat[2][i];
	fullMat = M;
	sparseMat.end();
}

void readData (formatData * fData_p)
{
	// create data environment
	fData_p->dataEnv = IloEnv();
	IloEnv * dataEnv = & fData_p->dataEnv;

	cout << "Start reading data from files... " << endl;

	// read number of stages
	readArray<IloInt> (fData_p->numStage, "data/numStage.dat");
	IloInt numStage = fData_p->numStage;

	// read number of samples drawn in the forward pass (number of candidate solns)
	readArray<IloInt> (fData_p->numFWsample, "data/numFWsample.dat");

	// read number of scenarios available at each stage
	fData_p->numScen = IloIntArray(*dataEnv, numStage);
	readArray<IloIntArray> (fData_p->numScen, "data/numScen.dat");

	// compute the total number of scenarios
	fData_p->totalScen = fData_p->numScen[0];
	for (int t = 1; t < numStage; ++t)
		fData_p->totalScen *= fData_p->numScen[t];

	// read initial state of binary variables
	fData_p->initState = IloNumArray(*dataEnv, numStage);
	readArray<IloNumArray> (fData_p->initState, "data/initState.dat");

	// read prescribed lower bound for value function at each stage
	fData_p->thetaLB = IloNumArray(*dataEnv, numStage);
	readArray<IloNumArray> (fData_p->thetaLB, "data/thetaLB.dat");

	// read the tolerance of constraint feasibility (needed when continuous variables exist)
	// fData_p->constrSlack = IloNumArray2(*dataEnv, 2);
	// readArray<IloNumArray2> (fData_p->constrSlack, "data/constrSlack.dat");

	// read objective coefficients for x (binary) variables 
	fData_p->x = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->x, "data/x.dat");

	// read objective coefficients for y1 (integral) variables
	fData_p->y1 = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->y1, "data/y1.dat");

	// read objective coefficients for y2 (continuous) variables
	fData_p->y2 = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->y2, "data/y2.dat");

	// read rhs b
	fData_p->b = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->b, "data/rhs.dat");

	IloInt dimX = fData_p->x[0].getSize();   // state variables - binary
	IloInt dimZ = dimX;						 // copy of state variables
	IloInt dimY1 = fData_p->y1[0].getSize(); // local integral variables
	IloInt dimY2 = fData_p->y2[0].getSize(); // local continuous variables
	IloInt numRows = fData_p->b[0].getSize();// number of constraints in constr1

	// read matrix A
	fData_p->A = IloNumArray3(*dataEnv, numStage);
	IloInt dimA[2] = {numRows, dimX};
	sparse2full3d(fData_p->A, dimA, "data/A.dat");
	cout << "A matrix complete." << endl;
	// cin.get();
	
	// read matrix B
	fData_p->B = IloNumArray3(*dataEnv, numStage);
	IloInt dimB[2] = {numRows, dimZ};
	sparse2full3d(fData_p->B, dimB, "data/B.dat");	
	cout << "B matrix complete." << endl;
	// cin.get();

	// read matrix W1
	fData_p->W1 = IloNumArray3(*dataEnv, numStage);
	IloInt dimW1[2] = {numRows, dimY1};
	sparse2full3d(fData_p->W1, dimW1, "data/W1.dat");
	cout << "W1 matrix complete." << endl;
	// cin.get();

	// read matrix W2
	fData_p->W2 = IloNumArray3(*dataEnv, numStage);
	IloInt dimW2[2] = {numRows, dimY2};
	sparse2full3d(fData_p->W2, dimW2, "data/W2.dat");
	cout << "W2 matrix complete." << endl;
	// cin.get();

	// read matrix S
	fData_p->S = IloNumArray2(*dataEnv);
	IloInt dimS[2] = {numRows, numRows};
	sparse2full2d(fData_p->S, dimS, "data/S.dat");
	cout << "S matrix complete." << endl;
	// cin.get();

	// read uncertain indices
	fData_p->uncertainData = IloIntArray(*dataEnv, 8);
	readArray<IloIntArray> (fData_p->uncertainData, "data/uncertainData.dat");
	cout << "uncertainData index complete." << endl;
	// cin.get();

	IloIntArray index = fData_p->uncertainData;
	if ( index[0] )
	{
		fData_p->xScen = IloNumArray3(*dataEnv, numStage);
		cout << "load xScen..." << endl;
		readArray<IloNumArray3> (fData_p->xScen, "data/xScen.dat");
	}
	if ( index[1] )
	{
		fData_p->y1Scen = IloNumArray3(*dataEnv, numStage);
		cout << "load y1Scen..." << endl;
		readArray<IloNumArray3> (fData_p->y1Scen, "data/y1Scen.dat");
	}
	if ( index[2] )
	{
		fData_p->y2Scen = IloNumArray3(*dataEnv, numStage);
		cout << "load y2Scen..." << endl;
		readArray<IloNumArray3> (fData_p->y2Scen, "data/y2Scen.dat");
	}
	if ( index[3] )
	{
		fData_p->AScen = IloNumArray4(*dataEnv, numStage);
		cout << "load AScen..." << endl;
		char fileName[100];
		IloInt dimA [2] = {numRows, dimX};
		for ( int t = 0; t < fData_p->numStage; ++t )
		{
			fData_p->AScen[t] = IloNumArray3(*dataEnv, fData_p->numScen[t]);
			sprintf(fileName, "data/AScen_%d.dat", t);
			sparse2full3d(fData_p->AScen[t], dimA, fileName);
			// cout << "period " << t << endl;
		}
	}
	if ( index[4] )
	{
		fData_p->BScen = IloNumArray4(*dataEnv, numStage);
		cout << "load BScen..." << endl;
		char fileName[100];
		IloInt dimB [2] = {numRows, dimZ};
		for ( int t = 0; t < fData_p->numStage; ++t )
		{
			fData_p->BScen[t] = IloNumArray3(*dataEnv, fData_p->numScen[t]);
			sprintf(fileName, "data/BScen_%d.dat", t);
			sparse2full3d(fData_p->BScen[t], dimB, fileName);
			// cout << "period " << t << endl;
		}
	}
	if ( index[5] )
	{
		fData_p->W1Scen = IloNumArray4(*dataEnv, numStage);
		cout << "load W1Scen..." << endl;
		char fileName[100];
		IloInt dimW1 [2] = {numRows, dimY1};
		for ( int t = 0; t < fData_p->numStage; ++t )
		{
			fData_p->W1Scen[t] = IloNumArray3(*dataEnv, fData_p->numScen[t]);
			sprintf(fileName, "data/W1Scen_%d.dat", t);
			sparse2full3d(fData_p->W1Scen[t], dimW1, fileName);
			// cout << "period " << t << endl;
		}
	}
	if ( index[6] )
	{
		fData_p->W2Scen = IloNumArray4(*dataEnv, numStage);
		cout << "load W2Scen..." << endl;
		char fileName[100];
		IloInt dimW2 [2] = {numRows, dimY2};
		for ( int t = 0; t < fData_p->numStage; ++t )
		{
			fData_p->W2Scen[t] = IloNumArray3(*dataEnv, fData_p->numScen[t]);
			sprintf(fileName, "data/W1Scen_%d.dat", t);
			sparse2full3d(fData_p->W2Scen[t], dimW2, fileName);
			// cout << "period " << t << endl;
		}
	}
	if ( index[7] )
	{
		fData_p->bScen = IloNumArray3(*dataEnv, numStage);
		cout << "load bScen..." << endl;
		readArray<IloNumArray3> (fData_p->bScen, "data/bScen.dat");
	}
} // End of readData

void buildModel (model * models, formatData * fData_p, const bool LP)
{
	cout << "Start to build one model for each stage..." << endl;

	int t, i, j, k;
	char varName[100];
	IloInt numStage = fData_p->numStage;

	// extract decision variable dimensions
	IloInt dimX = fData_p->x[0].getSize();   // state variables - binary
	IloInt dimZ = dimX;						 // copy of state variables
	IloInt dimY1 = fData_p->y1[0].getSize(); // current stage integral variables
	IloInt dimY2 = fData_p->y2[0].getSize(); // current stage continuous variables
	IloInt numRows = fData_p->b[0].getSize();// number of constraints in constr1

	for ( t = 0; t < numStage; ++t )
	{
		cout << "start building model at stage " << t << endl;

		// create cplex environment and initialize model
		models[t].env = IloEnv();
		IloEnv currentEnv = models[t].env;
		models[t].mod = IloModel(currentEnv);

		// create variables and add them to model
		// x variable
		if ( LP )
			models[t].x = IloNumVarArray(currentEnv, dimX, 0.0, 1.0, ILOFLOAT);
		else
			models[t].x = IloNumVarArray(currentEnv, dimX, 0.0, 1.0, ILOINT);
		for ( i = 0; i < dimX; ++i )
		{
			sprintf(varName, "x_%d", i+1);
			models[t].x[i].setName(varName);
		}
		models[t].mod.add(models[t].x);

		// z variable
		models[t].z = IloNumVarArray(currentEnv, dimZ, 0.0, 1.0, ILOFLOAT);
		for ( i = 0; i < dimZ; ++i )
		{
			sprintf(varName, "z_%d", i+1);
			models[t].z[i].setName(varName);
		}
		models[t].mod.add(models[t].z);

		// y1 y2 variables
		if ( LP )
			models[t].y1 = IloNumVarArray(currentEnv, dimY1, 0.0, 100, ILOFLOAT);
		else
			models[t].y1 = IloNumVarArray(currentEnv, dimY1, 0.0, 100, ILOINT);
		models[t].y2 = IloNumVarArray(currentEnv, dimY2, 0.0, 200, ILOFLOAT);
		for ( i = 0; i < dimY1; ++i )
		{
			sprintf(varName, "y1_%d", i+1);
			models[t].y1[i].setName(varName);
		}
		for ( i = 0; i < dimY2; ++i )
		{
			sprintf(varName, "y2_%d", i+1);
			models[t].y2[i].setName(varName);
		}
		models[t].mod.add(models[t].y1);
		models[t].mod.add(models[t].y2);

		// theta variable
		models[t].theta = IloNumVar(currentEnv, fData_p->thetaLB[t], 1e8);
		sprintf(varName, "theta_%d", t+1);
		models[t].theta.setName(varName);
		models[t].mod.add(models[t].theta);

		// slack variables
		models[t].s = IloNumVarArray(currentEnv, numRows, 0.0, 2000, ILOFLOAT);
		for ( i = 0; i < numRows; ++i )
		{
			sprintf(varName, "s_%d", i+1);
			models[t].s[i].setName(varName);
		}
		models[t].mod.add(models[t].s);

		// slack variables 2 （needed when continuous state variables exist）
		// models[t].s2 = IloNumVarArray(currentEnv, fData_p->constrSlack[0], fData_p->constrSlack[1], ILOFLOAT);
		// for ( i = 0; i < numRows; ++i )
		// {
		// 	sprintf(varName, "s2_%d", i+1);
		// 	models[t].s2[i].setName(varName);
		// }
		// models[t].mod.add(models[t].s2);

		// cout << "All decision variables added" << endl;

		// create objective function
		IloExpr objExpr(currentEnv);
		objExpr = IloScalProd(models[t].x, fData_p->x[t]);
		objExpr += IloScalProd(models[t].y1, fData_p->y1[t]);
		objExpr += IloScalProd(models[t].y2, fData_p->y2[t]);
		if ( t < numStage - 1 )
			objExpr += models[t].theta;
		char objName[100];
		sprintf(objName, "objective_%d", t);
		models[t].obj = IloObjective(currentEnv, objExpr, IloObjective::Minimize, objName);
		models[t].mod.add(models[t].obj);
		objExpr.end();
		// cout << "objective added." << endl;

		// create constraints
		// Add constraints A_tx_t + B_tz_t + W1_ty1_t + W2_ty2_t + S s_t == b_t
		models[t].constr1 = IloRangeArray(currentEnv);
		for (i = 0; i < numRows; ++i)
		{
			IloExpr expr(currentEnv);
			expr = IloScalProd(models[t].x, fData_p->A[t][i]);
			expr += IloScalProd(models[t].z, fData_p->B[t][i]);
			expr += IloScalProd(models[t].y1, fData_p->W1[t][i]);
			expr += IloScalProd(models[t].y2, fData_p->W2[t][i]);
			expr += IloScalProd(models[t].s, fData_p->S[i]);
			// expr += models[t].s2[i];
			models[t].constr1.add(expr == fData_p->b[t][i]);
			expr.end();
		} // end of rows
		models[t].mod.add(models[t].constr1);
		// cout << "Constraints added." << endl;

		// Add constraints z_t = x_{t-1}, rhs initialized as 0
		models[t].constr2 = IloRangeArray(currentEnv);
		for (i = 0; i < dimZ; ++i)
			models[t].constr2.add(models[t].z[i] == fData_p->initState[i]);
		models[t].mod.add(models[t].constr2);
		// cout << "Constraints z_t = x_{t-1} added." << endl;

		// Initialize cut constraints
		models[t].cuts = IloRangeArray(currentEnv);

		// create cplex algorithms
		models[t].cplex = IloCplex(models[t].mod);

		// set model algorithm
		models[t].cplex.setParam(IloCplex::EpGap, 0.01);

		// models[t].cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);

		// set model solve output
		models[t].cplex.setOut(models[t].env.getNullStream());

		// set model warning message
		models[t].cplex.setWarning(models[t].env.getNullStream());

		// write model to lp file
		char fileName[100];
		sprintf(fileName, "m%d.lp", t);
		models[t].cplex.exportModel(fileName);

	} // End of t-loop for model construction
} // End of buildModel

void getSamplePaths (forwardPath & samplePaths, formatData * fData_p) 
	// const bool sampling)
{
	// cout << "Sampling forward paths from scenarios..." << endl;
	IloEnv env = fData_p->dataEnv;
	IloInt numStage = fData_p->numStage;
	IloInt numPaths = fData_p->numFWsample;

	// IloInt numPaths = 1;
	// if (sampling)
		// numPaths = fData_p->numFWsample;
	// else
		// for ( int t = 0; t < numStage; ++t )
			// numPaths *= numScen[t]
	IloIntArray index = fData_p->uncertainData;

	if ( index[0] )
		samplePaths.x = IloNumArray3(env);
	if ( index[1] )
		samplePaths.y1 = IloNumArray3(env);
	if ( index[2] )
		samplePaths.y2 = IloNumArray3(env);
	if ( index[3] )
		samplePaths.A = IloNumArray4(env);
	if ( index[4] )
		samplePaths.B = IloNumArray4(env);
	if ( index[5] )
		samplePaths.W1 = IloNumArray4(env);
	if ( index[6] )
		samplePaths.W2 = IloNumArray4(env);
	if ( index[7] )
		samplePaths.b = IloNumArray3(env);


	for (int p = 0; p < numPaths; ++p)
	{	
		// initialize each sample path
		if ( index[0] )
		{
			samplePaths.x.add(IloNumArray2(env));
			samplePaths.x[p].add(fData_p->xScen[0][0]);
		}
		if ( index[1] )
		{
			samplePaths.y1.add(IloNumArray2(env));
			samplePaths.y1[p].add(fData_p->y1Scen[0][0]);
		}
		if ( index[2] )
		{
			samplePaths.y2[p].add(IloNumArray2(env));
			samplePaths.y2[p].add(fData_p->y2Scen[0][0]);
		}
		if ( index[3] )
		{
			samplePaths.A.add(IloNumArray3(env));
			samplePaths.A[p].add(fData_p->AScen[0][0]);
		}
		if ( index[4] )
		{
			samplePaths.B.add(IloNumArray3(env));
			samplePaths.B[p].add(fData_p->BScen[0][0]);
		}
		if ( index[5] )
		{
			samplePaths.W1.add(IloNumArray3(env));
			samplePaths.W1[p].add(fData_p->W1Scen[0][0]);
		}
		if ( index[6] )
		{
			samplePaths.W2.add(IloNumArray3(env));
			samplePaths.W2[p].add(fData_p->W2Scen[0][0]);
		}
		if ( index[7] )
		{
			samplePaths.b.add(IloNumArray2(env));
			samplePaths.b[p].add(fData_p->bScen[0][0]);
		}

		// draw sample from each stage
		for (int t = 1; t < numStage; ++t)
		{
			double pdf = 1.0/fData_p->numScen[t];
			double U = genrand64_real1(); // generate a random number between 0 and 1
			// cout << "scenario chosen: " << int(U/pdf) << endl;
			int chosen = int(U/pdf);
			// int chosen = p;
			// if ( sampling )
				// chosen = int(U/pdf);
			if ( index[0] )
				samplePaths.x[p].add(fData_p->xScen[t][chosen]);
			if ( index[1] )
				samplePaths.y1[p].add(fData_p->y1Scen[t][chosen]);
			if ( index[2] )
				samplePaths.y2[p].add(fData_p->y2Scen[t][chosen]);
			if ( index[3] )
				samplePaths.A[p].add(fData_p->AScen[t][chosen]);
			if ( index[4] )
				samplePaths.B[p].add(fData_p->BScen[t][chosen]);
			if ( index[5] )
				samplePaths.W1[p].add(fData_p->W1Scen[t][chosen]);
			if ( index[6] )
				samplePaths.W2[p].add(fData_p->W2Scen[t][chosen]);
			if ( index[7] )
				samplePaths.b[p].add(fData_p->bScen[t][chosen]);
		} // End of stage for-loop
	} // End of sample paths for loop
} // End of getSamplePaths

void forward (model * models, formatData * fData_p,
	const forwardPath samplePaths, IloNumArray3 & candidateSol,
	IloNumArray & ub_c, IloNumArray & ub_l, IloNumArray & ub_r, const bool LP)
{
	cout << "Start the forward process..." << endl;

	int p, t, i;

	IloInt numRows = models[0].constr1.getSize();
	IloInt sampleSize = fData_p->numFWsample;
	IloIntArray index = fData_p->uncertainData;
	IloAlgorithm::Status solStatus;

	// create an array to record the objective function value for each sample path
	IloNumArray sampleObj(fData_p->dataEnv, sampleSize);
	// for ( p = 0; p < sampleSize; ++p) sampleObj[p] = 0.0;

	// Find candidate solutions for each sample path
	for ( p = 0; p < sampleSize; ++p )
	{
		// cout << "================================" << endl;
		// cout << "Compute solution for sample path p =  " << p << endl;
		candidateSol.add(IloNumArray2(fData_p->dataEnv));

		for ( t = 0; t < fData_p->numStage; ++t )
		{
			// cout << "Stage t = " << t << endl;
			if ( t > 0 )
			{
				if ( index[0] )
					models[t].obj.setLinearCoefs(models[t].x, samplePaths.x[p][t]);
				if ( index[1] )
					models[t].obj.setLinearCoefs(models[t].y1, samplePaths.y1[p][t]);
				if ( index[2] )
					models[t].obj.setLinearCoefs(models[t].y2, samplePaths.y2[p][t]);
				for ( i = 0; i < numRows; ++i )
				{
					if ( index[3] )
						models[t].constr1[i].setLinearCoefs(models[t].x, samplePaths.A[p][t][i]);
					if ( index[4] )
						models[t].constr1[i].setLinearCoefs(models[t].z, samplePaths.B[p][t][i]);
					if ( index[5] )
						models[t].constr1[i].setLinearCoefs(models[t].y1, samplePaths.W1[p][t][i]);
					if ( index[6] )
						models[t].constr1[i].setLinearCoefs(models[t].y2, samplePaths.W2[p][t][i]);
				}
				if ( index[7] )
					models[t].constr1.setBounds(samplePaths.b[p][t], samplePaths.b[p][t]);

				// update the state variables z_t = x_{t-1}
				models[t].constr2.setBounds(candidateSol[p][t-1], candidateSol[p][t-1]);
				if ( t == 1 )
				{
					for ( i = 0; i < models[t].constr2.getSize(); ++i )
						models[t].constr2[i].setBounds(0.0, 0.0);
				}

				// char fileName[100];
				// sprintf(fileName, "model_%d.lp", t);
				// models[t].cplex.exportModel(fileName);
			} // End of update problem models[t]

			// solve the current MIP/LP
			if ( models[t].cplex.solve() )
			{
				//get solution status
				solStatus = models[t].cplex.getStatus();
				// cout << "Solution status: " << solStatus << endl;

				if ( solStatus == IloAlgorithm::Optimal )
				{
					// record the solution
					IloNumArray vals(fData_p->dataEnv);
					models[t].cplex.getValues(vals, models[t].x);
					
					if ( LP )
					{
						for ( i = 0; i < vals.getSize(); ++i )
						{
							if (abs(vals[i]) < 2 * 1e-4)
							{
								// cout << k1 << ": " << vals[k1] << endl;
								vals[i] = 0.0;
								// cout << models[t].x[k1] << ": " << models[t].cplex.getValue(models[t].x[k1]) << endl; 
							}
						}
					}
					else
					{
						for ( i = 0; i < vals.getSize(); ++i )
							vals[i] = round(vals[i]);
					}
					// cout << "stage objective: " << models[t].cplex.getObjValue() << endl;
					// cout << vals << endl;
					// cin.get();
					candidateSol[p].add(vals);

					// IloNumArray unit(fData_p->dataEnv, vals.getSize()/3);
					// int j, k;
					// int FC = 6;
					// int ODI = 12;
					// int numBin [2] = {6, 9};
					// int dim = vals.getSize();
					// IloNumArray2 T(fData_p->dataEnv, FC);
					// for ( i = 0; i < 6; ++i )
					// {
					// 	T[i] = IloNumArray(fData_p->dataEnv, dim/3/ODI);
					// 	if ( i < 2 )
					// 	{
					// 		for ( j = 0; j < numBin[0]; ++j )
					// 			T[i][i*numBin[0] + j] = pow(2,j);
					// 	}
					// 	else
					// 	{
					// 		for ( j = 0; j < numBin[1]; ++j )
					// 			T[i][2 * numBin[0] + (i-2)*numBin[1] + j] = pow(2,j);
					// 	}
					// }

					// // cout << T << endl;
					// // cin.get();

					// for ( i = 0; i < ODI; ++i )
					// {
					// 	IloNumArray2 P(fData_p->dataEnv, FC);
					// 	IloNumArray2 B(fData_p->dataEnv, FC);
					// 	IloNumArray2 C(fData_p->dataEnv, FC);

					// 	for ( j = 0; j < FC; ++j )
					// 	{
					// 		B[j] = IloNumArray(fData_p->dataEnv, dim/3/ODI);
					// 		C[j] = IloNumArray(fData_p->dataEnv, dim/3/ODI);
					// 		P[j] = IloNumArray(fData_p->dataEnv, dim/3/ODI);

					// 		if ( j < 2 )
					// 		{
					// 			for ( k = 0; k < numBin[0]; ++k )
					// 			{
					// 				B[j][j*numBin[0] + k] = vals[i*dim/3/ODI + j*numBin[0] + k];
					// 				C[j][j*numBin[0] + k] = vals[dim/3 + i*dim/3/ODI + j*numBin[0] + k];
					// 				P[j][j*numBin[0] + k] = vals[dim/3*2 + i*dim/3/ODI + j*numBin[0] + k];
					// 			}
					// 		}
					// 		else
					// 		{
					// 			for ( k = 0; k < numBin[1]; ++k )
					// 			{
					// 				B[j][j*numBin[0] + k] = vals[i*dim/3/ODI + 2*numBin[0] + (j-2)*numBin[1] + k];
					// 				C[j][j*numBin[0] + k] = vals[dim/3 + i*dim/3/ODI + 2*numBin[0] + (j-2)*numBin[1] + k];
					// 				P[j][j*numBin[0] + k] = vals[dim/3*2 + i*dim/3/ODI + 2*numBin[0] + (j-2)*numBin[1] + k];
					// 			}
					// 		}

					// 		// cout << "B" << i << "," << j << ": " << IloScalProd(T[j], B[j]) << " ||  ";
					// 		// cout << "C" << i << "," << j << ": " << IloScalProd(T[j], C[j]) << " ||  ";
					// 		// cout << "P" << i << "," << j << ": " << IloScalProd(T[j], P[j]) << endl;
					// 	}
					// }

					
					// IloNumArray val2(fData_p->dataEnv);
					// models[t].cplex.getValues(val2, models[t].y2);
					// cout << val2 << endl;
					// cin.get();
					
					// Update objective function value of current sample path:
					// sampleObj[p] += models[t].ObjValue - theta_{t+1}^*
					if ( t < fData_p->numStage - 1 )
					{
						IloNum costToGo = models[t].cplex.getValue(models[t].theta);
						sampleObj[p] = sampleObj[p] + models[t].cplex.getObjValue() - costToGo;
					}
					else
						sampleObj[p] = sampleObj[p] + models[t].cplex.getObjValue();
				}
				else
				{
					cout << "Solution status is not optimal... " << endl;
					cout << "Sample " << p << ", Stage " << t << endl;
					throw ("Not Optimal...");
				}
			}
			else
			{
				cout << "Solution status is infeasible..." << endl;
				cout << "Sample " << p << ", Stage " << t << endl;
				models[t].cplex.exportModel("fw.lp");
				throw ("FW infeasible...");
			}

		} // End of loop over stages

	} //End of loop over sub-samples

	// calculate CI for upper bounds
	IloNum center = IloSum(sampleObj)/sampleSize;
	IloNum sumSq = 0.0;
	for ( p = 0; p < sampleSize; ++p)
		sumSq += pow( (sampleObj[p] - center), 2);
	IloNum halfLength = ZALPHA * sqrt(sumSq / (sampleSize - 1)) / sqrt(sampleSize);
	ub_c.add(center);
	ub_l.add(center - halfLength);
	ub_r.add(center + halfLength);
	// cout << sampleObj << endl;
	cout << "95\% CI for the upper bound: [ " << center - halfLength <<  ", " << center + halfLength << " ]." << endl;

	// free momory
	sampleObj.end();
} // End of forward pass to generate candidate solutions

void backward (model * models, formatData * fData_p, const IloNumArray3 candidateSol,
	IloNumArray & lb, const bool cutFlag[4], bool LP, const IloInt iter)
//	cutFlag = [bendersFlag, impvdBendersFlag, lagrangianFlag，integerFlag]
{
	cout << "================================" << endl;
	cout << "Start the backward process..." << endl;

	IloInt p, t, k, i;
	IloNum rhs;
	IloInt sampleSize = fData_p->numFWsample;
	IloInt numRows = models[0].constr1.getSize();
	IloInt dimX = models[0].constr2.getSize();
	IloIntArray index = fData_p->uncertainData;
	IloAlgorithm::Status solStatus;

	for ( t = fData_p->numStage - 1; t > 0; --t ) // for each stage
	{
		cout << "======================" << endl;
		cout << "Current stage t =  " << t << endl;

		for ( p = 0; p < sampleSize; ++p )  // 
		{
			cout << "Evaluating candidate solution: " << p << endl;
			cout << "\e[A";
			// cout << "Total number of scenarios at this stage: " << fData_p->numScen[t] << endl;
			// create arrays to store MIP, Lagrangian with optimal LP dual, and LP optimal values
			IloNumArray scenLPobj(fData_p->dataEnv);
			IloNumArray scenMIPobj(fData_p->dataEnv);
			IloNumArray scenMIPBestBound(fData_p->dataEnv);
			IloNumArray scenLGOPTobj(fData_p->dataEnv);
			IloNumArray scenLGobj(fData_p->dataEnv);
			IloNumArray dualVar(fData_p->dataEnv, dimX);
			IloNumArray LPdualAvg(fData_p->dataEnv, dimX);
			IloNumArray LGdualAvg(fData_p->dataEnv, dimX);

			for ( k = 0; k < fData_p->numScen[t]; ++k )  // for each scenario at stage t
			{
				// update the problem
				if ( index[0] )
					models[t].obj.setLinearCoefs(models[t].x, fData_p->xScen[t][k]);
				if ( index[1] )
					models[t].obj.setLinearCoefs(models[t].y1, fData_p->y1Scen[t][k]);
				if ( index[2] )
					models[t].obj.setLinearCoefs(models[t].y2, fData_p->y2Scen[t][k]);
				for ( i = 0; i < numRows; ++i )
				{
					if ( index[3] )
						models[t].constr1[i].setLinearCoefs(models[t].x, fData_p->AScen[t][k][i]);
					if ( index[4] )
						models[t].constr1[i].setLinearCoefs(models[t].z, fData_p->BScen[t][k][i]);
					if ( index[5] )
						models[t].constr1[i].setLinearCoefs(models[t].y1, fData_p->W1Scen[t][k][i]);
					if ( index[6] )
						models[t].constr1[i].setLinearCoefs(models[t].y2, fData_p->W2Scen[t][k][i]);
				}
				if ( index[7] )
					models[t].constr1.setBounds(fData_p->bScen[t][k], fData_p->bScen[t][k]);

				// update the state variables z_t = x_{t-1}
				models[t].constr2.setBounds( candidateSol[p][t-1], candidateSol[p][t-1] );
				if ( t == 1 )
				{
					for ( i = 0; i < models[t].constr2.getSize(); ++i )
						models[t].constr2[i].setBounds(0.0, 0.0);
				}

				// char fileName[100];
				// sprintf(fileName, "model_%d_%d_%d.lp", t, iter, k);
				// models[t].cplex.exportModel(fileName);
				
				// solve nodal MIP/LP
				// If MIP, objective value can be used in integer cut or serve as an ub for LGobj
				if ( models[t].cplex.solve() ) // feasible
				{	
					// get solution status
					solStatus = models[t].cplex.getStatus();

					if ( solStatus == IloAlgorithm::Optimal )
					{
						// record objective function value
						if ( LP )
						{
							scenLPobj.add(models[t].cplex.getObjValue());
							models[t].cplex.getDuals(dualVar, models[t].constr2);
							for (int k1 = 0; k1 < dualVar.getSize(); ++k1 )
							{
								if (abs(dualVar[k1]) < 2 * 1e-4)
									dualVar[k1] = 0.0;
							}
							// cout << "dualVar: " << dualVar << endl;
							for (i = 0; i < dimX; ++i)
								LPdualAvg[i] += dualVar[i] / fData_p->numScen[t];
						}
						else
						{
							scenMIPobj.add(models[t].cplex.getObjValue()); // incumbent, upper bound
							scenMIPBestBound.add(models[t].cplex.getBestObjValue()); // best lower bound
							// cout << "mip objective value: " << models[t].cplex.getObjValue() << endl;
							// models[t].cplex.exportModel("IP.lp");
						}


						// // record the solution
						// IloNumArray vals(fData_p->dataEnv);
						// models[t].cplex.getValues(vals, models[t].x);
						// if ( ! LP )
						// 	for ( i = 0; i < vals.getSize(); ++i )
						// 		vals[i] = round(vals[i]);

						// int j, m;
						// int FC = 6;
						// int ODI = 12;
						// int numBin [2] = {6, 9};
						// int dim = vals.getSize();
						// IloNumArray2 T(fData_p->dataEnv, FC);
						// for ( i = 0; i < 6; ++i )
						// {
						// 	T[i] = IloNumArray(fData_p->dataEnv, dim/2/ODI);
						// 	if ( i < 2 )
						// 	{
						// 		for ( j = 0; j < numBin[0]; ++j )
						// 			T[i][i*numBin[0] + j] = pow(2,j);
						// 	}
						// 	else
						// 	{
						// 		for ( j = 0; j < numBin[1]; ++j )
						// 			T[i][2 * numBin[0] + (i-2)*numBin[1] + j] = pow(2,j);
						// 	}
						// }

						// // cout << T << endl;
						// // cin.get();

						// for ( i = 0; i < ODI; ++i )
						// {
						// 	// IloNumArray2 P(fData_p->dataEnv, FC);
						// 	IloNumArray2 B(fData_p->dataEnv, FC);
						// 	IloNumArray2 C(fData_p->dataEnv, FC);

						// 	for ( j = 0; j < FC; ++j )
						// 	{
						// 		B[j] = IloNumArray(fData_p->dataEnv, dim/2/ODI);
						// 		C[j] = IloNumArray(fData_p->dataEnv, dim/2/ODI);
						// 		if ( j < 2 )
						// 		{
						// 			for ( m = 0; m < numBin[0]; ++m )
						// 			{
						// 				// cout << dim/2 << " || " << i*dim/2/ODI + j*numBin[0] + m << endl;
						// 				B[j][j*numBin[0] + m] = vals[        i*dim/2/ODI + j*numBin[0] + m];
						// 				C[j][j*numBin[0] + m] = vals[dim/2 + i*dim/2/ODI + j*numBin[0] + m];
						// 			}
						// 		}
						// 		else
						// 		{
						// 			for ( m = 0; m < numBin[1]; ++m )
						// 			{
						// 				// cout << dim/2 << " || " << i*dim/2/ODI + 2*numBin[0] + (j-2)*numBin[1] + m << endl;
						// 				B[j][2*numBin[0] + (j-2)*numBin[1] + m] = vals[        i*dim/2/ODI + 2*numBin[0] + (j-2)*numBin[1] + m];
						// 				C[j][2*numBin[0] + (j-2)*numBin[1] + m] = vals[dim/2 + i*dim/2/ODI + 2*numBin[0] + (j-2)*numBin[1] + m];
						// 			}
						// 		}

						// 		if ( IloScalProd(T[j], B[j]) < IloScalProd(T[j], C[j]) )
						// 		{
						// 			cout << "B" << i << "," << j << ": " << IloScalProd(T[j], B[j]) << " ||  ";
						// 			cout << "C" << i << "," << j << ": " << IloScalProd(T[j], C[j]) << endl;
						// 			cout << T[j] << endl;
						// 			cout << B[j] << endl;
						// 			cout << C[j] << endl;
						// 		}
						// 	}
						// 	B.end();
						// 	C.end();
						// }
						// // cin.get();




					}
					else // not optimal
					{
						cout << "Solution status: " << solStatus << endl;
						throw ("Stage MIP/LP not optimal...");
					}
				}
				else // infeasible
				{
					cout << "Solution status: " << models[t].cplex.getStatus() << endl;
					models[t].cplex.exportModel("bw.lp");
					cout << models[t].cplex.isDualFeasible() << endl;
					cout << models[t].cplex.isPrimalFeasible() << endl;
					throw ("Stage MIP/LP has no solution...");
				}

				if ( (! LP) && ( cutFlag[0] + cutFlag[1] + cutFlag[2] ) ) // Benders + ImprovedBenders + Lagrangian
				{
					// solve models[t] LP relaxation
					IloConversion relaxVarX(models[t].env, models[t].x, ILOFLOAT);
					IloConversion relaxVarY(models[t].env, models[t].y1, ILOFLOAT);
					models[t].mod.add(relaxVarX);
					models[t].mod.add(relaxVarY);
					// cout << "Relaxation added." << endl;
					
					if ( models[t].cplex.solve() ) // LP relaxation feasible
					{
						// get solution status
						solStatus = models[t].cplex.getStatus();
						// cout << "LP relaxation: " << solStatus << endl;

						if ( solStatus == IloAlgorithm::Optimal )
						{
							// record objective function value
							scenLPobj.add(models[t].cplex.getObjValue());
							// record optimal dual variables
							models[t].cplex.getDuals(dualVar, models[t].constr2);
							for ( i = 0; i < dualVar.getSize(); ++i )
							{
								if (abs(dualVar[i]) < 2 * 1e-4)
									dualVar[i] = 0.0;
							}
							// cout << "dualVar: " << dualVar << endl;
							// cout << models[t].constr2 << endl;
							for (i = 0; i < dimX; ++i)
								LPdualAvg[i] += dualVar[i] / fData_p->numScen[t];
							
							// change model back to MIP
							relaxVarX.end();
							relaxVarY.end();

							// continue to solve Lagrangian if needed
							if ( cutFlag[1] + cutFlag[2] ) // ImprovedBenders + Lagrangian
							{
								// create a new model with models[t].env
								model modelLGD;
								modelLGD.env = models[t].env;
								modelLGD.mod = IloModel(modelLGD.env);

								// add variables and constr1
								modelLGD.mod.add(models[t].x);
								modelLGD.mod.add(models[t].z);
								modelLGD.mod.add(models[t].y1);
								modelLGD.mod.add(models[t].y2);
								modelLGD.mod.add(models[t].s);
								// modelLGD.mod.add(models[t].s2);
								modelLGD.mod.add(models[t].theta);
								modelLGD.mod.add(models[t].constr1);
								modelLGD.mod.add(models[t].cuts);

								// create new objective function with additional term -\pi_LP' z
								modelLGD.obj = IloMinimize(modelLGD.env);
								IloExpr expr = models[t].obj.getExpr();
								expr -= IloScalProd(dualVar, models[t].z);
								
								// int k1 = dualVar.getSize();
								// for ( int k2 = k1/3*2; k2 < k1; ++k2 )
								// 	expr -= dualVar[k2] * models[t].z[k2];
								
								modelLGD.obj.setExpr(expr);
								modelLGD.mod.add(modelLGD.obj);

								// create cplex algorithm for this model
								modelLGD.cplex = IloCplex(modelLGD.mod);
								// char fileName[100];
								// sprintf(fileName, "LGDmodel_%d.lp", t);
								// cplexLGD.exportModel(fileName);
								modelLGD.cplex.setOut(modelLGD.env.getNullStream());

								// solve IP to get improved Benders' cut
								if ( cutFlag[1] && modelLGD.cplex.solve() )
										scenLGobj.add(modelLGD.cplex.getBestObjValue());
								
								// update lagrangian multipliers
								if ( cutFlag[2] )
								{
									double val = LGsolve(modelLGD, dualVar, models[t].z, candidateSol[p][t], scenMIPobj[k]);
									scenLGOPTobj.add(val);
									for (i = 0; i < dimX; ++i)
										LGdualAvg[i] += dualVar[i] / fData_p->numScen[t];
									cout << "Lagrangian dual solved." << endl;
									// cout << "IP solution " << scenMIPobj << endl;
									// cout << "LGD solution " << scenLGOPTobj << endl;
									cout << "========================" << endl;
								}

								// release memory
								expr.end();
								modelLGD.obj.end();
								modelLGD.mod.end();
								modelLGD.cplex.end();
							}
						}
						else // not optimal
						{
							cout << "Solution status: " << solStatus << endl;
							throw ("LP not optimal...");
						}
					}
					else // infeasible / unbounded
					{
						cout << "Solution status: " << models[t].cplex.getStatus() << endl;
						models[t].cplex.exportModel("LP.lp");
						throw ("LP has no solution...");
					}
			
				}
			} // End of loop over all scenarios in stage t
			// cin.get();

			// construct and add different types of cuts
			if ( cutFlag[0] + cutFlag[1] ) // Benders + ImprovedBenders
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				
				// int k1 = LPdualAvg.getSize();
				// for ( int k2 = k1/3*2; k2 < k1; ++k2 )
				// 	expr -= LPdualAvg[k2] * models[t-1].x[k2];
				
				if ( t > 1)
					expr -= IloScalProd(LPdualAvg, models[t-1].x);
				// cout << expr << endl;
				// cout << LPdualAvg << endl;
				if ( cutFlag[0] )
				{
					// cout << "test" << t << p << k << endl;
					// cout << LPdualAvg.getSize() << candidateSol[p][t].getSize() << endl;		
					// cout << IloScalProd(LPdualAvg, candidateSol[p][t]) << endl;
					// cout << IloSum(scenLPobj) / fData_p->numScen[t] << endl;
					
					// int k3 = 0;
					// for ( int k2 = k1/3*2; k2 < k1; ++k2 )
					// 	k3 += LPdualAvg[k2] * candidateSol[p][t-1][k2];
					// rhs = IloSum(scenLPobj) / fData_p->numScen[t] - k3;
					
					if ( t == 1 )
					{
						rhs = IloSum(scenLPobj) / fData_p->numScen[t];
						// cout << scenLPobj << endl;
					}
					else
						rhs = IloSum(scenLPobj) / fData_p->numScen[t] - IloScalProd(LPdualAvg, candidateSol[p][t-1]);
					
					// cout << rhs << endl;
					// cout << "Benders, " ;
				}
				else
				{	
					// int k3 = 0;
					// for ( int k2 = 0; k2 < k1/3*2; ++k2 )
					// 	k3 += LPdualAvg[k2] * candidateSol[p][t-1][k2];
					// rhs = IloSum(scenLGobj) / fData_p->numScen[t] + k3;
					rhs = IloSum(scenLGobj) / fData_p->numScen[t];
					// cout << "Improved Benders', ";
				}
				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				// cout << expr << endl;
				// cout << rhs << endl;
				expr.end();
			}

			if ( cutFlag[2] ) // Lagrangian
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				expr -= IloScalProd(LGdualAvg, models[t-1].x);
				rhs = IloSum(scenLGOPTobj) / fData_p->numScen[t];
				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				// cout << "Lagrangian, ";
				expr.end();
			}

			// construct and add integer L-shaped cut
			if ( cutFlag[3] )
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				IloNum mipObjAvg = IloSum(scenMIPBestBound) / fData_p->numScen[t];
				if ( t > 1 )
				{	
					for ( int j = 0; j < dimX; ++j )
					{
						if ( abs( candidateSol[p][t][j] ) < EPSILON )
							expr += (mipObjAvg - fData_p->thetaLB[t-1]) * models[t-1].x[j];
						else
							expr -= (mipObjAvg - fData_p->thetaLB[t-1]) * models[t-1].x[j];
					}
					rhs = fData_p->thetaLB[t-1] - (mipObjAvg - fData_p->thetaLB[t-1]) * (IloSum(candidateSol[p][t]) -1);
				}
				else
					rhs = mipObjAvg;

				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				// cout << expr << ">= " << rhs << endl;
				// cout << "Integer L-shaped, ";
				expr.end();
			}

			// cout << "cuts are added." << endl;

			// free momery
			scenLPobj.end();
			dualVar.end();
			LPdualAvg.end();
			if ( ! LP )
			{
				scenMIPobj.end();
				scenMIPBestBound.end();
				scenLGOPTobj.end();
				scenLGobj.end();
				LPdualAvg.end();
			}
		} // End of loop over unique candidate solutions
		// cin.get();
		// cout << models[t-1].cuts << endl;
	} // End of loop over stages

	cout << "================================" << endl;
	cout << "Solving the master problem: " << endl;

	// solve models[0] as a MIP/LP
	if ( models[0].cplex.solve() ) // feasible
	{
		// get solution status
		solStatus = models[0].cplex.getStatus();
		// cout << "solution status obtained." << endl;

		if ( solStatus == IloAlgorithm::Optimal )
		{
			// update lower bound lb as the optimal objective function value
			if ( LP )
				lb.add(models[0].cplex.getObjValue());
			else
				lb.add(models[0].cplex.getBestObjValue());
			
			// char fileName[100];
			// cout << t << endl;
			// sprintf(fileName, "model_%d_%d.lp", t, iter);
			// models[t].cplex.exportModel(fileName);

			// record the solution
			IloNumArray vals(fData_p->dataEnv);
			models[0].cplex.getValues(vals, models[0].x);
			if ( ! LP )
				for ( i = 0; i < vals.getSize(); ++i )
					vals[i] = round(vals[i]);
			cout << "master objective: " << models[0].cplex.getObjValue() << endl;
			// cout << vals << endl;
			// cin.get();

			// int j, k;
			// int FC = 6;
			// int ODI = 12;
			// int numBin [2] = {6, 9};
			// int dim = vals.getSize();
			// IloNumArray2 T(fData_p->dataEnv, FC);
			// for ( i = 0; i < 6; ++i )
			// {
			// 	T[i] = IloNumArray(fData_p->dataEnv, dim/2/ODI);
			// 	if ( i < 2 )
			// 	{
			// 		for ( j = 0; j < numBin[0]; ++j )
			// 			T[i][i*numBin[0] + j] = pow(2,j);
			// 	}
			// 	else
			// 	{
			// 		for ( j = 0; j < numBin[1]; ++j )
			// 			T[i][2 * numBin[0] + (i-2)*numBin[1] + j] = pow(2,j);
			// 	}
			// }

			// cout << T << endl;
			// // cin.get();

			// for ( i = 0; i < ODI; ++i )
			// {
			// 	// IloNumArray2 P(fData_p->dataEnv, FC);
			// 	IloNumArray2 B(fData_p->dataEnv, FC);
			// 	IloNumArray2 C(fData_p->dataEnv, FC);

			// 	for ( j = 0; j < FC; ++j )
			// 	{
			// 		B[j] = IloNumArray(fData_p->dataEnv, dim/2/ODI);
			// 		C[j] = IloNumArray(fData_p->dataEnv, dim/2/ODI);
			// 		if ( j < 2 )
			// 		{
			// 			for ( k = 0; k < numBin[0]; ++k )
			// 			{
			// 				B[j][j*numBin[0] + k] = vals[i*dim/2/ODI + j*numBin[0] + k];
			// 				C[j][j*numBin[0] + k] = vals[dim/2 + i*dim/2/ODI + j*numBin[0] + k];
			// 			}
			// 		}
			// 		else
			// 		{
			// 			for ( k = 0; k < numBin[1]; ++k )
			// 			{
			// 				B[j][j*numBin[0] + k] = vals[i*dim/2/ODI + 2*numBin[0] + (j-2)*numBin[1] + k];
			// 				C[j][j*numBin[0] + k] = vals[dim/2 + i*dim/2/ODI + 2*numBin[0] + (j-2)*numBin[1] + k];
			// 			}
			// 		}

			// 		cout << "B" << i << "," << j << ": " << IloScalProd(T[j], B[j]) << " ||  ";
			// 		cout << "C" << i << "," << j << ": " << IloScalProd(T[j], C[j]) << endl;
			// 	}
			// }
			// // cin.get();
		}
		else // not optimal
		{
			cout << "First-stage problem status: " << solStatus << endl;
			throw ("first-stage problem not optimal...");
		}
	}
	else // infeasible
	{
		cout << "First stage problem status: " << models[0].cplex.getStatus() << endl;
		models[0].cplex.exportModel("root.lp");
		throw ("Relaxed first-stage problem infeasible...");
	}
} //


double LGsolve (model & LGmodel, IloNumArray & dualVar, IloNumVarArray z, const IloNumArray state, const double ub )
{
	int dimZ = dualVar.getSize();
	int i, j;
	int numIter = 5*1e3;
	double normGrad, stepSize, objVal, gap;
	bool optimalFlag = 0;
	IloNumArray gradient(LGmodel.env, dimZ);

	for ( i = 0; i < numIter; ++i )
	{
		LGmodel.cplex.solve();
		objVal = LGmodel.cplex.getObjValue();
		double constant = IloScalProd(dualVar, state);
		// cout << objVal + temp << endl;
		// cout << ub;
		double diff = ub - objVal - constant;
		gap = abs(diff / ub);
		cout << "Current gap: " << gap << " || Absolute difference: " << diff << endl;
		// optimality check
		if ( gap < 1e-4 )
		{
			optimalFlag = 1;
			cout << "Lagrangian problem is solved to optimality." << endl;
			break;
		}
		else if ( i < numIter - 1 )
		{
			cout << "\e[A";
			stepSize = 1.0 / sqrt(i+1);
		
			// obtain and revise gradient
			LGmodel.cplex.getValues(gradient, z);
			for ( j = 0; j < dimZ; ++j)
				gradient[j] -= state[j];
			normGrad = IloScalProd(gradient, gradient);

			// update dualVar and objective function
			for ( j = 0; j < dimZ; ++j )
			{
				dualVar[j] -= stepSize / (normGrad + 1e-3) * gradient[j];
				LGmodel.obj.setLinearCoef(z[j], -dualVar[j]);
			}
		}
	}
	return objVal;
} // End of subgradient method

// bool checkCut(IloRangeArray cuts, IloRange newcut)
// {
// 	int k = cuts.getSize();
// 	for (int i = 0; i < k ; ++i)
// 	{
// 		IloExpr expr = cuts[i].getExpr();
// 		IloNum LB = cuts[i].getLB();
// 		IloExpr newExpr = newcut.getExpr();
// 		IloNum newLB = newcut.getLB();
// 		if ( (expr - newExpr == 0) && (abs(LB - newLB) < 1e-6) )
// 		{
// 			return 1;
// 			break;
// 		}
// 	}
// 	return 0;
// }

double avg ( vector<float> & v )
{
    double sumV = 0.0;
    for ( unsigned i=0; i < v.size(); i++)
        sumV += v[i];
            
    return sumV / v.size();
}
// End of mean funtion

double std_dev ( vector<float> & v )
{
        double sumSquare = 0.0;
        double var =0.0;
        
        double mean = avg(v);
        
        for ( unsigned j = 0; j < v.size(); ++j )
            sumSquare += pow((v[j] - mean),2);
       
        return var = sqrt(sumSquare / v.size());
}
// End of standard deviation funtion

void usage (char *progname)
{
	cerr << "Usage:  " << progname << " arg1 arg2 arg3 arg4 arg5 [arg6]" << endl;
	cerr << "At least two parameters must be specified." << endl;
	cerr << "arg1: 0 -- turn off Benders' cuts;" << endl;
	cerr << "      1 -- turn on Benders' cuts." << endl;
	cerr << "arg2: 0 -- turn off Strengthened Benders' cuts;" << endl;
	cerr << "      1 -- turn on Strengthened Benders' cuts." << endl;
	cerr << "arg3: 0 -- turn off Lagrangian cuts;" << endl;
	cerr << "      1 -- turn on Lagrangian cuts." << endl;
	cerr << "arg4: 0 -- turn off L-shaped cuts;" << endl;
	cerr << "      1 -- turn on L-shaped cuts." << endl;
	cerr << "arg5: 0 -- solve the integer problem;" << endl;
	cerr << "      1 -- solve the LP-relaxation." << endl;
	cerr << "arg6: [optional] used as the seed of the random number generator." << endl;
	cerr << "      If not provided, system will generate one automatically." << endl;
} // END usage
