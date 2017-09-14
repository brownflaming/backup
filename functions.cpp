#include <ilcplex/ilocplex.h>
#include <cmath>
#include <unordered_set>
#include <string>
#include <sstream>
#include <vector>
#include <cstddef>  // std::size_t
#include <numeric>
#include <set>
#include <limits>
// #include <random>

#include "global.h"
#include "functions.h"
#include "mt64.h"


using namespace std;
//using std::tr1::unordered_set;

void readData (formatData * fData_p)
{
	// create data environment
	fData_p->dataEnv = IloEnv();
	IloEnv * dataEnv = & fData_p->dataEnv;

	cout << "Start reading data from files... " << endl;

	// read number of stages
	readArray<IloInt> (fData_p->numStage, "data/numStage.dat");
	IloInt numStage = fData_p->numStage;
	// -----------------------------------------
	// cout << "Time horizion: " << fData_p->numFWsample << endl;
	// -----------------------------------------
	
	// read number of scenarios available at each stage
	fData_p->numScen = IloIntArray(*dataEnv, numStage);
	readArray<IloIntArray> (fData_p->numScen, "data/numScen.dat");
	// -----------------------------------------
	cout << "Number of scenarios at each stage: " << fData_p->numScen << endl;
	// -----------------------------------------
	
	// compute the total number of scenarios
	fData_p->totalScen = fData_p->numScen[0];
	for (int t = 1; t < numStage; ++t)
		fData_p->totalScen *= fData_p->numScen[t];
	// -----------------------------------------
	cout << "Total number of scenarios: " << fData_p->totalScen << endl;
	// -----------------------------------------

	// read number of samples drawn in the forward pass (number of candidate solns)
	readArray<IloInt> (fData_p->numFWsample, "data/numFWsample.dat");
	// -----------------------------------------
	cout << "Number of FW samples: " << fData_p->numFWsample << endl;
	// -----------------------------------------
	
	// read initial state of binary variables
	fData_p->initState = IloNumArray(*dataEnv, numStage);
	readArray<IloNumArray> (fData_p->initState, "data/initState.dat");

	// read prescribed lower bound for value function at each stage
	fData_p->valueLB = IloNumArray(*dataEnv, numStage);
	readArray<IloNumArray> (fData_p->valueLB, "data/valueLB.dat");

	// read objective coefficients for x (binary) variables 
	fData_p->xCoef = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->xCoef, "data/xCoef.dat");
	// -----------------------------------------
	// cout << fData_p->xCoef[0] << fData_p->xCoef[2] << endl;
	// -----------------------------------------

	// read objective coefficients for y1 (integral) variables
	fData_p->y1Coef = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->y1Coef, "data/y1Coef.dat");

	// read objective coefficients for y2 (continuous) variables
	fData_p->y2Coef = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->y2Coef, "data/y2Coef.dat");

	// read matrix A in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t >= b_t
	fData_p->Amatrix = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->Amatrix, "data/Amatrix.dat");
	// -----------------------------------------
	// cout << fData_p->Amatrix[0] << endl;
	// -----------------------------------------
	
	// read matrix B in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t >= b_t
	fData_p->Bmatrix = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->Bmatrix, "data/Bmatrix.dat");
	// -----------------------------------------
	// cout << "Bmatrix" << endl;
	// -----------------------------------------
	
	// read matrix W1 in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t >= b_t
	fData_p->W1matrix = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->W1matrix, "data/W1matrix.dat");
	// -----------------------------------------
	// cout << "W1matrix" << endl;
	// -----------------------------------------
	
	// read matrix W2 in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t >= b_t
	fData_p->W2matrix = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->W2matrix, "data/W2matrix.dat");
	// -----------------------------------------
	// cout << "W2matrix" << endl;
	// -----------------------------------------
	
	// read rhs b in constraint A_tx_t + W_ty_t + B_tz_t <= b_t
	fData_p->bRhs = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->bRhs, "data/bRhs.dat");
	// -----------------------------------------
	// cout << "rhs" << endl;
	// -----------------------------------------

	// read uncertain indices in b (every stage have the same uncertain parameter)
	fData_p->uncertainIndex = IloIntArray(*dataEnv, numStage);
	readArray<IloIntArray> (fData_p->uncertainIndex, "data/uncertainIndex.dat");
	// -----------------------------------------
	// cout << "index" << endl;
	// -----------------------------------------

	// read generated rhs scenarios at each stage
	// dim1: numStage; dim2: numScen[t]; dim3: uncertainIndex.getSize()
	fData_p->scenarios = IloNumArray3(*dataEnv, numStage);
	readArray<IloNumArray3> (fData_p->scenarios, "data/scenarios.dat");
	// -----------------------------------------
	// cout << "scenarios" << endl;
	// -----------------------------------------

	// -----------------------------------------
	// cout << "y2Scenarios" << endl;
	// -----------------------------------------

	return;
} // End of readData

void buildModel (Model * models, formatData * fData_p)
{
	cout << "Start to build one model for each stage..." << endl;

	int t, i, j, k;
	IloInt numStage = fData_p->numStage;

	// extract decision variable dimensions
	IloInt x_dim = fData_p->xCoef[0].getSize();   // state variables
	IloInt z_dim = x_dim;						  // copy of state variables
	IloInt y1_dim = fData_p->y1Coef[0].getSize(); // current stage integral variables
	IloInt y2_dim = fData_p->y2Coef[0].getSize(); // current stage continuous variables
	IloInt numRows = fData_p->bRhs[0].getSize();  // number of constraints
	// -------------------------------------------------------
	cout << "number of constraints: " << numRows << endl;
	cout << "x(z) dim:" << x_dim << endl;
	cout << "y1 dim:" << y1_dim << endl;
	cout << "y2 dim:" << y2_dim << endl;
	cout << "A dim: " << fData_p->Amatrix.getSize() << 
			"x" << fData_p->Amatrix[0].getSize() << endl;
	cout << "B dim: " << fData_p->Bmatrix.getSize() <<
			"x" << fData_p->Bmatrix[0].getSize() << endl;
	cout << "W1 dim: " << fData_p->W1matrix.getSize() << 
			"x" << fData_p->W1matrix[0].getSize() << endl;
	cout << "W2 dim: " << fData_p->W2matrix.getSize() << 
			"x" << fData_p->W2matrix[0].getSize() << endl;
	// -------------------------------------------------------	

	for (t = 0; t < numStage; ++t)
	{
		// create cplex environment and initialize model
		models[t].env = IloEnv();
		IloEnv currentEnv = models[t].env;
		models[t].mod = IloModel(currentEnv);

		// create variables and add them to model
		models[t].x = IloNumVarArray(currentEnv, x_dim, 0.0, 1.0, ILOINT);
		models[t].y1 = IloNumVarArray(currentEnv, y1_dim, 0.0, IloInfinity, ILOINT);
		models[t].y2 = IloNumVarArray(currentEnv, y2_dim, 0.0, IloInfinity, ILOFLOAT);
		models[t].z = IloNumVarArray(currentEnv, z_dim, 0.0, 1.0, ILOFLOAT);
		models[t].theta = IloNumVar(currentEnv, 0.0, IloInfinity);
		
		char varName[100];
		for ( i = 0; i < x_dim; ++i )
		{
			sprintf(varName, "x_%d", i+1);
			models[t].x[i].setName(varName);
		}
		
		for ( i = 0; i < y1_dim; ++i )
		{
			sprintf(varName, "y1_%d", i+1);
			models[t].y1[i].setName(varName);
		}
		
		for ( i = 0; i < y2_dim; ++i )
		{
			sprintf(varName, "y2_%d", i+1);
			models[t].y2[i].setName(varName);
		}

		for ( i = 0; i < z_dim; ++i )
		{
			sprintf(varName, "z_%d", i+1);
			models[t].z[i].setName(varName);
		}

		sprintf(varName, "theta");
		models[t].theta.setName(varName);

		models[t].mod.add(models[t].x);
		models[t].mod.add(models[t].y1);
		models[t].mod.add(models[t].y2);		
		models[t].mod.add(models[t].z);
		models[t].mod.add(models[t].theta);

		cout << "variables" << endl;

		// create objective function and add to model
		models[t].obj = IloObjective(currentEnv);
		IloExpr objExpr(currentEnv);
		objExpr = IloScalProd(models[t].x, fData_p->xCoef[t]) + 
					IloScalProd(models[t].y1, fData_p->y1Coef[t]) +
					IloScalProd(models[t].y2, fData_p->y2Coef[t]) +
					models[t].theta;
		models[t].obj.setExpr(objExpr);
		objExpr.end();
		models[t].obj.setSense(IloObjective::Minimize);
		char objName[100];
		sprintf(objName, "objective_%d", t);
		models[t].obj.setName(objName);
		models[t].mod.add(models[t].obj);

		cout << "objective" << endl;

		// create constraints
		// Add constraints A_tx_t + B_tz_t + W1_ty1_t + W2_ty2_t >= b_t
		models[t].constr1 = IloRangeArray(currentEnv);
		for (i = 0; i < numRows; ++i)
		{
			IloExpr expr(currentEnv);
			expr = IloScalProd(models[t].x, fData_p->Amatrix[i]) + 
					IloScalProd(models[t].y1, fData_p->W1matrix[i]) +
					IloScalProd(models[t].y2, fData_p->W2matrix[i]) +
					IloScalProd(models[t].z, fData_p->Bmatrix[i]);
			models[t].constr1.add(expr >= fData_p->bRhs[t][i]);
			expr.end();
		} // end of rows
		models[t].mod.add(models[t].constr1);

		cout << "constraints 1" << endl;

		// Add constraints z_t = x_{t-1}, rhs initialized as 0
		models[t].constr2 = IloRangeArray(currentEnv);
		for (i = 0; i < z_dim; ++i)
		{
			IloExpr expr(currentEnv);
			expr = models[t].z[i];
			if ( t == 0 ) // initialize the state variable for models[0]
				models[t].constr2.add(expr == fData_p->initState[i]);
			else  // all others left as zero for future update
				models[t].constr2.add(expr == 0);
			expr.end();
		}
		models[t].mod.add(models[t].constr2);

		cout << "constraints 2" << endl;

		// Initialize cut constraints
		models[t].cuts = IloRangeArray(currentEnv);

		// relaxation conversion
		models[t].relaxVarX = IloConversion(models[t].env, models[t].x, ILOFLOAT);
		models[t].relaxVarY = IloConversion(models[t].env, models[t].y1, ILOFLOAT);

		// create cplex algorithms
		models[t].cplex = IloCplex(models[t].mod);

		// models[t].cplex.setParam(IloCplex::PrePass, 2);
		// models[t].cplex.setParam(IloCplex::RepeatPresolve, 3); 
		models[t].cplex.setParam(IloCplex::Param::Advance, 0);
		models[t].cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, MIPTOL);
		// models[t].cplex.setParam(IloCplex::Param::RandomSeed, rand() % CPX_BIGINT);

		// set model algorithm
		// models[t].cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Barrier);

		// set model solve output
		models[t].cplex.setOut(models[t].env.getNullStream());

		// set model warning message
		models[t].cplex.setWarning(models[t].env.getNullStream());

		// write model to lp file
		// char fileName[100];
		// sprintf(fileName, "model_%d.lp", t);
		// models[t].cplex.exportModel(fileName);

	} // End of t-loop for model construction
	return;
} // End of function buildModel


void getSamplePaths (IloNumArray3 & samplePaths, IloIntArray2 & scenarioIndex, formatData * fData_p, const bool sampling)
{
	// cout << "Sampling forward paths from scenarios..." << endl;

	IloEnv * currentEnv = &(fData_p->dataEnv);
	IloInt numStage = fData_p->numStage;

	if (sampling)
	{
		IloInt numPaths = fData_p->numFWsample;
		IloInt exsitingSampleSize = scenarioIndex.getSize();

		for (int p = 0; p < numPaths; ++p)
		{
			// cout << "sample path p = " << p << endl;
			
			// initialize each sample path
			samplePaths.add(IloNumArray2(*currentEnv));
			samplePaths[p].add(fData_p->scenarios[0][0]);

			scenarioIndex.add(IloIntArray(*currentEnv));

			// coefSamplePaths.add(IloNumArray2(*currentEnv));
			// coefSamplePaths[p].add(fData_p->y2Scenarios[0][0]);
			
			// draw sample from each stage
			for (int t = 1; t < numStage; ++t)
			{
				double pdf = 1.0/fData_p->numScen[t];
				double U = genrand64_real1();
				int chosen = int(U/pdf);
				samplePaths[p].add(fData_p->scenarios[t][chosen]);
				// coefSamplePaths[p].add(fData_p->y2Scenarios[t][chosen]);
				scenarioIndex[exsitingSampleSize + p].add(chosen);
			} // End of stage for-loop
		} // End of sample paths for loop
	}
	else
	{
		fData_p->numFWsample = fData_p->totalScen;
		IloInt numPaths = fData_p->totalScen;
		IloInt exsitingSampleSize = scenarioIndex.getSize();

		for (int p = 0; p < numPaths; ++p )
		{
			samplePaths.add(IloNumArray2(*currentEnv));
			samplePaths[p].add(fData_p->scenarios[0][0]);

			scenarioIndex.add(IloIntArray(*currentEnv));
			scenarioIndex[exsitingSampleSize + p].add(0);

			// coefSamplePaths.add(IloNumArray2(*currentEnv));
			// coefSamplePaths[p].add(fData_p->y2Scenarios[0][0]);

			int f = fData_p->totalScen;
			for ( int t = 1; t < numStage; ++t )
			{
				f = f / fData_p->numScen[t];
				int chosen = int(p / f) % fData_p->numScen[t];
				samplePaths[p].add(fData_p->scenarios[t][chosen]);
				// coefSamplePaths[p].add(fData_p->y2Scenarios[t][chosen]);
				scenarioIndex[exsitingSampleSize + p].add(chosen);
			}
		}
	}
	// cout << scenarioIndex << endl;

	return;
} // End of getSamplePathss

void forward (Model * models, formatData * fData_p, const IloNumArray3 samplePaths,
			  IloNumArray3 & candidateSol, IloNumArray & ub_c, IloNumArray & ub_l, IloNumArray & ub_r)
{
	cout << "Start the forward process..." << endl;

	int p, t, i, k, index;
	IloInt constr1Size = models[0].constr1.getSize();
	IloInt constr2Size = models[0].constr2.getSize();
	IloInt sampleSize = fData_p->numFWsample;
	IloInt uncertainArraySize = fData_p->uncertainIndex.getSize();
	IloAlgorithm::Status solStatus;

	// create an array to record the objective function value for each sample path
	IloNumArray sampleObj(fData_p->dataEnv, sampleSize);
	for ( p = 0; p < sampleSize; ++p)
		sampleObj[p] = 0.0;

	// cout << "Number of FW sample paths: " << sampleSize << endl;
	// Find candidate solutions for each sample path
	for ( p = 0; p < sampleSize; ++p )
	{
		// cout << "================================" << endl;
		cout << "Compute solution for sample path p =  " << p << endl;
		if ( p != sampleSize )
			cout << "\e[A";

		for ( t = 0; t < fData_p->numStage; ++t )
		{
			candidateSol.add(IloNumArray2(fData_p->dataEnv));
			// update objective coefficients
			// models[t].obj.setLinearCoefs(models[t].y2, coefSamplePaths[p][t]);
			
			// update b_t with samplePaths[p][t]
			for ( i = 0; i < uncertainArraySize; ++i )
			{
				index = fData_p->uncertainIndex[i];
				models[t].constr1[index].setLB(samplePaths[p][t][i]);
			}

			// update the state variables z_t = x_{t-1}
			if ( t > 0 )
			{
				int nsol = candidateSol[t-1].getSize();
				//for ( k = 0; k < NcandidateSol; ++k )
				//{
				// int select = rand() % NcandidateSol;
				// cout << "Available solutions: " << NcandidateSol << " selected: " << select + 1 << endl;
				models[t].constr2.setBounds(candidateSol[t-1][nsol-1], candidateSol[t-1][nsol-1]);

				//}	
			}
			// End of update problem models[t]

			// models[t].cplex.exportModel("FW.lp");


			// solve the current MIP model
			if ( models[t].cplex.solve() )
			{
				//get solution status
				solStatus = models[t].cplex.getStatus();
				// cout << "Solution status: " << solStatus << endl;

				if ( solStatus == IloAlgorithm::Optimal )
				{
					IloNumArray vals(fData_p->dataEnv);
					if ( (models[t].cplex.isMIP()) && (models[t].cplex.getParam(IloCplex::Param::MIP::Pool::Intensity) == 4) )
					{
						int nsol = models[t].cplex.getSolnPoolNsolns();
						cout << "number of solutions in the pool: " << nsol << endl;
						set<int> opt;
						double best = std::numeric_limits<double>::infinity();
						// cin.get();

						for ( i = 0; i < nsol; ++i)
						{
							double z = models[t].cplex.getObjValue(i);
							cout << "objective function value of current solution: " << z << endl;
							if ( z < best - 1e-2 )
							{
								best = z;
								opt.clear();
								opt.insert(i);
							}
							else if ( z < best + 1e-2 )
								opt.insert(i);
						}
						// cin.get();

						for ( k = 0; k < nsol; ++k )
						{
							if ( opt.find(k) == opt.end() )
							{
								models[t].cplex.getValues(vals, models[t].x, k);
								for ( i = 0; i < vals.getSize(); ++i )
									vals[i] = round(vals[i]);
								candidateSol[t].add(vals);
							}
						}

						for ( auto k = opt.begin(); k != opt.end(); ++k )
						{
							models[t].cplex.getValues(vals, models[t].x, *k);
							for ( i = 0; i < vals.getSize(); ++i )
								vals[i] = round(vals[i]);
							
							cout << "optimal value: " << models[t].cplex.getObjValue(*k) << endl;
							// cout << vals << endl;

							candidateSol[t].add(vals);
						}
					}
					else
					{
						// record the solution
						models[t].cplex.getValues(vals, models[t].x);
						// cout << vals << endl;
						// cout << candidateSol.getSize() << endl;
						// cout << candidateSol << endl;
						candidateSol[t].add(vals);
					}

					IloNum costToGo = models[t].cplex.getValue(models[t].theta);
					// cin.get();

					// Update objective function value of current sample path:
					// sampleObj[p] += models[t].ObjValue - theta_{t+1}^*
					sampleObj[p] = sampleObj[p] + models[t].cplex.getObjValue() - costToGo;
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
				char fileName[100];
				sprintf(fileName, "inf_model_%d.lp", t);
				models[t].cplex.exportModel(fileName);
				cout << "Solution status is infeasible..." << endl;
				cout << "Sample " << p << ", Stage " << t << endl;
				throw ("Infeasible...");
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
	// cout << "Objective value for forward sample paths: " << sampleObj << endl;
	/*
	cout << ub_l << endl;
	cout << ub_c << endl;
	cout << ub_r << endl;
	*/
	cout << "95\% CI for the upper bound: [ " << center - halfLength <<  ", " << center + halfLength << " ]." << endl;

	// free momory
	sampleObj.end();
} // End of forward pass to generate candidate solutions

void backward (Model * models, formatData * fData_p, const IloNumArray3 candidateSol, IloNumArray & lb,
		const bool bendersFlag,	const bool impvdBendersFlag,
		const bool integerFlag,	const bool lagrangianFlag)
{
	// cout << "================================" << endl;
	cout << "Start the backward process...." << endl;

	int p, t, k, i, j, index;
	IloNum rhs;
	IloInt sampleSize = fData_p->numFWsample;
	IloInt uncertainArraySize = fData_p->uncertainIndex.getSize();
	IloInt constr1Size = models[0].constr1.getSize();
	IloInt constr2Size = models[0].constr2.getSize();
	IloAlgorithm::Status solStatus;

	for ( t = fData_p->numStage - 1; t > 0; --t ) // for each stage
	{
		// cout << "======================" << endl;
		// cout << "Current stage t =  " << t << endl;

		// Find out unique candidate solutions at current stage among candidateSol[p][t] for p = 1,..., numFWsample
		// Create a 2-d array to store candidateSol[p][t] for p = 1,..., numFWsample
		// for ( p = 0; p < fData_p->numFWsample; ++p )
		// 	uniqueSet.insert(toString(candidateSol[p][t-1]));

		// // convert uniqueSet(unordered_set) to uniqueSol(2-d vector)
		// for ( auto it = uniqueSet.begin(); it != uniqueSet.end(); ++it )
		// 	uniqueSol.push_back(getDigit(*it));

		// for ( auto it = uniqueSol.begin(); it != uniqueSol.end(); ++it )  // 
		// for ( p = 0; p < fData_p->numFWsample; ++p )
		int nsol = candidateSol[t-1].getSize();
		// cout << "number of solutions: " <<  NcandidateSol << endl;
		for ( j = 0; j < nsol; ++j )
		{
			// cout << "Total number of scenarios at this stage: " << fData_p->numScen[t] << endl;
			// create arrays to store optimal values of original MIP, Lagrangian with optimal LP dual, Lagrangian, and LP
			IloNumArray scenMIPobj(fData_p->dataEnv);
			IloNumArray scenLGRobj(fData_p->dataEnv);
			IloNumArray scenLGRobj_update(fData_p->dataEnv);
			IloNumArray scenLPobj(fData_p->dataEnv);
			IloNumArray multiplier(fData_p->dataEnv, constr2Size);
			vector<float> dualAvg (constr2Size, 0.0);  // LP dual corresponding to constraints z(t) = x(t-1)
			vector<float> multiplierAvg (constr2Size, 0.0);  // LG dual corresponding to constraints z(t) = x(t-1)
			
			// IloNumArray scenFYobj(fData_p->dataEnv);
			// vector<float> dualAvg2 (constr2Size, 0.0);
			// IloNumArray val_y1(models[t].env);
			// IloNumArray val_y2(models[t].env);

			for ( k = 0; k < fData_p->numScen[t]; ++k )  // for each scenario at stage t
			{
				// cout << "Solving scenario " << k << " in stage " << t << endl;
				
				// update y2 coefficient with y2Scenarios[t][k]
				// models[t].obj.setLinearCoefs(models[t].y2, fData_p->y2Scenarios[t][k]);

				// update b_t with scenarios[t][k]
				for ( i = 0; i < uncertainArraySize; ++i )
					models[t].constr1[fData_p->uncertainIndex[i]].setLB(fData_p->scenarios[t][k][i]);
				// --------------------------------------------
				// cout << "RHS updated..." << endl;
				// --------------------------------------------

				// update the state variables z_t = x_{t-1}
				models[t].constr2.setBounds( candidateSol[t-1][j], candidateSol[t-1][j] );

				// --------------------------------------------
				// cout << "state variables updated..." << endl;
				// --------------------------------------------
				
				// solve models[t] as a MIP
				// cout << "Problem is MIP? " << models[t].cplex.isMIP() << endl;
				// models[t].cplex.exportModel("BW.lp");

				if ( lagrangianFlag + integerFlag )
				{
					// solve node IP, can be used in integer cut or serve as an ub for LGR
					if ( models[t].cplex.solve() ) // feasible
					{	
						// get solution status
						solStatus = models[t].cplex.getStatus();
						// cout << "MIP solved! Solution status: " << solStatus << endl;

						if ( solStatus == IloAlgorithm::Optimal )
						{
							scenMIPobj.add(models[t].cplex.getObjValue());

							/*
							fix local solutions y1, y2, solve the problem again to obtain dual multipliers
							
							models[t].cplex.getValues(val_y1, models[t].y1);
							models[t].cplex.getValues(val_y2, models[t].y2);
							
							if ( yFlag )
							{
								IloModel FY(models[t].env);
								FY.add(models[t].x);
								FY.add(models[t].y1);
								FY.add(models[t].y2);
								FY.add(models[t].z);
								FY.add(models[t].theta);
								FY.add(models[t].obj);
								FY.add(models[t].constr1);
								FY.add(models[t].constr2);
								FY.add(models[t].cuts);

								IloRangeArray cy1(models[t].env);
								IloRangeArray cy2(models[t].env);
								for ( i = 0; i < val_y1.getSize(); ++i )
									cy1.add(models[t].y1[i] == val_y1[i]);
								// for ( i = 0; i < val_y2.getSize(); ++i )
								// 	cy2.add(models[t].y2[i] == val_y2[i]);
								FY.add(cy1);
								FY.add(cy2);

								IloConversion relaxVarX(models[t].env, models[t].x, ILOFLOAT);
								IloConversion relaxVarY(models[t].env, models[t].y1, ILOFLOAT);
								FY.add(relaxVarX);
								FY.add(relaxVarY);

								IloCplex cplexFY(FY);
								cplexFY.setOut(models[t].env.getNullStream());

								if ( cplexFY.solve() )
								{
									solStatus = cplexFY.getStatus();
									if (solStatus == IloAlgorithm::Optimal )
									{
										scenFYobj.add(cplexFY.getObjValue());
										cplexFY.getDuals(multiplier, models[t].constr2);
										for (i = 0; i < multiplier.getSize(); ++i)
											dualAvg2[i] += multiplier[i] / fData_p->numScen[t];

										relaxVarX.end();
										relaxVarY.end();
										FY.end();
										cplexFY.end();
									}
									else // not optimal
									{
										cout << "Solution status: " << solStatus << endl;
										throw ("FY not optimal...");
									}
								}
								else // infeasible
								{
									char fileName[100];
									sprintf(fileName, "inf_FY.lp");
									cplexFY.exportModel(fileName);
									cout << "Solution status: " << cplexFY.getStatus() << endl;
									throw ("Problem with fixed y has no solution...");
								}
							}
							*/
						}
						else // not optimal
						{
							cout << "Solution status: " << solStatus << endl;
							throw ("MIP not optimal...");
						}
					}
					else // infeasible
					{
						cout << "Solution status: " << models[t].cplex.getStatus() << endl;
						throw ("MIP has no solution...");
					}
				}

				if ( bendersFlag + impvdBendersFlag + lagrangianFlag )
				{
					// solve models[t] LP relaxation
					IloModel LP(models[t].env);
					LP.add(models[t].x);
					LP.add(models[t].y1);
					LP.add(models[t].y2);
					LP.add(models[t].z);
					LP.add(models[t].theta);
					LP.add(models[t].obj);
					LP.add(models[t].constr1);
					LP.add(models[t].constr2);
					LP.add(models[t].cuts);
					IloConversion relaxVarX(models[t].env, models[t].x, ILOFLOAT);
					IloConversion relaxVarY(models[t].env, models[t].y1, ILOFLOAT);
					LP.add(relaxVarX);
					LP.add(relaxVarY);
					IloCplex cplexLP(LP);
					cplexLP.setOut(models[t].env.getNullStream());

					// solve models[t] as a LP
					// cout << "problem changed to LP? " << !(cplexLP.isMIP()) << endl;
					// cplexLP.exportModel("BWLP.lp");

					if ( cplexLP.solve() ) // feasible
					{
						// get solution status
						solStatus = cplexLP.getStatus();
						
						if ( solStatus == IloAlgorithm::Optimal )
						{
							// cout << "LP relaxation solved." << endl;
							// record objective function value
							scenLPobj.add(cplexLP.getObjValue());
							// record optimal dual multipliers of constraints z(t) = x(t-1)
							cplexLP.getDuals(multiplier, models[t].constr2);
							// cout << multiplier << endl;
							// for (int idx = 0; idx < 5; ++idx)
							// 	cout << (*it)[idx] << ", ";
							// cout << endl;
							for (i = 0; i < multiplier.getSize(); ++i)
								dualAvg[i] += multiplier[i] / fData_p->numScen[t];
							
							// change model back to MIP
							relaxVarX.end();
							relaxVarY.end();
							LP.end();
							cplexLP.end();

							// continue to solve Lagrangian if needed
							if ( impvdBendersFlag + lagrangianFlag )
							{
								// create a new model with models[t].env
								IloModel modelLGR(models[t].env);

								// add variables and constr1
								modelLGR.add(models[t].x);
								modelLGR.add(models[t].y1);
								modelLGR.add(models[t].y2);
								modelLGR.add(models[t].z);
								modelLGR.add(models[t].theta);
								modelLGR.add(models[t].constr1);
								modelLGR.add(models[t].cuts);

								// create new objective function with additional term -\pi_LP' z
								IloObjective objLGR = IloMinimize(models[t].env);
								IloExpr expr(models[t].env);
								expr = IloScalProd(fData_p->xCoef[t], models[t].x) + 
										IloScalProd(fData_p->y1Coef[t], models[t].y1) + 
										IloScalProd(fData_p->y2Coef[t], models[t].y2) +
										models[t].theta - 
										IloScalProd(multiplier, models[t].z);
								objLGR.setExpr(expr);
								modelLGR.add(objLGR);

								// create cplex algorithm for this model
								IloCplex cplexLGR(modelLGR);
								cplexLGR.setParam(IloCplex::EpGap, MIPTOL);
								// char fileName[100];
								// sprintf(fileName, "LGmodel_%d.lp", t);
								// cplexLGR.exportModel("LG.lp");
								cplexLGR.setOut(models[t].env.getNullStream());

								// solve IP to get improved Benders' cut
								if ( impvdBendersFlag )
								{
									// cout << "solving relaxed MIP... " << endl;
									if ( cplexLGR.solve() )
										scenLGRobj.add(cplexLGR.getObjValue());
								}
								
								// update lagrangian multipliers for a few iterations
								if ( lagrangianFlag )
								{
									if ( LEVEL_SWITCH )
									{
										IloNumArray2 pi(models[t].env);
										IloNumArray2 z(models[t].env);
										IloNumArray vpi(models[t].env);
										IloNumArray fval(models[t].env);
										IloNumArray2 grad(models[t].env);
										IloNumArray pi_new(models[t].env);

										if ( cplexLGR.solve() )
										{
											pi.add(multiplier);
											z.add(IloNumArray(models[t].env));
											cplexLGR.getValues(z[0], models[t].z);
											vpi.add(cplexLGR.getObjValue());
											// cout << vpi << endl;
											fval.add(vpi[0] + IloScalProd(pi[0], candidateSol[t-1][j]));
											grad.add(IloNumArray(models[t].env));
											for ( i = 0; i < candidateSol[t-1][j].getSize(); ++i )
												grad[0].add(candidateSol[t-1][j][i] - z[0][i]);

											double bestUB = models[t].cplex.getObjValue();//scenMIPobj[k];
											double bestLB, level, lambda = 0.35;
											double Delta = numeric_limits<double>::infinity();
											int best, iter = 0;

											// cout << "Level method iteration: " << iter << endl;
											while ( (Delta / bestUB > TOLLEVEL) && (IloScalProd(grad[iter], grad[iter]) > EPSILON) )
											{
												iter += 1;

												best = find_max(fval);
												bestLB = fval[best];
												
												Delta = bestUB - bestLB;
												level = bestUB - lambda * Delta;

												// double approxUB = pwl(fval, grad, pi, bestUB);
												
												pi_new = projection(pi[best], fval, grad, pi, level);
												cout << " current Delta: " << Delta << endl;
												// cout << " optimal value of the scenario MIP: " << bestUB << endl;
												// cout << fval << endl;
												// cout << " approx. upper bound: " << approxUB << endl;

												for ( i = 0; i < pi_new.getSize(); ++i )
													objLGR.setLinearCoef(models[t].z[i], -pi_new[i]);
												
												if ( cplexLGR.solve() )
												{
													pi.add(pi_new);
													z.add(IloNumArray(models[t].env));
													cplexLGR.getValues(z[iter], models[t].z);
													vpi.add(cplexLGR.getObjValue());
													// cout << pi_new.getSize() << " " << candidateSol[t-1][j].getSize() << endl;
													fval.add(vpi[iter] + IloScalProd(pi_new, candidateSol[t-1][j]));
													grad.add(IloNumArray(models[t].env));
													for ( i = 0; i < candidateSol[t-1][j].getSize(); ++i )
														grad[iter].add(candidateSol[t-1][j][i] - z[iter][i]);
												}
												else
													cout << "No solution for current pi. " << endl;


												// cout << "Level method iteration: " << iter <<
												        // "  Current Delta: " << Delta << endl;
											}
										}
										for (i = 0; i < pi_new.getSize(); ++i)
											multiplierAvg[i] += pi_new[i] / fData_p->numScen[t];
										scenLGRobj_update.add(fval[fval.getSize()-1]);

										pi.end();
										z.end();
										vpi.end();
										fval.end();
										grad.end();
										pi_new.end();
									}
									else
									{
										LGupdate(modelLGR, cplexLGR, objLGR,
												multiplier, scenLGRobj_update,
												models[t], candidateSol[t-1][j], scenMIPobj[k]);
										for (i = 0; i < multiplier.getSize(); ++i)
											multiplierAvg[i] += multiplier[i] / fData_p->numScen[t];
									}
									// cout << "LG update finished" << endl;
									// cout << "IP solution " << scenMIPobj << endl;
									// cout << "LGR solution " << scenLGRobj_update << endl;
									// cout << "========================" << endl;
								}

								// release memory
								// cout << "relaxed MIP solved..." << endl;
								expr.end();
								objLGR.end();
								modelLGR.end();
								cplexLGR.end();
							}
						}
						else // not optimal
						{
							cout << "Solution status: " << solStatus << endl;
							throw ("LP not optimal...");
						}
					}
					else // infeasible
					{
						char fileName[100];
						sprintf(fileName, "inf_LP.lp");
						cplexLP.exportModel(fileName);
						cout << "Solution status: " << models[t].cplex.getStatus() << endl;
						throw ("LP has no solution...");
					}
			
				}
			} // End of loop over all scenarios in stage t

			// cout << "MIP obj:" << scenMIPobj << endl;
			// cout << "fixY obj:" << scenFYobj << endl;
			// cout << "LP obj:" << scenLPobj << endl;

			// construct and add integer L-shaped cut
			if ( integerFlag )
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				IloNum mipObjAvg = IloSum(scenMIPobj) / fData_p->numScen[t];
				// cout << "average mip objective: " << mipObjAvg << endl;
				for ( i = 0; i < constr2Size; ++i )
				{
					if ( abs( candidateSol[t-1][j][i] ) < EPSILON )
						expr += (mipObjAvg - fData_p->valueLB[t-1]) * models[t-1].x[i];
					else
						expr -= (mipObjAvg - fData_p->valueLB[t-1]) * models[t-1].x[i];
				}
				rhs = mipObjAvg - (mipObjAvg - fData_p->valueLB[t-1]) * IloSum(candidateSol[t-1][j]);
				// cout << fData_p->valueLB[t-1] - (mipObjAvg - fData_p->valueLB[t-1]) * (IloSum(candidateSol[t-1][j]) - 1) << endl;
				// cout << rhs << endl;
				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				// cout << "L-shpaed cut added." << endl;
				expr.end();
			}

			// construct and add (Strengthened) Benders' cuts
			if ( bendersFlag + impvdBendersFlag )
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				for ( i = 0; i < models[t-1].x.getSize(); ++i )
					expr -= dualAvg[i] * models[t-1].x[i];
				if ( impvdBendersFlag )
					rhs = IloSum(scenLGRobj) / fData_p->numScen[t];
				else
				{
					double v = 0.0;
					for ( i = 0; i < constr2Size; ++i )
						v += dualAvg[i] * candidateSol[t-1][j][i];
					rhs = IloSum(scenLPobj) / fData_p->numScen[t] - v;
				}
				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				// -----------------------------------------
				// cout << expr << " >= " << rhs << endl;
				// cout << IloSum(scenLPobj) / fData_p->numScen[t] << endl;
				// cout << IloSum(scenLGRobj) / fData_p->numScen[t] << endl;
				// for (i = 0; i < constr2Size; ++i)
					// cout << dualAvg[i] << "   " << (*it)[i] << endl;
				// cout << inner_product(dualAvg.begin(), dualAvg.end(), (*it).begin(), 0.0) << endl;
				// -----------------------------------------
				// if ( impvdBendersFlag )
				// 	cout << "Strengthened Benders' cut added." << endl;
				// else
				// 	cout << "Benders' cut added." << endl;
				expr.end();
			}

			/*
			if ( yFlag )
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				for ( i = 0; i < models[t-1].x.getSize(); ++i )
					expr -= dualAvg2[i] * models[t-1].x[i];
				
				double v = 0.0;
				for ( i = 0; i < constr2Size; ++i )
					v += dualAvg2[i] * candidateSol[t-1][j][i];
				rhs = IloSum(scenLPobj) / fData_p->numScen[t] - v;
				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				// -----------------------------------------
				// cout << expr << " >= " << rhs << endl;
				// cout << IloSum(scenLPobj) / fData_p->numScen[t] << endl;
				// cout << IloSum(scenLGRobj) / fData_p->numScen[t] << endl;
				// for (i = 0; i < constr2Size; ++i)
					// cout << dualAvg[i] << "   " << (*it)[i] << endl;
				// cout << inner_product(dualAvg.begin(), dualAvg.end(), (*it).begin(), 0.0) << endl;
				// -----------------------------------------
				// if ( impvdBendersFlag )
				// 	cout << "Strengthened Benders' cut added." << endl;
				// else
				// 	cout << "Benders' cut added." << endl;
				expr.end();
			}
			*/

			if ( lagrangianFlag )
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				for ( i = 0; i < models[t-1].x.getSize(); ++i )
					expr -= multiplierAvg[i] * models[t-1].x[i];
				double v = 0.0;
				for ( i = 0; i < constr2Size; ++i )
					v += multiplierAvg[i] * candidateSol[t-1][j][i];
				rhs = IloSum(scenLGRobj_update) / fData_p->numScen[t] - v ;
				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				// cout << "Lagrangian cut added." << endl;
				expr.end();
			}

			// free momery
			scenMIPobj.end();
			if ( bendersFlag + impvdBendersFlag + lagrangianFlag )
			{
				scenLPobj.end();
				scenLGRobj.end();
				scenLGRobj_update.end();
				multiplier.end();
				dualAvg.clear();
				multiplierAvg.clear();
			}

		} // End of loop over unique candidate solutions
	} // End of loop over stages

	// cout << "================================" << endl;
	cout << "Solving the master problem." << endl;

	// solve models[0] as a MIP
	if ( models[0].cplex.solve() ) // feasible
	{
		// get solution status
		solStatus = models[0].cplex.getStatus();
		// cout << "solution status obtained." << endl;

		if ( solStatus == IloAlgorithm::Optimal )
		{
			// update lower bound lb as the optimal objective function value
			lb.add(models[0].cplex.getObjValue());
			// IloNumArray vals(fData_p->dataEnv, models[0].x.getSize());
			// models[0].cplex.getValues(vals, models[0].x);
			// for ( i = 0; i < vals.getSize(); ++i )
			// 	vals[i] = round(vals[i]);
			// vals.end();
		}
		else // not optimal
		{
			cout << "First-stage problem status: " << solStatus << endl;
			throw ("First-stage problem not optimal...");
		}
	}
	else // infeasible
	{
		cout << "First stage problem status: " << models[0].cplex.getStatus() << endl;
		throw ("First-stage problem infeasible...");
	}

} // End of backward pass to refine value functions


int find_max ( IloNumArray x )
{
	int i, imax, size = x.getSize();
	if ( size == 0 )
		cout << " No element exists in the array, please check!! " << endl;
	else
	{
		double vmax = -numeric_limits<double>::infinity();
		for ( i = 0; i < size; ++i )
		{
			if ( x[i] > vmax )
			{
				imax = i;
				vmax = x[i];
			}
		}
	}
	return imax;
}

IloNumArray projection ( IloNumArray xstar, IloNumArray fval, IloNumArray2 grad, IloNumArray2 x0, double level )
{
	IloEnv env;
	IloInt dim = xstar.getSize();

	IloModel m(env);

	IloNumVarArray x(env, dim, -IloInfinity, IloInfinity, ILOFLOAT);
	m.add(x);

	IloObjective obj(env);
	IloExpr objExpr(env);
	objExpr = IloScalProd(x, x) - 2 * IloScalProd(x, xstar);
	obj.setExpr(objExpr);
	obj.setSense(IloObjective::Minimize);
	m.add(obj);

	IloRangeArray constr(env);
	int i, nrow = x0.getSize();
	for (i = 0; i < nrow; ++i)
	{
		IloExpr expr(env);
		expr = fval[i] + IloScalProd(grad[i], x) - IloScalProd(grad[i], x0[i]);
		// cout << fval[i] << " " << IloScalProd(grad[i], x0[i]) << " " << level << endl;
		constr.add(expr >= level);
	}
	m.add(constr);

	IloCplex cplex(m);
	cplex.setOut(env.getNullStream());

	IloNumArray opt(xstar.getEnv());
	if ( cplex.solve() )
		cplex.getValues(opt, x);
	else
	{
		cout << " No solution exists for the projection problem. Please check!! " << endl;
		cout << " current level: " << level << endl;
		cplex.exportModel("projection.lp");
	}
	env.end();

	return opt;
}

double pwl(IloNumArray fval, IloNumArray2 grad, IloNumArray2 x0, double opt )
{
	IloEnv env;
	IloInt dim = x0[0].getSize();

	IloModel m(env);

	IloNumVarArray x(env, dim, -IloInfinity, IloInfinity, ILOFLOAT);
	IloNumVar t(env, -IloInfinity, 1000000, ILOFLOAT);
	m.add(x);
	m.add(t);

	IloObjective obj(env);
	IloExpr objExpr(env);
	objExpr = t;
	obj.setExpr(objExpr);
	obj.setSense(IloObjective::Maximize);
	m.add(obj);

	IloRangeArray constr(env);
	int i, nrow = x0.getSize();
	for (i = 0; i < nrow; ++i)
	{
		IloExpr expr(env);
		expr = t - fval[i] - IloScalProd(grad[i], x) + IloScalProd(grad[i], x0[i]);
		constr.add(expr <= 0);
	}
	m.add(constr);

	IloCplex cplex(m);
	cplex.setOut(env.getNullStream());

	if ( ! cplex.solve() )
	{
		cout << " No solution exists for the relaxation problem. Please check!! " << endl;
		cplex.exportModel("approx.lp");
	}

	double approxUB = cplex.getObjValue();
	env.end();

	return approxUB;
}

void LGupdate ( IloModel & modelLGR, IloCplex & cplexLGR, IloObjective & objLGR,  
		IloNumArray & multiplier, IloNumArray & scenObj,
		const Model model, const IloNumArray state, const double ub )
{
	IloNumArray gradient(multiplier.getEnv(), multiplier.getSize());
	int i, j;
	int numIter = 1000;
	double norm_grad, stepSize, objVal, gap;
	bool optimalFlag = 0;

	for ( i = 0; i < numIter; ++i)
	{
		// cout << "maximum scale of current lambda: " << IloMax( IloMax(multiplier), IloAbs(IloMin(multiplier)) ) << endl;
		cplexLGR.solve();
		objVal = cplexLGR.getObjValue();
		double temp = IloScalProd(multiplier, state);
		gap = (ub - objVal - temp) / ub;
		cout << "iteration: " << i <<
			"   Current gap: " << gap <<
			"   Absolute difference: " << ub - objVal - temp << endl;

		// optimality check
		if ( gap < 1e-4 )
		{
			optimalFlag = 1;
			// cout << "Lagrangian problem is solved to optimality." << endl;
			// cout << "optimal lambda: " << multiplier << endl;
			// cout << "maximum scale of lambda: " << IloMax( IloMax(multiplier), IloAbs(IloMin(multiplier)) ) << endl;
			break;
		}
		else if ( i < numIter - 1 )
		{
			cout << "\e[A";
			norm_grad = 0.0;
			stepSize = 20 / sqrt(i+1);
		
			// obtain and revise gradient
			cplexLGR.getValues(gradient, model.z);
			for ( j = 0; j < gradient.getSize(); ++j)
			{
				gradient[j] -= state[j];
				norm_grad += gradient[j] * gradient[j];
			}
			// update multiplier and objective function
			for ( j = 0; j < multiplier.getSize(); ++j )
			{
				multiplier[j] -= stepSize / (norm_grad + 1e-3) * gradient[j];
				objLGR.setLinearCoef(model.z[j], -multiplier[j]);
			}
		}
	}
	scenObj.add(objVal + IloScalProd(multiplier, state));
}

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
	cerr << "Usage:  " << progname << " arg1 arg2 arg3 arg4 [arg5]" << endl;
	cerr << "At least 4 parameters must be specified." << endl;
	cerr << "arg1: 0 -- turn off Benders' cuts;" << endl;
	cerr << "      1 -- turn on Benders' cuts." << endl;
	cerr << "arg2: 0 -- turn off Strengthened Benders' cuts;" << endl;
	cerr << "      1 -- turn on Strengthened Benders' cuts." << endl;
	cerr << "arg3: 0 -- turn off Lagragian cuts;" << endl;
	cerr << "      1 -- turn on Lagragian cuts." << endl;
	cerr << "arg4: 0 -- turn off Integer L-shaped cuts;" << endl;
	cerr << "      1 -- turn on Integer L-shaped cuts." << endl;
	cerr << "arg5: [optional] used as the seed of the random number generator." << endl;
	cerr << "      If not provided, system will generate one automatically." << endl;
} // END usage
