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

void readData (formatData * fData_p)
{
	// create data environment
	fData_p->dataEnv = IloEnv();
	IloEnv * dataEnv = & fData_p->dataEnv;

	cout << "Start reading data from files... " << endl;

	// read number of stages
	readArray<IloInt> (fData_p->numStage, "data/numStage.dat");
	IloInt numStage = fData_p->numStage;

	// read initial state of binary variables
	fData_p->initState = IloNumArray(*dataEnv, numStage);
	readArray<IloNumArray> (fData_p->initState, "data/initState.dat");

	// read prescribed lower bound for value function at each stage
	fData_p->valueLB = IloNumArray(*dataEnv, numStage);
	readArray<IloNumArray> (fData_p->valueLB, "data/valueLB.dat");

	// read number of samples drawn in the forward pass (number of candidate solns)
	readArray<IloInt> (fData_p->numFWsample, "data/numFWsample.dat");

	// read number of scenarios available at each stage
	fData_p->numScen = IloIntArray(*dataEnv, numStage);
	readArray<IloIntArray> (fData_p->numScen, "data/numScen.dat");

	// compute the total number of scenarios
	fData_p->totalScen = fData_p->numScen[0];
	for (int t = 1; t < numStage; ++t)
		fData_p->totalScen *= fData_p->numScen[t];

	// read objective coefficients for x (binary) variables 
	fData_p->xCoef = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->xCoef, "data/xCoef.dat");

	// read objective coefficients for y1 (integral) variables
	fData_p->y1Coef = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->y1Coef, "data/y1Coef.dat");

	// read objective coefficients for y2 (continuous) variables
	fData_p->y2Coef = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->y2Coef, "data/y2Coef.dat");

	// read matrix A in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t <= b_t
	fData_p->Amatrix = IloNumArray3(*dataEnv, numStage);
	readArray<IloNumArray3> (fData_p->Amatrix, "data/Amatrix.dat");
	
	// read matrix B in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t <= b_t
	fData_p->Bmatrix = IloNumArray3(*dataEnv, numStage);
	readArray<IloNumArray3> (fData_p->Bmatrix, "data/Bmatrix.dat");
	
	// read matrix W1 in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t <= b_t
	fData_p->W1matrix = IloNumArray3(*dataEnv, numStage);
	readArray<IloNumArray3> (fData_p->W1matrix, "data/W1matrix.dat");
	
	// read matrix W2 in constraint A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t <= b_t
	fData_p->W2matrix = IloNumArray3(*dataEnv, numStage);
	readArray<IloNumArray3> (fData_p->W2matrix, "data/W2matrix.dat");
	
	// read rhs b in constraint A_tx_t + W_ty_t + B_tz_t <= b_t
	fData_p->bRhs = IloNumArray2(*dataEnv, numStage);
	readArray<IloNumArray2> (fData_p->bRhs, "data/bRhs.dat");

	// read uncertain indices in b (every stage have the same uncertain parameter)
	fData_p->uncertainIndex = IloIntArray(*dataEnv, numStage);
	readArray<IloIntArray> (fData_p->uncertainIndex, "data/uncertainIndex.dat");

	// read generated rhs scenarios at each stage
	// dim1: numStage; dim2: numScen[t]; dim3: uncertainIndex.getSize()
	fData_p->scenarios = IloNumArray3(*dataEnv, numStage);
	readArray<IloNumArray3> (fData_p->scenarios, "data/scenarios.dat");

	// read generated y2 coefficient scenarios at each stage
	// dim1: numStage; dim2: numScen[t]; dim3: nType * SUBPERIOD
	fData_p->y2Scenarios = IloNumArray3(*dataEnv, numStage);
	readArray<IloNumArray3> (fData_p->y2Scenarios, "data/y2Scenarios.dat");

	return;
} // End of readData

void buildModel (Model * models, formatData * fData_p)
{
	cout << "Start to build one model for each stage..." << endl;

	int t, i, j, k;
	IloInt numStage = fData_p->numStage;

	// extract decision variable dimensions
	IloInt x_dim = fData_p->xCoef[0].getSize();   // state variables
	IloInt y1_dim = fData_p->y1Coef[0].getSize(); // current stage integral variables
	IloInt y2_dim = fData_p->y2Coef[0].getSize(); // current stage continuous variables
	IloInt z_dim = x_dim;						  // copy of state variables
	IloInt numRows = fData_p->bRhs[0].getSize();

	for (t = 0; t < numStage; ++t)
	{
		// create cplex environment and initialize model
		models[t].env = IloEnv();
		IloEnv currentEnv = models[t].env;
		models[t].mod = IloModel(currentEnv);

		// create variables and add them to model
		char varName[100];
		
		models[t].x = IloNumVarArray(currentEnv, x_dim, 0.0, 1.0, ILOINT);
		for ( i = 0; i < x_dim; ++i )
		{
			sprintf(varName, "x_%d", i+1);
			models[t].x[i].setName(varName);
		}
		models[t].mod.add(models[t].x);
		
		models[t].y1 = IloNumVarArray(currentEnv, y1_dim, 0.0, IloInfinity, ILOINT);
		models[t].y2 = IloNumVarArray(currentEnv, y2_dim, 0.0, IloInfinity, ILOFLOAT);
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
		models[t].mod.add(models[t].y1);
		models[t].mod.add(models[t].y2);
		
		
		models[t].z = IloNumVarArray(currentEnv, z_dim, 0.0, 1.0, ILOFLOAT);
		for ( i = 0; i < z_dim; ++i )
		{
			sprintf(varName, "z_%d", i+1);
			models[t].z[i].setName(varName);
		}
		models[t].mod.add(models[t].z);

		models[t].theta = IloNumVar(currentEnv, 0.0, IloInfinity);
		sprintf(varName, "theta");
		models[t].theta.setName(varName);
		models[t].mod.add(models[t].theta);

		// create objective function and add to model
		models[t].obj = IloObjective(currentEnv);
		IloExpr objExpr(currentEnv);

		objExpr = IloScalProd(models[t].x, fData_p->xCoef[t]);
		objExpr += IloScalProd(models[t].y1, fData_p->y1Coef[t]);
		objExpr += IloScalProd(models[t].y2, fData_p->y2Coef[t]);
		objExpr += models[t].theta;
		models[t].obj.setExpr(objExpr);
		objExpr.end();
		models[t].obj.setSense(IloObjective::Minimize);
		char objName[100];
		sprintf(objName, "objective_%d", t);
		models[t].obj.setName(objName);
		models[t].mod.add(models[t].obj);

		// create constraints
		// Add constraints A_tx_t + B_tz_t + W1_ty1_t + W2_ty2_t >= b_t
		models[t].constr1 = IloRangeArray(currentEnv);
	
		for (i = 0; i < numRows; ++i)
		{
			IloExpr expr(currentEnv);
			expr = IloScalProd(models[t].x, fData_p->Amatrix[t][i]);
			expr += IloScalProd(models[t].y1, fData_p->W1matrix[t][i]);
			expr += IloScalProd(models[t].y2, fData_p->W2matrix[t][i]);
			expr += IloScalProd(models[t].z, fData_p->Bmatrix[t][i]);
			models[t].constr1.add(expr >= fData_p->bRhs[t][i]);
			expr.end();
		} // end of rows

		models[t].mod.add(models[t].constr1);

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

		// Initialize cut constraints
		models[t].cuts = IloRangeArray(currentEnv);

		// create cplex algorithms
		models[t].cplex = IloCplex(models[t].mod);

		// set model algorithm
		models[t].cplex.setParam(IloCplex::RootAlg, IloCplex::Primal);

		// set model solve output
		models[t].cplex.setOut(models[t].env.getNullStream());

		// set model warning message
		models[t].cplex.setWarning(models[t].env.getNullStream());

		// write model to lp file
		char fileName[100];
		sprintf(fileName, "model_%d.lp", t);
		models[t].cplex.exportModel(fileName);

	} // End of t-loop for model construction
	return;
} // End of function buildModel

void getSamplePaths (IloNumArray3 & samplePaths, IloNumArray3 & coefSamplePaths, formatData * fData_p)
{
	// cout << "Sampling forward paths from scenarios..." << endl;
	IloEnv * currentEnv = &(fData_p->dataEnv);
	IloInt numPaths = fData_p->numFWsample;
	IloInt numStage = fData_p->numStage;
	// IloNumArray2 emptyPath(*currentEnv);

	for (int p = 0; p < numPaths; ++p)
	{
		// cout << "sample path p = " << p << endl;
		
		// initialize each sample path
		samplePaths.add(IloNumArray2(*currentEnv));
		samplePaths[p].add(fData_p->scenarios[0][0]);
		// cout << "good here" << endl;
		// cout << fData_p->y2Scenarios[0][0] << endl;
		coefSamplePaths.add(IloNumArray2(*currentEnv));
		coefSamplePaths[p].add(fData_p->y2Scenarios[0][0]);
		// cout << "good here" << endl;
		

		// draw sample from each stage
		for (int t = 1; t < numStage; ++t)
		{
			double pdf = 1.0/fData_p->numScen[t];
			// generate a random number between 0 and 1
			double U = genrand64_real1();
			// cout << "pdf: " << pdf << endl;
			// cout << "random number generated: " << U << endl;
			// cout << int(U/pdf) << ",";
			// cout << "scenario chosen: " << int(U/pdf) << endl;
			int chosen = int(U/pdf);
			samplePaths[p].add(fData_p->scenarios[t][chosen]);
			coefSamplePaths[p].add(fData_p->y2Scenarios[t][chosen]);

		} // End of stage for-loop
		// cout << " " << endl;
	} // End of sample paths for loop
	return;
} // End of getSamplePathss

void forward (Model * models, formatData * fData_p, const IloNumArray3 samplePaths, const IloNumArray3 coefSamplePaths,
			  IloNumArray3 & candidateSol, IloNumArray & ub_c, IloNumArray & ub_l, IloNumArray & ub_r)
{
	cout << "Start the forward process..." << endl;

	int p, t, i, index;

	IloInt constr1Size = models[0].constr1.getSize();
	IloInt constr2Size = models[0].constr2.getSize();
	IloInt sampleSize = fData_p->numFWsample;
	IloInt uncertainArraySize = fData_p->uncertainIndex.getSize();
	IloAlgorithm::Status solStatus;

	// create an array to record the objective function value for each sample path
	IloNumArray sampleObj(fData_p->dataEnv, sampleSize);

	for ( p = 0; p < sampleSize; ++p)
		sampleObj[p] = 0.0;

	cout << "number of sample paths: " << sampleSize << endl;
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
				// update objective coefficients
				models[t].obj.setLinearCoefs(models[t].y2, coefSamplePaths[p][t]);

				// update b_t with samplePaths[p][t]
				for ( i = 0; i < uncertainArraySize; ++i )
				{
					index = fData_p->uncertainIndex[i];
					models[t].constr1[index].setLB(samplePaths[p][t][i]);
				}

				// update the state variables z_t = x_{t-1}
				for ( i = 0; i < constr2Size; ++i )
					models[t].constr2[i].setBounds(candidateSol[p][t-1][i], candidateSol[p][t-1][i]);
			} // End of update problem models[t]

			// solve the current MIP model
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
					for ( i = 0; i < vals.getSize(); ++i )
						vals[i] = round(vals[i]);
					candidateSol[p].add(vals);
					
					IloNum costToGo = models[t].cplex.getValue(models[t].theta);
					// cout << "!!!!!!!!!! OBJECTIVE: " << models[t].cplex.getObjValue() << endl;

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
	// cout << "objective value for forward sample paths: " << sampleObj << endl;
	/*
	cout << ub_l << endl;
	cout << ub_c << endl;
	cout << ub_r << endl;
	*/
	cout << "95\% CI for the upper bound: [ " << center - halfLength <<  ", " << center + halfLength << " ]." << endl;

	// free momory
	sampleObj.end();
} // End of forward pass to generate candidate solutions

void backward (Model * models, formatData * fData_p, const IloNumArray3 candidateSol,
		IloNumArray & lb, unordered_set<string> & masterSol,
		const bool bendersFlag, const bool impvdBendersFlag, const bool integerFlag)
{
	cout << "================================" << endl;
	cout << "Start the backward process..." << endl;

	IloInt p, t, k, i, index;
	IloNum rhs;
	IloInt sampleSize = fData_p->numFWsample;
	IloInt uncertainArraySize = fData_p->uncertainIndex.getSize();
	IloInt constr1Size = models[0].constr1.getSize();
	IloInt constr2Size = models[0].constr2.getSize();
	IloAlgorithm::Status solStatus;

	unordered_set<string> uniqueSet;
	IntMatrix uniqueSol;

	for ( t = fData_p->numStage - 1; t > 0; --t ) // for each stage
	{
		cout << "======================" << endl;
		cout << "Current stage t =  " << t << endl;

		// find out unique candidate solutions at current stage among candidateSol[p][t] for p = 1,..., numFWsample

		// create a 2-d array to store candidateSol[p][t] for p = 1,..., numFWsample
		for ( p = 0; p < fData_p->numFWsample; ++p )
			uniqueSet.insert(toString(candidateSol[p][t-1]));

		// convert uniqueSet(unordered_set) to uniqueSol(2-d vector)
		for ( auto it = uniqueSet.begin(); it != uniqueSet.end(); ++it )
			uniqueSol.push_back(getDigit(*it));

		for ( auto it = uniqueSol.begin(); it != uniqueSol.end(); ++it )  // 
		{
			// cout << "Total number of scenarios at this stage: " << fData_p->numScen[t] << endl;
			// create arrays to store MIP, Lagrangian with optimal LP dual, and LP optimal values
			IloNumArray scenMIPobj(fData_p->dataEnv);
			IloNumArray scenLGRobj(fData_p->dataEnv);
			IloNumArray scenLPobj(fData_p->dataEnv);
			IloNumArray vals(fData_p->dataEnv, constr2Size);
			// IloNumArray dualAvg(fData_p->dataEnv, constr2Size);
			vector<float> dualAvg ( (*it).size(), 0.0);

			for ( k = 0; k < fData_p->numScen[t]; ++k )  // for each scenario at stage t
			{
				// cout << "scenario " << k << " in stage " << t << endl;
				// update y2 coefficient with y2Scenarios[t][k]
				models[t].obj.setLinearCoefs(models[t].y2, fData_p->y2Scenarios[t][k]);

				// update b_t with scenarios[t][k]
				for ( i = 0; i < uncertainArraySize; ++i )
				{
					index = fData_p->uncertainIndex[i];
					models[t].constr1[index].setLB(fData_p->scenarios[t][k][i]);
				}

				// update the state variables z_t = x_{t-1}
				for ( i = 0; i < constr2Size; ++i )
					models[t].constr2[i].setBounds( (*it)[i], (*it)[i] );

				// solve models[t] as a MIP
				// cout << "Problem is MIP? " << models[t].cplex.isMIP() << endl;
				// models[t].cplex.exportModel("mip.lp");

				if ( integerFlag )
				{
					if ( models[t].cplex.solve() ) // feasible
					{	
						// get solution status
						// cout << "solution status: optimal" << endl;
						solStatus = models[t].cplex.getStatus();

						if ( solStatus == IloAlgorithm::Optimal )
						{
							// record objective function value
							scenMIPobj.add(models[t].cplex.getObjValue());
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

				if ( bendersFlag + impvdBendersFlag )
				{
					// solve models[t] LP relaxation
					IloConversion relaxVarX(models[t].env, models[t].x, ILOFLOAT);
					IloConversion relaxVarY(models[t].env, models[t].y1, ILOFLOAT);
					models[t].mod.add(relaxVarX);
					models[t].mod.add(relaxVarY);

					// solve models[t] as a LP
					//cout << "problem changed to LP? " << !(models[t].cplex.isMIP()) << endl;
					// models[t].cplex.exportModel("lpr.lp");

					if ( models[t].cplex.solve() ) // feasible
					{
						// get solution status
						solStatus = models[t].cplex.getStatus();

						if ( solStatus == IloAlgorithm::Optimal )
						{
							// record objective function value
							scenLPobj.add(models[t].cplex.getObjValue());
							// record optimal dual multipliers
							models[t].cplex.getDuals(vals, models[t].constr2);
							for (i = 0; i < vals.getSize(); ++i)
								dualAvg[i] += vals[i] / fData_p->numScen[t];
							
							// change model back to MIP
							relaxVarX.end();
							relaxVarY.end();

							// continue to solve Lagrangian if impvdBendersFlag is 1
							if ( impvdBendersFlag )
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

								// create new objective function with additional term \pi_LP' z
								IloObjective newObj = IloMinimize(models[t].env);
								IloExpr expr(models[t].env);
								expr += IloScalProd(fData_p->xCoef[t], models[t].x);
								expr += IloScalProd(fData_p->y1Coef[t], models[t].y1);
								expr += IloScalProd(fData_p->y2Scenarios[t][k], models[t].y2);
								expr += models[t].theta;
								expr -= IloScalProd(vals, models[t].z);
								newObj.setExpr(expr);
								modelLGR.add(newObj);

								// create cplex algorithm for this model
								IloCplex cplexLGR(modelLGR);
								char fileName[100];
								sprintf(fileName, "LGmodel_%d.lp", t);
								// cplexLGR.exportModel(fileName);
								cplexLGR.setOut(models[t].env.getNullStream());

								// solve new model
								if ( cplexLGR.solve() )
								{
									double consObj = 0;
									for ( i = 0; i < vals.getSize(); ++i )
										consObj += vals[i] * (*it)[i];
								//	cout << consObj << endl;
									scenLGRobj.add(cplexLGR.getObjValue() + consObj);
								}
								// cout << "vals: " << vals << endl;
								// cout << "LP relaxation optimal value: " << scenLPobj << endl;
								// cout << "LG relaxation optimal value: " << scenLGRobj << endl;

								//system("read");

								//cout << models[t].mod << endl;;
								//system("read");
								//cout << modelLGR << endl;
								//system("read");

								// release memory
								expr.end();
								newObj.end();
								modelLGR.end();
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
						cout << "Solution status: " << models[t].cplex.getStatus() << endl;
						throw ("LP has no solution...");
					}
			
				}
			} // End of loop over all scenarios in stage t

			// cout << "MIP objective: " << scenMIPobj << "  MIP Average: " << IloSum(scenMIPobj) / fData_p->numScen[t] << endl;
			// cout << "LP objective: "  << scenLPobj  << "  LP Average: "  << IloSum(scenLPobj) / fData_p->numScen[t] << endl;

			// construct and add integer L-shaped cut
			if ( integerFlag )
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				IloNum mipObjAve = IloSum(scenMIPobj) / fData_p->numScen[t];
				for ( unsigned j = 0; j < (*it).size(); ++j )
				{
					if ( abs( (*it)[j]) < EPSILON )
						expr += (mipObjAve - fData_p->valueLB[t-1]) * models[t-1].x[j];
					else
						expr -= (mipObjAve - fData_p->valueLB[t-1]) * models[t-1].x[j];
				}
			
				rhs = fData_p->valueLB[t-1] - (mipObjAve - fData_p->valueLB[t-1]) * (accumulate( (*it).begin(), (*it).end(), 0) -1);
				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				cout << "L-shpaed cut added." << endl;
				expr.end();
			}

			// construct and add Benders cut
			if ( bendersFlag + impvdBendersFlag )
			{
				IloExpr expr(models[t-1].env);
				expr = models[t-1].theta;
				for ( i = 0; i < models[t-1].x.getSize(); ++i )
				{
					expr -= dualAvg[i] * models[t-1].x[i];
				}
				if ( impvdBendersFlag )
					rhs = IloSum(scenLGRobj) / fData_p->numScen[t] - inner_product(dualAvg.begin(), dualAvg.end(), (*it).begin(), 0);
				else
					rhs = IloSum(scenLPobj) / fData_p->numScen[t] - inner_product(dualAvg.begin(), dualAvg.end(), (*it).begin(), 0);
				models[t-1].cuts.add(expr >= rhs);
				models[t-1].mod.add(expr >= rhs);
				if ( impvdBendersFlag )
					cout << "Improved Benders' cut added." << endl;
				else
					cout << "Benders' cut added." << endl;
				expr.end();
			}

			// free momery
			scenMIPobj.end();
			if ( bendersFlag + impvdBendersFlag )
			{
				scenLPobj.end();
				scenLGRobj.end();
				vals.end();
				dualAvg.clear();
			}

		} // End of loop over unique solutions

		uniqueSet.clear();
		uniqueSol.clear();

	} // End of loop over stages

	cout << "================================" << endl;
	cout << "Solving the master problem: " << endl;
	//cout << models[0].mod << endl;

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
			IloNumArray vals(fData_p->dataEnv, models[0].x.getSize());
			models[0].cplex.getValues(vals, models[0].x);
			// cout << "solution value obtained." << endl;
			for ( i = 0; i < vals.getSize(); ++i )
				vals[i] = round(vals[i]);
			masterSol.insert(toString(vals));
			vals.end();
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
		throw ("Relaxed first-stage problem infeasible...");
	}

} // End of backward pass to refine value functions


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
	cerr << "Usage:  " << progname << " arg1 [arg2]" << endl;
	cerr << "At least one parameter must be specified." << endl;
	cerr << "arg1: 0 -- turn off Benders' cuts;" << endl;
	cerr << "      1 -- turn on Benders' cuts." << endl;
	cerr << "arg2: [optional] used as the seed of the random number generator." << endl;
	cerr << "      If not provided, system will generate one automatically." << endl;
} // END usage
