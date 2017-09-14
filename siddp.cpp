#include <ilcplex/ilocplex.h>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <chrono>
#include <unordered_set>
#include <algorithm> // std::max

#include "global.h"
#include "functions.h"
#include "mt64.h"


ILOSTLBEGIN

int main (int argc, char *argv[])
{
	try 
	{
		if ( argc != 6 && argc != 7 )
		{
			usage (argv[0]);
			throw (-1);
		}

		bool bendersFlag, impvdBendersFlag, lagrangianFlag, integerFlag;
		bool sampling = atoi(argv[5]);

		unsigned long long seed;
		if ( argc == 7 )
			seed = atoi(argv[6]);
		else
			seed = chrono::system_clock::now().time_since_epoch().count();
		init_genrand64(seed);

		// create a new formated data structure and read data from file
		formatData fData;
		formatData * fData_p = & fData;
		readData (fData_p);

		cout << "All data has been read into fData." << endl;
		cout << "==================================================" << endl;
		cout << "total number of stages: " << fData.numStage << endl;
		cout << "number of scenarios at each stage: " << fData.numScen << endl;
		cout << "total number of scenarios: " << fData.totalScen << endl;
		cout << "sample size in the forward pass: " << fData.numFWsample << endl;
		cout << "==================================================" << endl;

		// construct models from fData
		Model * models = new Model[fData.numStage];
		buildModel(models, fData_p);
		cout << "Model construction completed..." << endl;

		// start siddp method
		// clock_t startTime = clock();
		chrono::time_point<chrono::system_clock> start, end;
		chrono::duration<double> elapsed_seconds;
		double runtime;
		vector< double > timestamp;

		cout << "==================================================" << endl;
		start = chrono::system_clock::now();

		cout << "Starting SDDP procedure ... " << endl;
		IloInt initSampleSize = fData.numFWsample;
		IloNumArray lb(fData.dataEnv); // double array to record lowerbound
		IloNumArray ub_c(fData.dataEnv); // double array to record ub center
		IloNumArray ub_l(fData.dataEnv); // double array to record ub lower interval
		IloNumArray ub_r(fData.dataEnv); // double array to record ub upper interval

		int iteration = 0; // iteration counter

		// create an array to store the sampled paths (rhs) in the forward pass
		IloNumArray3 samplePaths(fData.dataEnv);
		// create an array to store the sampled paths (coeff) in the forward pass
		// IloNumArray3 coefSamplePaths(fData.dataEnv);
		// create an array to store the candidate solutions corrsp. to the sampled paths
		IloNumArray3 candidateSol(fData.dataEnv);

		// create an array to store scenario indices
		IloIntArray2 scenarioIndex(fData.dataEnv);

		// start the loop until some stopping creterior is hit
		bool stable = 0;		
		do
		{
			cout << "================================" << endl;
			cout << "Iteration: " << iteration + 1 << endl;
			if (iteration == 0)
			{
				bendersFlag = 1;
				impvdBendersFlag = 0;
				lagrangianFlag = 0;
				integerFlag = 0;
				// yFlag = 0;
				for ( int t = 0; t < fData.numStage; ++t )
				{
					models[t].mod.add(models[t].relaxVarX);
					models[t].mod.add(models[t].relaxVarY);
				}
			}
			if (iteration == LPITER)
			{
				bendersFlag = atoi(argv[1]);
				impvdBendersFlag = atoi(argv[2]);
				lagrangianFlag = atoi(argv[3]);
				integerFlag = atoi(argv[4]);
				// yFlag = atoi(argv[6]);
				
				if ( impvdBendersFlag + lagrangianFlag + integerFlag )
				{
					for ( int t = 0; t < fData.numStage; ++t )
					{
						models[t].relaxVarX.end();
						models[t].relaxVarY.end();

						// models[t].cplex.setParam(IloCplex::Param::MIP::Pool::Intensity, 4);
						// models[t].cplex.setParam(IloCplex::Param::MIP::Limits::Populate, CPX_BIGINT);
						// models[t].cplex.setParam(IloCplex::Param::MIP::Pool::AbsGap, 1e-2);
						// models[t].cplex.setParam(IloCplex::Param::MIP::Pool::Replace, 2);
					}
				}
			}

			//getSamplePaths(samplePaths, fData_p, unif, generator);
			getSamplePaths(samplePaths, scenarioIndex, fData_p, sampling);
			
			cout << "Sampling completed." << endl;
			cout << "---------------------" << endl;

			forward(models, fData_p, samplePaths, candidateSol, ub_c, ub_l, ub_r);

			cout << "Forward pass completed." << endl;
			cout << "---------------------" << endl;
			
			backward(models, fData_p, candidateSol, lb, bendersFlag, impvdBendersFlag, integerFlag, lagrangianFlag);

			cout << "Backward pass completed." << endl;
			

			cout << "Current bounds on objVal: [" << lb[iteration] << ", " << ub_r[iteration]  << "]" << endl;
			
			candidateSol.clear();
			samplePaths.clear();
			// coefSamplePaths.clear();

			/* Stopping criteria heuristic
			    After a certain number of iterations, if lower bound starts to stablize,
				i.e., variance of the last five lower bounds is small enough
			*/
			/*
			if ( iteration > 12 )
			{
				vector<float> recentLB;
				for ( int i = 1; i < 10; ++i )
					recentLB.push_back(lb[iteration-i]);

				double stdReLB = std_dev(recentLB);
				if ( stdReLB < TOLOPT )
				{
					if ( fData.numFWsample == 100 )
					{
						cout << "Lower bound has stablized." << endl;
						break;
					}
					else
					{
						fData.numFWsample = 100;
					}				
				}
				else
				{
					fData.numFWsample = initSampleSize;
				}
			}
			*/
			iteration += 1;
			end = chrono::system_clock::now();
			elapsed_seconds = end - start;
			runtime = elapsed_seconds.count();
			timestamp.push_back(runtime);

		} while ( (iteration < MAXITER) && (runtime < 18000) );

		cout << "================================" << endl;
		cout << "Increasing sample size." << endl;
		fData.numFWsample = 300;
		getSamplePaths(samplePaths, scenarioIndex, fData_p, 1);	
		forward(models, fData_p, samplePaths, candidateSol, ub_c, ub_l, ub_r);
		lb.add(lb[iteration-1]);
		cout << "================================" << endl;
		cout << "Bounds on objVal: [" << lb[iteration] << ", " << ub_r[iteration] << "]" << endl;
		printf("Gap: %.2f%%. \n", (ub_r[iteration] - lb[iteration])/lb[iteration] * 100);
		
		candidateSol.clear();
		samplePaths.clear();
		// coefSamplePaths.clear();

		end = chrono::system_clock::now();
		elapsed_seconds = end - start;
		runtime = elapsed_seconds.count();
		printf("Total running time %.2f seconds.\n", runtime);
		
		// ofstream summary ("summary_118bus.csv", ios::out | ios::app);
		// if ( summary.is_open() )
		// {
		// 	summary << fData.numStage << "," << fData.numScen[1] << "," <<
		// 			   impvdBendersFlag << "," << integerFlag << "," <<
		// 			   lagrangianFlag << "," << (LEVEL_SWITCH && lagrangianFlag) << "," <<
		// 			   lb[LPITER - 1] << "," << lb[iteration] << "," << 
		// 			   ub_r[iteration] << "," << 
		// 			   (ub_r[iteration] - lb[LPITER-1])/lb[iteration] << ","  <<
		// 			   (ub_r[iteration] - lb[iteration])/lb[iteration] << ","  <<
		// 			   int(runtime) << endl;
		// }

		// ofstream summary_time ("summary_time.csv", ios::out | ios::app);
		// if ( summary_time.is_open() )
		// {
		// 	summary_time << fData.numStage << "," << fData.numScen[1] << "," <<
		// 			   impvdBendersFlag << "," << integerFlag << "," <<
		// 			   lagrangianFlag << "," << LEVEL_SWITCH;
		// 	for (auto i = timestamp.begin(); i != timestamp.end(); ++i)
		// 		summary_time << "," << *i ;
		// 	summary_time << endl;
		// }

		// ofstream summary_bound ("summary_118bus_bound.csv", ios::out | ios::app);
		// if ( summary_bound.is_open() )
		// {
		// 	summary_bound << fData.numStage << "," << fData.numScen[1] << "," <<
		// 			   impvdBendersFlag << "," << integerFlag << "," <<
		// 			   lagrangianFlag << "," << (LEVEL_SWITCH && lagrangianFlag) ;
		// 	for (int i = 0; i < lb.getSize(); ++i)
		// 		summary_bound << "," << lb[i];
		// 	summary_bound << endl;
		// }

		// free memory
		delete [] models;
		lb.end();
		ub_l.end();
		ub_c.end();
		ub_r.end();
		samplePaths.end();
		candidateSol.end();
	}
	catch (const IloException& e)
	{
		cerr << "Exception caught: " << e << endl;
  	}
   	catch (char const* status)
   	{
   		cerr << "Exception caught: " << status;
   	}
   	catch (...)
   	{
      	cerr << "Unknown exception caught!" << endl;
   	}

   	return 0;
}


// TODO:
// 1. generate data input (DONE)
// 2. sparse representation of input data
// 3. Subgradient algorithm is not working properly [DONE]
