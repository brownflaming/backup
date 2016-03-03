#include <ilcplex/ilocplex.h>
#include <cmath>
#include <ctime>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <chrono>
#include <random>
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
		if ( argc != 5 && argc != 6 )
		{
			usage (argv[0]);
			throw (-1);
		}

		bool cutFlag[4];
		cutFlag[0] = atoi(argv[1]);  // Benders' cut
		cutFlag[1] = atoi(argv[2]);  // Improved Benders' cut
		cutFlag[2] = atoi(argv[3]);  // Lagrangian cut
		cutFlag[3] = atoi(argv[4]);  // L-shaped cut
		//if ( impvdBendersFlag + lagrangianFlag )
		// 	integerFlag = 0;

		unsigned long long seed;
		if ( argc == 6 )
			seed = atoi(argv[5]);
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
		cout << "subsample set size in the forward pass: " << fData.numFWsample << endl;
		cout << "==================================================" << endl;

		// construct models from fData
		model * models = new model[fData.numStage];
		buildModel(models, fData_p);
		cout << "Model construction completed." << endl;
		cout << "==================================================" << endl;

		// start siddp method
		// clock_t startTime = clock();
		chrono::time_point<chrono::system_clock> start, end;
		chrono::duration<double> elapsed_seconds;
		double runtime, runtime_fw, runtime_bw;
		start = chrono::system_clock::now();

		cout << "Starting SDDP procedure ... " << endl;
		IloInt initSampleSize = fData.numFWsample;
		IloNumArray lb(fData.dataEnv); // double array to record lowerbound
		IloNumArray ub_c(fData.dataEnv); // double array to record ub center
		IloNumArray ub_l(fData.dataEnv); // double array to record ub lower interval
		IloNumArray ub_r(fData.dataEnv); // double array to record ub upper interval

		IloInt iteration = 0; // iteration counter

		forwardPath samplePaths;
		// forwardPath * samplePaths_p = & samplePaths;
		IloNumArray3 candidateSol(fData.dataEnv);

		// start the loop until some stopping creterior is hit		
		do
		{
			iteration += 1;
			cout << "********************************" << endl;
			cout << "********************************" << endl;			
			cout << "Iteration: " << iteration << endl;
			getSamplePaths(samplePaths, fData_p);	
			cout << "Forward sample paths obtained." << endl;	
			cout << "================================" << endl;
			forward(models, fData_p, samplePaths, candidateSol, ub_c, ub_l, ub_r);
			cout << "Forward pass completed." << endl;
			cout << "================================" << endl;
			// cout << candidateSol << endl;
			// cin.get();
			
			/*
			if ( (! integerFlag) && (! lagrangianFlag) && (iteration > 2) )
			{
				if ( lb[iteration-2]-lb[iteration-3] < 0.05 )
				{
					integerFlag = 1;
					impvdBendersFlag = 0;
				}
			}
			cout << "L-shaped cuts: " << integerFlag << endl;
			*/

			backward(models, fData_p, candidateSol, lb, cutFlag, iteration);
			cout << "Backward pass completed." << endl;
			cout << "================================" << endl;

			cout << "lb: "   << lb   << endl;
			cout << "ub_c: " << ub_c << endl;
			cout << "ub_r: " << ub_r << endl;
			cout << "number of iterations finished: " << iteration << endl;

			candidateSol.clear();
			samplePaths = {};


			/* Stopping criteria heuristic
			    After a certain number of iterations, if lower bound starts to stablize,
				i.e., variance of the last five lower bounds is small enough
			*/
			
			if ( iteration >= 5 )
			{
				vector<float> recentLB;
				for ( int i = 1; i < 4; ++i )
					recentLB.push_back(lb[iteration-i]);

				double stdReLB = std_dev(recentLB);
				if ( stdReLB < TOLOPT )
				{
					if ( fData.numFWsample == 50 )
					{
						cout << "Lower bound has stablized." << endl;
						break;
					}
					else
						fData.numFWsample = 50;			
				}
				else
					fData.numFWsample = initSampleSize;
			}
			
			end = chrono::system_clock::now();
			elapsed_seconds = end - start;
			runtime = elapsed_seconds.count();

		} while ( (iteration < MAXITER) && (runtime < 18000) );

		fData.numFWsample = 100;
		getSamplePaths(samplePaths, fData_p);
		forward(models, fData_p, samplePaths, candidateSol, ub_c, ub_l, ub_r);

		end = chrono::system_clock::now();
		elapsed_seconds = end - start;
		runtime = elapsed_seconds.count();
		printf("Total running time %.2f seconds.\n", runtime);

		ofstream output ("0215_result.txt", ios::out | ios::app);
		if ( output.is_open() )
		{
			output << "==================================================" << endl;
			output << "==================================================" << endl;
			output << "time horizon: " << fData.numStage << endl;
			output << "Benders cut: " << cutFlag[0] << endl;
			output << "Improved Benders cut: " << cutFlag[1] << endl;
			output << "Lagrangian cut: " << cutFlag[2] << endl;
			output << "Integer Optimality cut: " << cutFlag[3] << endl;
			output << "FW sample paths: " << initSampleSize << endl;
			output << "total iterations: " << iteration << endl;
			output << "total time elapsed: " << runtime << " seconds." << endl;
			output << "lower bounds improvement: " << lb   << endl;
			output << "upper bounds improvement: " << ub_c << endl;
			output << "left 95\% CI for the upper bound: "  << ub_l << endl;
			output << "right 95\% CI for the upper bound: " << ub_r << endl;
		}

		ofstream table ("0226.csv", ios::out | ios::app);
		if ( table.is_open() )
		{
			table << initSampleSize << ", " << 
				lb[iteration-1] << ", " <<
				iteration << ", " <<
				ub_l[iteration] << ", " <<
				ub_r[iteration] << ", " <<
				(ub_r[iteration] - lb[iteration-1])/ub_r[iteration] << ", " <<
				runtime / iteration << ", " <<
				runtime << endl;
		}

		// free memory
		delete [] models;
		lb.end();
		ub_l.end();
		ub_c.end();
		ub_r.end();
		candidateSol.end();
		samplePaths = {};

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