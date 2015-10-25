#include <ilcplex/ilocplex.h>
#include <cmath>

ILOSTLBEGIN

typedef IloArray<IloIntVarArray> IloIntVarArray2;
typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumArray3> IloNumArray4;
typedef IloArray<IloNumArray2> IloNumArray3;

template <class T>
inline void readArray (T & target, const char * fileName)
{
	ifstream data(fileName);
	if ( !data ) throw(-1);
	data >> target;
	data.close();
}

int main ()
{

	// construct saa problem
	IloEnv env;
	try 
	{
		// load data //
		cout << "Reading data from files..." << endl;
		IloInt numStage, numGen, numSub, numChi, numNode;
		IloNumArray2 xCoef(env);     // [numStage][numGen];
		IloNumArray4 yCoef(env);     // [numStage][numScen][numSub][numGen]
		IloNumArray2 zCoef(env);     // [numStage][numSub]
		IloNumArray maxOutput(env);  // [numGen]
		IloNumArray maxUnit(env);    // [numGen]
		IloNumArray3 demand(env);    // [T][numChi][numSub]

		readArray<IloInt> (numStage, "data/numStage.dat");
		readArray<IloInt> (numNode, "data/numNode.dat");
		readArray<IloInt> (numGen, "data/numGen.dat");
		readArray<IloInt> (numSub, "data/numSub.dat");
		readArray<IloInt> (numChi, "data/numChi.dat");
		readArray<IloNumArray2> (xCoef, "data/xCoef.dat");
		readArray<IloNumArray4> (yCoef, "data/yCoef.dat");
		readArray<IloNumArray2> (zCoef, "data/zCoef.dat");
		readArray<IloNumArray> (maxOutput, "data/maxOutput.dat");
		readArray<IloNumArray> (maxUnit, "data/maxUnit.dat");
		readArray<IloNumArray3> (demand, "data/demand.dat");
		
		cout << "All data loaded. Constructing model..." << endl;

		IloModel mod(env);
		
		// initialize decision variables
		int n, i, k, t;
		IloIntVarArray2 x(env, numNode);  // expansion decision
		IloNumVarArray3 y(env, numNode);  // generation decision
		IloNumVarArray2 z(env, numNode);  // penalties
		for ( n = 0; n < numNode; ++n )
		{
			x[n] = IloIntVarArray(env, numGen, 0.0, IloInfinity);
			mod.add(x[n]);
			y[n] = IloNumVarArray2(env, numSub);
			for ( k = 0; k < numSub; ++k )
			{
				y[n][k] = IloNumVarArray(env, numGen, 0.0, IloInfinity, ILOFLOAT);
				mod.add(y[n][k]);
			}
			z[n] = IloNumVarArray(env, numSub, 0.0, IloInfinity, ILOFLOAT);
			mod.add(z[n]);
		}

		cout << "Decision variables constructed." << endl;

		// construct objective function
		IloObjective obj = IloMinimize(env);
		IloExpr objExpr(env);
		for ( t = 0; t < numStage; ++t )
		{
			for ( n = (pow(numChi, t) - 1) / (numChi - 1) ; n < (pow(numChi, t+1) - 1) / (numChi - 1); ++n )
			{
				objExpr += IloScalProd(xCoef[t], x[n]);
				if ( t == 0 )
				{
					for ( k = 0; k < numSub; ++k )
						objExpr += IloScalProd(yCoef[t][0][k], y[n][k])
				}
				else
				{
					for ( k = 0; k < numSub; ++k )
						objExpr += IloScalProd(yCoef[t][(n-1)%numChi][k], y[n][k]);
				}
				objExpr += IloScalProd(zCoef[t], z[n]);
			}
		}
		obj.setExpr(objExpr);
		mod.add(obj);

		cout << "Objective function added to model." << endl;

		// construct constraints
		IloRangeArray constr(env);
		for ( t = 0; t < numStage; ++t )
		{
			for ( n = (pow(numChi, t) - 1) / (numChi - 1) ; n < (pow(numChi, t+1) - 1) / (numChi - 1); ++n )
			{
				// find the path to n
				int path2n[t+1];
				int parent = n;
				for ( int s = t; s >= 0; --s )
				{
					path2n[s] = parent;
					parent = (parent - 1) / numChi;
				}
				// constraints "\sum_{m \in P(n)} x_m  >= A_n * y_nk for all n and k"
				for ( i = 0; i < numGen; ++i ) // for each type of technology
				{
					for  ( k = 0; k < numSub; ++k ) // for each subperiod
					{
						IloExpr expr(env);
						for ( int m = 0; m < sizeof(path2n)/sizeof(path2n[0]); ++m )
							expr += x[path2n[m]][i];
						expr -= (1.0 / maxOutput[i]) * y[n][k][i];
						constr.add(expr >= 0);
						expr.end();
					}
				}
				
				// constraints "\sum_{m \in P(n)} x_m  <= UMAX for all n in S_T"
				if ( t == (numStage - 1) )
				{
					for ( i = 0; i < numGen; ++i )
					{
						IloExpr expr(env);
						for ( int m = 0; m < sizeof(path2n)/sizeof(path2n[0]); ++m )
							expr += x[path2n[m]][i];
						constr.add(expr <= maxUnit[i]);
						expr.end();
					}
				}
			
				// constraints "\sum_{i \in I} y_{nki} == d_{nk} for all n, k"
				if ( t > 0 )
				{	
					for ( k = 0; k < numSub; ++k )
					{
						IloExpr expr(env);
						for ( i = 0; i < numGen; ++i )
							expr += y[n][k][i];
						expr += z[n][k];
						constr.add(expr == demand[t][(n-1)%numChi][k]);
						expr.end();
					}
				}
			} //end of loop over n
		}// end of loop over t

		mod.add(constr);

		cout << "Constraints added to model." << endl;

		// create cplex algorithm
		IloCplex cplex(mod);
		cplex.setParam(IloCplex::TiLim, 7200);
		
		// write model to file
		char fileName[100];
		sprintf(fileName, "tree_model.lp");
		cplex.exportModel(fileName);

		cout << "Model written in file tree_model.lp." << endl;
		
		cplex.solve();

		cout << "=============================================" << endl;
		IloAlgorithm::Status solStatus;
		solStatus = cplex.getStatus();
		cout << "solution status: " << solStatus << endl;	
		if ( solStatus == IloAlgorithm::Optimal )
			cout << "Optimal value: " << cplex.getObjValue() << endl;
	}
	catch (const IloException & e)
	{
		cerr << "Exception caught: " << e << endl;
	}
	catch (char const* status)
	{
		cerr << "Exception caught: " << status << endl;
	}
	catch (...)
	{
		cerr << "Unknown exception caught!" << endl;
	}

	return 0;
}
