#include <ilcplex/ilocplex.h>
#include <cmath>

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumArray2> IloNumArray3;
typedef IloArray<IloIntVarArray> IloIntVarArray2;


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
		IloInt numStage, numStock, numChi, assetLimit;
		IloNum buyTran, sellTran;
		IloNumArray initState(env);
		readArray<IloInt> (numStage, "data/numStage.dat");
		readArray<IloInt> (numStock, "data/numStock.dat");
		readArray<IloInt> (numChi, "data/numChi.dat");
		readArray<IloNum> (buyTran, "data/buyTran.dat");
		readArray<IloNum> (sellTran, "data/sellTran.dat");
		readArray<IloNumArray> (initState, "data/initState.dat");
		readArray<IloInt> (assetLimit, "data/assetLimit.dat");

		IloNumArray3 scenario(env); //[numStage][numChi][numStock]
		readArray<IloNumArray3> (scenario, "data/scenario.dat");

		cout << "All data loaded. Constructing model..." << endl;

		IloInt numNode = (pow(numChi, numStage) - 1)/(numChi - 1);
		IloInt NT = pow(numChi, numStage-1);

		IloModel mod(env);
		
		// initialize decision variables
		int n, i, k, t;
		IloNumVarArray2 x(env, numNode);  // asset allocation
		IloNumVarArray2 buy(env, numNode);  // buy amount
		IloNumVarArray2 sell(env, numNode); // sell amount
		IloIntVarArray2 pos(env, numNode);  // position at each node
		char varName[100];
		for ( n = 0; n < numNode; ++n )
		{
			x[n] = IloNumVarArray(env, numStock, 0.0, IloInfinity);
			for ( i = 0; i < numStock; ++i )
			{
				sprintf(varName, "x_%d%d", n, i);
				x[n][i].setName(varName);
			}

			buy[n] = IloNumVarArray(env, numStock, 0.0, IloInfinity);
			for ( i = 0; i < numStock; ++i )
			{
				sprintf(varName, "b_%d%d", n, i);
				buy[n][i].setName(varName);
			}

			sell[n] = IloNumVarArray(env, numStock, 0.0, IloInfinity);
			for ( i = 0; i < numStock; ++i )
			{
				sprintf(varName, "s_%d%d", n, i);
				sell[n][i].setName(varName);
			}

			pos[n] = IloIntVarArray(env, numStock, 0.0, 1.0);
			for ( i = 0; i < numStock; ++i )
			{
				sprintf(varName, "p_%d%d", n, i);
				pos[n][i].setName(varName);
			}


			mod.add(x[n]);
			mod.add(buy[n]);
			mod.add(sell[n]);
			mod.add(pos[n]);
		}
		

		cout << "Decision variables constructed." << endl;

		// construct objective function
		IloObjective obj = IloMaximize(env);
		IloExpr objExpr(env);
		for ( n = numNode - NT; n < numNode; ++n )
		{
			for ( k = 0; k < numChi; ++k )
				objExpr += IloScalProd(x[n], scenario[numStage][k]);
		}
		objExpr = objExpr / (NT * numChi);
		obj.setExpr(objExpr);
		mod.add(obj);

		cout << "Objective function added to model." << endl;

		// construct constraints
		IloRangeArray constr(env);
		IloExpr expr(env);
		for ( t = 0; t < numStage; ++t )
		{
			for ( n = (pow(numChi, t) - 1) / (numChi - 1) ; n < (pow(numChi, t+1) - 1) / (numChi - 1); ++n )
			{
				int k = (n-1) % numChi;
				// constraints " x(n) = r(n)x(a(n)) + buy(n) - sell(n)"
				for ( i = 0; i < numStock; ++i )
				{
					if ( n == 0 )
						expr = x[n][i] - initState[i] - buy[n][i] + sell[n][i];
					else
						expr = x[n][i] - scenario[t][k][i] * x[(n-1)/numChi][i] - buy[n][i] + sell[n][i];
					constr.add(expr == 0);
				}
				expr.clear();

				// self-financing constraints 
				for ( i = 0; i < numStock; ++i )
				{
					expr += (1 + buyTran) * buy[n][i];
					expr -= (1 - sellTran) * sell[n][i];
				}
				constr.add(expr == 0);
				expr.clear();

				// linking constraints
				for ( i = 0; i < numStock; ++i )
				{
					expr = x[n][i] - 1000 * pos[n][i];
					constr.add(expr <= 0);
				}
				expr.clear();
				for ( i = 0; i < numStock; ++i )
					expr += pos[n][i];
				constr.add(expr <= assetLimit);
				expr.clear();

			} //end of loop over n
		}// end of loop over t
		expr.end();
		mod.add(constr);

		cout << "Constraints added to model." << endl;

		// create cplex algorithm
		IloCplex cplex(mod);
		
		// write model to file
		// char fileName[100];
		// sprintf(fileName, "tree_model.lp");
		// cplex.exportModel(fileName);
		
		cplex.solve();

		cout << "=============================================" << endl;
		IloAlgorithm::Status solStatus;
		solStatus = cplex.getStatus();
		cout << "solution status: " << solStatus << endl;	
		if ( solStatus == IloAlgorithm::Optimal )
		{
			cout << "Optimal value: " << cplex.getObjValue() << endl;
			/*
			cout << scenario << endl;
			IloNumArray val(env);
			for ( n = 0; n < numNode; ++n )
			{
				cout << "====================================" << endl;
				cout << "node " << n << endl;
				cplex.getValues(val, x[n]);
				cout << "position: " << val << endl;
				cplex.getValues(val, buy[n]);
				cout << "buy decision: " << val << endl;
				cplex.getValues(val, sell[n]);
				cout << "sell decision: " << val << endl;
			}
			*/
		}
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
