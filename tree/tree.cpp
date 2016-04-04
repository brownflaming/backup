#include <ilcplex/ilocplex.h>
#include <cmath>

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumVarArray> IloNumVarArray2;
typedef IloArray<IloNumVarArray2> IloNumVarArray3;
typedef IloArray<IloNumArray2> IloNumArray3;
typedef IloArray<IloNumArray3> IloNumArray4;


template <class T>
inline void readArray (T & target, const char * fileName)
{
	ifstream input(fileName);
	if ( !input ) throw(-1);
	input >> target;
	input.close();
}

int main (int argc, char *argv[])
{
	// construct tree problem
	IloEnv env;
	try
	{
		int LP = atoi(argv[1]);

		IloInt numStage;
		IloInt numBranch;
		readArray<IloInt> (numStage, "data/stage.dat");
		readArray<IloInt> (numBranch, "data/branch.dat");
		cout << "test" << endl;
		IloInt ODI = 12;
		IloInt CLASS = 6;
		IloInt numNode = (pow(numBranch, numStage) - 1)/(numBranch - 1);

		IloNumArray PRICE1 = IloNumArray(env, 6, 500.0, 340.0, 200.0, 160.0, 130.0, 100.0);
		IloNumArray PRICE2 = IloNumArray(env, 6, 800.0, 540.0, 320.0, 260.0, 210.0, 160.0);
		IloNumArray CLNRATE = IloNumArray(env, 6, 0.1, 0.1, 0.05, 0.05, 0.0, 0.0);
		IloNumArray SEAT = IloNumArray(env, 2, 24, 216);
		IloInt K = (SEAT[0] + SEAT[1]) * 5;

		IloNumArray4 DEMAND(env);  // [t][k][i][j]
		readArray<IloNumArray4> (DEMAND, "data/demand.dat");

		IloModel mod(env);
		int n, i, j, k, t;

		// initialize decision variables
		IloNumVarArray3 B(env, numNode);
		IloNumVarArray3 C(env, numNode);
		// IloNumVarArray3 P(env, numNode);
		IloNumVarArray3 b(env, numNode);
		IloNumVarArray3 c(env, numNode);
		// IloNumVarArray3 z(env, numNode);
		// IloNumVarArray3 zp(env, numNode);
		// IloNumVarArray3 zd(env, numNode);
		char varName[100];
		for ( n = 0; n < numNode; ++n )
		{
			B[n] = IloNumVarArray2(env, ODI);
			C[n] = IloNumVarArray2(env, ODI);
			// P[n] = IloNumVarArray2(env, ODI);
			b[n] = IloNumVarArray2(env, ODI);
			c[n] = IloNumVarArray2(env, ODI);
			// z[n] = IloNumVarArray2(env, ODI);
			// zp[n] = IloNumVarArray2(env, ODI);
			// zd[n] = IloNumVarArray2(env, ODI);

			for ( i = 0; i < ODI; ++i )
			{
				if ( ! LP )
				{
					B[n][i] = IloNumVarArray(env, CLASS, 0, 511, ILOINT);
					C[n][i] = IloNumVarArray(env, CLASS, 0, 511, ILOINT);
					// P[n][i] = IloNumVarArray(env, CLASS, 0, 511, ILOINT);
					b[n][i] = IloNumVarArray(env, CLASS, 0, 511, ILOINT);
					c[n][i] = IloNumVarArray(env, CLASS, 0, 511, ILOINT);
					// z[n][i] = IloNumVarArray(env, CLASS, 0, 1, ILOINT);
				}
				else
				{
					B[n][i] = IloNumVarArray(env, CLASS, 0, 511);
					C[n][i] = IloNumVarArray(env, CLASS, 0, 511);
					// P[n][i] = IloNumVarArray(env, CLASS, 0, 511);
					b[n][i] = IloNumVarArray(env, CLASS, 0, 511);
					c[n][i] = IloNumVarArray(env, CLASS, 0, 511);
					// z[n][i] = IloNumVarArray(env, CLASS, 0, 1);	
				}
				// zp[n][i] = IloNumVarArray(env, CLASS, 0, IloInfinity);
				// zd[n][i] = IloNumVarArray(env, CLASS, 0, IloInfinity);

				for ( j = 0; j < CLASS; ++j )
				{ 
					sprintf(varName, "B_%d_%d_%d", n, i, j);
					B[n][i][j].setName(varName);
					sprintf(varName, "C_%d_%d_%d", n, i, j);
					C[n][i][j].setName(varName);
					// sprintf(varName, "P_%d_%d_%d", n, i, j);
					// P[n][i][j].setName(varName);
					sprintf(varName, "b_%d_%d_%d", n, i, j);
					b[n][i][j].setName(varName);
					sprintf(varName, "c_%d_%d_%d", n, i, j);
					c[n][i][j].setName(varName);
					// sprintf(varName, "z_%d_%d_%d", n, i, j);
					// z[n][i][j].setName(varName);
					// sprintf(varName, "zp_%d_%d_%d", n, i, j);
					// zp[n][i][j].setName(varName);
					// sprintf(varName, "zd_%d_%d_%d", n, i, j);
					// zd[n][i][j].setName(varName);
				}
				mod.add(B[n][i]);
				mod.add(C[n][i]);
				// mod.add(P[n][i]);
				mod.add(b[n][i]);
				mod.add(c[n][i]);
				// mod.add(z[n][i]);
				// mod.add(zp[n][i]);
				// mod.add(zd[n][i]);
			}
		}
		cout << "Decision variables constructed." << endl;

		// construct objective function
		IloObjective obj = IloMaximize(env);
		IloExpr objExpr(env);
		IloExpr term(env);
		for ( n = 1; n < numNode; ++n )
		{
			int stage = ceil(log(numBranch + (numBranch - 1) * n)/log(numBranch)) - 1;
			cout << stage << endl;
			for ( i = 0; i < ODI; ++i )
			{
				if ( i < 6 )
				{
					term += IloScalProd(PRICE1, b[n][i]);
					term -= IloScalProd(PRICE1, c[n][i]);
					term = term / pow(numBranch, stage);
					objExpr += term;
					term.clear();
				}
				else
				{
					term += IloScalProd(PRICE2, b[n][i]);
					term -= IloScalProd(PRICE2, c[n][i]);
					term = term / pow(numBranch, stage);
					objExpr += term;
					term.clear();
				}
			}
			obj.setExpr(objExpr);
			mod.add(obj);
		}
		cout << "Objective function constructed." << endl;

		// construct constraints
		IloRangeArray constr(env);
		IloExpr expr(env);
		for ( i = 0; i < ODI; ++i )
			for ( j = 0; j < CLASS; j++ )
			{
				expr = B[0][i][j];
				constr.add(expr == 0);
				expr.clear();
				expr = C[0][i][j];
				constr.add(expr == 0);
				expr.clear();
			}

		for ( n = 0; n < numNode; ++n )
			for ( i = 0; i < ODI; ++i )
				for ( j = 0; j < 2; ++j )
					constr.add(B[n][i][j] <= 63);
		

		for ( n = 1; n < numNode; ++n )
		{
			cout << n << ", " << numNode << endl;
			int stage = ceil(log(numBranch + (numBranch - 1) * n)/log(numBranch)) - 1;
			cout << stage << endl;
			int p = (n-1)/numBranch;
			for ( i = 0; i < ODI; ++i )
			{
				for ( j = 0; j < CLASS; ++j )
				{
					expr = B[n][i][j] - B[p][i][j] - b[n][i][j];
					constr.add(expr == 0);
					expr.clear();

					expr = C[n][i][j] - C[p][i][j] - c[n][i][j];
					constr.add(expr == 0);
					expr.clear();

					// expr = b[n][i][j] + zd[n][i][j];
					expr = b[n][i][j];
					constr.add(expr <= DEMAND[stage][(n-1)%numBranch][i][j]);
					expr.clear();

					// expr = B[n][i][j] + zp[n][i][j] - P[p][i][j] / (1.0 - CLNRATE[j]) - 0.5;
					// constr.add(expr <= 0);
					// expr.clear();
					// expr = B[n][i][j] + zp[n][i][j] - P[p][i][j] / (1.0 - CLNRATE[j]) + 0.5;
					// constr.add(expr >= 0);
					// expr.clear();

					expr = C[n][i][j] - CLNRATE[j] * B[n][i][j] - 0.5;
					constr.add(expr <= 0);
					expr.clear();
					expr = C[n][i][j] - CLNRATE[j] * B[n][i][j] + 0.5;
					constr.add(expr >= 0);
					expr.clear();

					// expr = zp[n][i][j] + K * z[n][i][j];
					// constr.add(expr <= K);
					// expr.clear();

					// expr = zd[n][i][j] - K * z[n][i][j];
					// constr.add(expr <= 0);
					// expr.clear();

				}
			}
		}

		int begin = (pow(numBranch, numStage - 1) - 1) / (numBranch - 1);

		for ( n = begin; n < numNode; ++n )
		{
			// business class

			expr = B[n][0][0] + B[n][0][1] + B[n][6][0] + B[n][6][1] +  B[n][8][0] + B[n][8][1]
				 -(C[n][0][0] + C[n][0][1] + C[n][6][0] + C[n][6][1] +  C[n][8][0] + C[n][8][1]);
			constr.add(expr <=  SEAT[0]);
			expr.clear();

			expr = B[n][1][0] + B[n][1][1] + B[n][7][0] + B[n][7][1] +	B[n][9][0] + B[n][9][1]
				 -(C[n][1][0] + C[n][1][1] + C[n][7][0] + C[n][7][1] +	C[n][9][0] + C[n][9][1]);
			constr.add(expr <=  SEAT[0]);
			expr.clear();

			expr = B[n][2][0] + B[n][2][1] + B[n][7][0] + B[n][7][1] + B[n][10][0] + B[n][10][1]
				 -(C[n][2][0] + C[n][2][1] + C[n][7][0] + C[n][7][1] + C[n][10][0] + C[n][10][1]);
			constr.add(expr <=  SEAT[0]);
			expr.clear();

			expr = B[n][3][0] + B[n][3][1] + B[n][6][0] + B[n][6][1] + B[n][11][0] + B[n][11][1]
				 -(C[n][3][0] + C[n][3][1] + C[n][6][0] + C[n][6][1] + C[n][11][0] + C[n][11][1]);
			constr.add(expr <=  SEAT[0]);
			expr.clear();

			expr = B[n][4][0] + B[n][4][1] + B[n][9][0] + B[n][9][1] + B[n][11][0] + B[n][11][1]
				 -(C[n][4][0] + C[n][4][1] + C[n][9][0] + C[n][9][1] + C[n][11][0] + C[n][11][1]);
			constr.add(expr <=  SEAT[0]);
			expr.clear();

			expr = B[n][5][0] + B[n][5][1] + B[n][8][0] + B[n][8][1] + B[n][10][0] + B[n][10][1]
				 -(C[n][5][0] + C[n][5][1] + C[n][8][0] + C[n][8][1] + C[n][10][0] + C[n][10][1]);
			constr.add(expr <=  SEAT[0]);
			expr.clear();

			// economic class

			expr = B[n][0][2] + B[n][0][3] + B[n][0][4] + B[n][0][5] + 
				   B[n][6][2] + B[n][6][3] + B[n][6][4] + B[n][6][5] +
				   B[n][8][2] + B[n][8][3] + B[n][8][4] + B[n][8][5] -
				  (C[n][0][2] + C[n][0][3] + C[n][0][4] + C[n][0][5] + 
				   C[n][6][2] + C[n][6][3] + C[n][6][4] + C[n][6][5] +
				   C[n][8][2] + C[n][8][3] + C[n][8][4] + C[n][8][5]);
			constr.add(expr <=  SEAT[1]);
			expr.clear();

			expr = B[n][1][2] + B[n][1][3] + B[n][1][4] + B[n][1][5] + 
				   B[n][7][2] + B[n][7][3] + B[n][7][4] + B[n][7][5] +
				   B[n][9][2] + B[n][9][3] + B[n][9][4] + B[n][9][5] -
				  (C[n][1][2] + C[n][1][3] + C[n][1][4] + C[n][1][5] + 
				   C[n][7][2] + C[n][7][3] + C[n][7][4] + C[n][7][5] +
				   C[n][9][2] + C[n][9][3] + C[n][9][4] + C[n][9][5]);
			constr.add(expr <=  SEAT[1]);
			expr.clear();
			
			expr = B[n][2][2] + B[n][2][3] + B[n][2][4] + B[n][2][5] + 
				   B[n][7][2] + B[n][7][3] + B[n][7][4] + B[n][7][5] +
				   B[n][10][2] + B[n][10][3] + B[n][10][4] + B[n][10][5] -
				  (C[n][2][2] + C[n][2][3] + C[n][2][4] + C[n][2][5] + 
				   C[n][7][2] + C[n][7][3] + C[n][7][4] + C[n][7][5] +
				   C[n][10][2] + C[n][10][3] + C[n][10][4] + C[n][10][5]);
			constr.add(expr <=  SEAT[1]);

			expr = B[n][3][2] + B[n][3][3] + B[n][3][4] + B[n][3][5] + 
				   B[n][6][2] + B[n][6][3] + B[n][6][4] + B[n][6][5] +
				   B[n][11][2] + B[n][11][3] + B[n][11][4] + B[n][11][5] -
				  (C[n][3][2] + C[n][3][3] + C[n][3][4] + C[n][3][5] + 
				   C[n][6][2] + C[n][6][3] + C[n][6][4] + C[n][6][5] +
				   C[n][11][2] + C[n][11][3] + C[n][11][4] + C[n][11][5]);
			constr.add(expr <=  SEAT[1]);
			expr.clear();

			expr = B[n][4][2] + B[n][4][3] + B[n][4][4] + B[n][4][5] + 
				   B[n][9][2] + B[n][9][3] + B[n][9][4] + B[n][9][5] +
				   B[n][11][2] + B[n][11][3] + B[n][11][4] + B[n][11][5] -
				  (C[n][4][2] + C[n][4][3] + C[n][4][4] + C[n][4][5] + 
				   C[n][9][2] + C[n][9][3] + C[n][9][4] + C[n][9][5] +
				   C[n][11][2] + C[n][11][3] + C[n][11][4] + C[n][11][5]);
			constr.add(expr <=  SEAT[1]);
			expr.clear();

			expr = B[n][5][2] + B[n][5][3] + B[n][5][4] + B[n][5][5] + 
				   B[n][8][2] + B[n][8][3] + B[n][8][4] + B[n][8][5] +
				   B[n][10][2] + B[n][10][3] + B[n][10][4] + B[n][10][5] -
				  (C[n][5][2] + C[n][5][3] + C[n][5][4] + C[n][5][5] + 
				   C[n][8][2] + C[n][8][3] + C[n][8][4] + C[n][8][5] +
				   C[n][10][2] + C[n][10][3] + C[n][10][4] + C[n][10][5]);
			constr.add(expr <=  SEAT[1]);
			expr.clear();
		}

		// int begin = (pow(numBranch, numStage - 2) - 1) / (numBranch - 1);
		// int end = (pow(numBranch, numStage - 1) - 1) / (numBranch - 1);

		// cout << begin << end << endl;
		// for ( n = begin; n < end; ++n )
		// {
		// 	// business class

		// 	expr = P[n][0][0] + P[n][0][1] +
		// 	       P[n][6][0] + P[n][6][1] + 
		// 	       P[n][8][0] + P[n][8][1];
		// 	constr.add(expr <=  SEAT[0]);
		// 	expr.clear();

		// 	expr = P[n][1][0] + P[n][1][1] +
		// 	       P[n][7][0] + P[n][7][1] +
		// 	       P[n][9][0] + P[n][9][1];
		// 	constr.add(expr <=  SEAT[0]);
		// 	expr.clear();

		// 	expr = P[n][2][0] + P[n][2][1] +
		// 	       P[n][7][0] + P[n][7][1] +
		// 	       P[n][10][0] + P[n][10][1];
		// 	constr.add(expr <=  SEAT[0]);
		// 	expr.clear();

		// 	expr = P[n][3][0] + P[n][3][1] +
		// 	       P[n][6][0] + P[n][6][1] +
		// 	       P[n][11][0] + P[n][11][1];
		// 	constr.add(expr <=  SEAT[0]);
		// 	expr.clear();

		// 	expr = P[n][4][0] + P[n][4][1] +
		// 	       P[n][9][0] + P[n][9][1] +
		// 	       P[n][11][0] + P[n][11][1];
		// 	constr.add(expr <=  SEAT[0]);
		// 	expr.clear();

		// 	expr = P[n][5][0] + P[n][5][1] +
		// 	       P[n][8][0] + P[n][8][1] +
		// 	       P[n][10][0] + P[n][10][1];
		// 	constr.add(expr <=  SEAT[0]);
		// 	expr.clear();

		// 	// economic class

		// 	expr = P[n][0][2] + P[n][0][3] + P[n][0][4] + P[n][0][5] + 
		// 		   P[n][6][2] + P[n][6][3] + P[n][6][4] + P[n][6][5] +
		// 		   P[n][8][2] + P[n][8][3] + P[n][8][4] + P[n][8][5];
		// 	constr.add(expr <=  SEAT[1]);
		// 	expr.clear();

		// 	expr = P[n][1][2] + P[n][1][3] + P[n][1][4] + P[n][1][5] + 
		// 		   P[n][7][2] + P[n][7][3] + P[n][7][4] + P[n][7][5] +
		// 		   P[n][9][2] + P[n][9][3] + P[n][9][4] + P[n][9][5];
		// 	constr.add(expr <=  SEAT[1]);
		// 	expr.clear();
			
		// 	expr = P[n][2][2] + P[n][2][3] + P[n][2][4] + P[n][2][5] + 
		// 		   P[n][7][2] + P[n][7][3] + P[n][7][4] + P[n][7][5] +
		// 		   P[n][10][2] + P[n][10][3] + P[n][10][4] + P[n][10][5];
		// 	constr.add(expr <=  SEAT[1]);

		// 	expr = P[n][3][2] + P[n][3][3] + P[n][3][4] + P[n][3][5] + 
		// 		   P[n][6][2] + P[n][6][3] + P[n][6][4] + P[n][6][5] +
		// 		   P[n][11][2] + P[n][11][3] + P[n][11][4] + P[n][11][5];
		// 	constr.add(expr <=  SEAT[1]);
		// 	expr.clear();

		// 	expr = P[n][4][2] + P[n][4][3] + P[n][4][4] + P[n][4][5] + 
		// 		   P[n][9][2] + P[n][9][3] + P[n][9][4] + P[n][9][5] +
		// 		   P[n][11][2] + P[n][11][3] + P[n][11][4] + P[n][11][5];
		// 	constr.add(expr <=  SEAT[1]);
		// 	expr.clear();

		// 	expr = P[n][5][2] + P[n][5][3] + P[n][5][4] + P[n][5][5] + 
		// 		   P[n][8][2] + P[n][8][3] + P[n][8][4] + P[n][8][5] +
		// 		   P[n][10][2] + P[n][10][3] + P[n][10][4] + P[n][10][5];
		// 	constr.add(expr <=  SEAT[1]);
		// 	expr.clear();
		// }

		mod.add(constr);
		expr.end();


		cout << "Constraints constructed." << endl;

		// create cplex and solve
		IloCplex cplex(mod);
		cplex.exportModel("tree.lp");
		cplex.solve();

		cout << "=============================================" << endl;
		IloAlgorithm::Status solStatus;
		IloNum objVal;
		solStatus = cplex.getStatus();
		objVal = cplex.getObjValue();
		cout << "solution status: " << solStatus << endl;	
		if ( solStatus == IloAlgorithm::Optimal )
		{
			cout << "Optimal value: " << cplex.getObjValue() << endl;

			// IloNumArray Pvals(env);
			// IloNumArray Bvals(env);
			// IloNumArray Cvals(env);
			// IloNumArray bvals(env);
			// IloNumArray cvals(env);
			// for ( n = 0; n < numNode; ++n )
			// 	for ( i = 0; i < ODI; ++i )
			// 	{
			// 		cplex.getValues(Pvals, P[n][i]);
			// 		cplex.getValues(Bvals, B[n][i]);
			// 		cplex.getValues(Cvals, C[n][i]);
			// 		cplex.getValues(bvals, b[n][i]);
			// 		cplex.getValues(cvals, c[n][i]);
			// 		for ( j = 0; j < CLASS; ++j )
			// 		{
			// 			if (Pvals[j] != 0)
			// 				cout << P[n][i][j].getName() << " = " << Pvals[j] << ", ";
			// 			if (Bvals[j] != 0)
			// 				cout << B[n][i][j].getName() << " = " << Bvals[j] << ", ";
			// 			if (Cvals[j] != 0)
			// 				cout << C[n][i][j].getName() << " = " << Cvals[j] << ", ";
			// 			if (bvals[j] != 0)
			// 				cout << b[n][i][j].getName() << " = " << bvals[j] << ", ";
			// 			if (cvals[j] != 0)
			// 				cout << c[n][i][j].getName() << " = " << cvals[j];
			// 			cout << " " << endl;
			// 		}
			// 	}



			// ofstream table ("tree.csv", ios::out | ios::app);
			// if ( table.is_open() )
			// {
			// 	table << numStage << ", " <<
			// 	NT <<
			// 	objVal << ", " << 
			// 	cplex.getTime() << endl;
			// }
		}
	}
	catch (const IloException & e)
	{
		cerr << "Exception caught: " << e << endl;
	}
	catch (char const * status)
	{
		cerr << "Exception caught: " << status << endl;
	}
	catch (...)
	{
		cerr << "Unknown exception caught! " << endl;
	}

	return 0;
}