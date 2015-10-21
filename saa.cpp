#include <ilcplex/ilocplex.h>
#include <cmath>

ILOSTLBEGIN

int main ()
{
	// load data //
	
	

	// construct saa problem
	IloEnv env;
	try 
	{
		IloModel mod(env);
		
		// initialize decision variables
		int n, i, k;
		// expansion decision
		IloIntVarArray2 x(env, NODE);
		for ( n = 0; n < NODE; ++n )
			x[n] = IloIntVarArray(env, GENERATOR);

		// generation decision
		IloNumVarArray3 y(env, NODE);
		for ( n = 0; n < NODE; ++n )
		{
			y[n] = IloNumVarArray2(env, GENERATOR);
			for ( i = 0; i < GENERATOR; ++i )
				y[n][i] = IloNumVarArray(env, SUBPERIOD);
		}

		// penalties
		IloNumVarArray2 z(env, NODE);
		for ( n = 0; n < NODE; ++n )
			z[n] = IloNumVarArray(env, SUBPERIOD);

		// construct objective function
		IloObjective obj = IloMinimize(env);

		// initialize constraints
		IloRangeArray constr(env);

		// construct constraints "\sum_{m \in P(n)} x_m  >= A_n * y_nk for all n and k"
		// construct constraints "\sum_{m \in P(n)} x_m  <= UMAX       for all n in S_T"
		// construct constraints "\sum_{i \in I} y_{nki} == d_{nk}     for all n, k"




	}

}
