#ifndef GLOBALS_H_INCLUDED
#define GLOBALS_H_INCLUDED

/*Maximum number of iterations*/
const int MAXITER=500;

/*Tolerance for optimality : stopping for SDDP iterations*/
const double TOLOPT = 1e-3;

/*Confidence interval parameter : z_alpha */
const double ZALPHA = 1.96;

/*Define a small number for numerical convenience*/
const double EPSILON = 1e-4;

/*Flag for sampling: 
 * 0 no sampling in policy construction and policy evaluation
 * 1 sampling with 1 trial solution per iteration*/
//const int Sampling=1;

/*Set the accuracy level for the continuous variable approximation*/
//const double precision = 0.1;

/* Improved optimality cut : 
 * 0 : simple cut
 * 1 : improved cut*/
//const int impFlag = 0;

// Define 3-d array
typedef IloArray<IloNumArray2> IloNumArray3;
typedef std::vector< std::vector<int> > IntMatrix;

struct formatData
{
	IloEnv dataEnv;   // An environment associated with all data
	IloInt numStage;  // T: number of stages
	IloNumArray initState; // x_0: initial state of binary variables
	IloNumArray valueLB;	// valueLB: lower bounds for cost-to-go functions
	IloInt numFWsample;  // M: number of samples used in the foward pass
	IloIntArray numScen; // N_t: number of scenarios st each stage.
	IloInt totalScen; 	// N: total # of scenarios = numScen[0] * numScen[1] * ... * numScen[T]
	IloNumArray2 xCoef; // c_t: coefficients for x at each stage
	IloNumArray2 yCoef; // d_t: coefficients for y at each stage 
	IloNumArray3 Amatrix; // A_t at each stage
	IloNumArray3 Bmatrix; // B_t at each stage
	IloNumArray3 Wmatrix; // W_t at each stage
	IloNumArray2 bRhs;  // b_t: rhs at each stage (some components are uncertain)
	IloIntArray uncertainIndex; // an array that contains the indices for uncertain parameters in b_t
	IloNumArray3 scenarios; // available scenarios at each stage (dim1: T; dim2: numScen[t])
};

struct Model
{	
	IloCplex cplex;
	IloEnv env;
	IloModel mod;
	IloNumVarArray x;		// current stage integer variables x_t
	IloNumVarArray y;		// current stage continuous variables y_t
	IloNumVarArray z;		// copy of previous integer variables z_t = x_{t-1}
	IloNumVar theta;		// future cost approximation
	IloObjective obj;		// c_tx_t + d_ty_t + theta_{t+1}
	IloRangeArray constr1; 	// A_tx_t + W_ty_t + B_tz_t <= b_t
	IloRangeArray constr2; 	// z_t = x_{t-1}
	IloRangeArray cuts;		// a range array that records added cuts
};


#endif // GLOBALS_H_INCLUDED
