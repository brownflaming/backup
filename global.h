#ifndef GLOBALS_H_INCLUDED
#define GLOBALS_H_INCLUDED

/*Maximum number of iterations*/
const int MAXITER=1000;

/*Tolerance for optimality : stopping for SDDP iterations*/
const double TOLOPT = 1;

/*Confidence interval parameter : z_alpha */
const double ZALPHA = 1.96;

/*Define a small number for numerical convenience*/
const double EPSILON = 1e-4;

// Define 3-d array
typedef IloArray<IloNumArray2> IloNumArray3;
typedef IloArray<IloNumArray3> IloNumArray4;
typedef std::vector< std::vector<int> > IntMatrix;

struct formatData
{
	IloEnv dataEnv;   // An environment associated with all data
	IloInt numStage;  // T: number of stages
	IloInt numFWsample;  // M: number of samples used in the foward pass
	IloIntArray numScen; // N_t: number of scenarios st each stage.
	IloInt totalScen; 	// N: total # of scenarios = numScen[0] * numScen[1] * ... * numScen[T]
	IloNumArray initState; // x_0: initial state of binary variables
	IloNumArray thetaLB;	// thetaLB: lower bounds for cost-to-go functions
	IloNumArray2 constrSlack;  // lower and upper bound for slakc variables s2 [2][numRows]
	
	IloNumArray2 x; // c_t: coefficients for x at each stage
	IloNumArray2 y1; // d1_t: coefficients for y1 at each stage
	IloNumArray2 y2; // d2_t: coefficients for y2 at each stage
	IloNumArray3 A; // A_t at each stage
	IloNumArray3 B; // B_t at each stage
	IloNumArray3 W1; // W1_t at each stage
	IloNumArray3 W2; // W2_t at each stage
	IloNumArray2 S; // matrix for the slack variables
	IloNumArray2 S2; // matrix for the slack variables to accommodate infeasibility
	IloNumArray2 b;  // b_t: rhs at each stage (some components are uncertain)
	
	IloIntArray uncertainData; // 8-d binary array indicating if [x, y1, y2, A, B, W1, W2, b] has uncertainty, resp.
	IloNumArray3 xScen; // scenarios for x objective coefficients [t][numScen[t]][dimX]
	IloNumArray3 y1Scen; //scenarios for y1 objective coefficients [t][numScen[t]][dimY1]
	IloNumArray3 y2Scen; //scenarios for y2 objective coefficients [t][numScen[t]][dimY2]
	IloNumArray4 AScen;  //scenarios for A matrix [t][numScen[t]][numRows][dimX]
	IloNumArray4 BScen;  //scenarios for B matrix [t][numScen[t]][3][#nonzeros]
	IloNumArray4 W1Scen; //scenarios for W1 matrix [t][numScen[t]][numRows][dimY1]
	IloNumArray4 W2Scen; //scenarios for W2 matrix [t][numScen[t]][numRows][dimY2]
	IloNumArray3 bScen;  //scenarios for b [t][numScen[t]][numRows]
};

struct forwardPath
{
	IloNumArray3 x;   // [p][t][dimX]
	IloNumArray3 y1;  // [p][t][dimY1]
	IloNumArray3 y2;  // [p][t][dimY2]
	IloNumArray4 A;   // [p][t][numRows][dimX]
	IloNumArray4 B;   // [p][t][numRows][dimZ]
	IloNumArray4 W1;  // [p][t][numRows][dimY1]
	IloNumArray4 W2;  // [p][t][numRows][dimY2]
	IloNumArray3 b;   // [p][t][numRows]
};

struct model
{	
	IloCplex cplex;
	IloEnv env;
	IloModel mod;
	IloNumVarArray x;		// current stage state variables x_t
	IloNumVarArray z;		// copy of previous state variables z_t = x_{t-1}
	IloNumVarArray y1;		// current stage integral variables y1_t
	IloNumVarArray y2;      // current stage continuous variables y2_t
	IloNumVar theta;		// future cost approximation
	IloNumVarArray s;       // slack variables
	IloNumVarArray s2;		// slack variables 2
	IloObjective obj;		// c_tx_t + d1_ty1_t + d2_ty2_t + theta_{t+1}
	IloRangeArray constr1; 	// A_tx_t + W1_ty1_t + W2_ty2_t + B_tz_t + S s_t + S2 s2_t = b_t
	IloRangeArray constr2; 	// z_t = x_{t-1}
	IloRangeArray cuts;		// a range array that records added cuts
};


#endif // GLOBALS_H_INCLUDED
