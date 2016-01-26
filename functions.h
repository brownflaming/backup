#ifndef FUNCTIONS_H_INCLUDED
#define FUNCTIONS_H_INCLUDED

#include "global.h"

template <class T>
inline std::string toString( const T & value)
{
    std::ostringstream streamOut;
    streamOut << value;
    return streamOut.str();
}

inline std::vector<int> getDigit (const std::string & str)
{
	std::vector<int> vec;
	size_t found = str.find_first_of("01");
	while ( found != std::string::npos )
	{
		vec.push_back( int(str[found] - '0') );
		found = str.find_first_of("01",found+1);
	}
	return vec;
}

template <class T>
inline void readArray (T & target, const char * fileName)
{
  	std::ifstream data(fileName);
	if ( !data ) throw(-1);
	data >> target;
	data.close();
}

void readData (formatData * fData_p);

void buildModel (Model * models, formatData * fData_p);

void getSamplePaths (IloNumArray3& samplePaths, IloNumArray3 & coefSamplePaths, formatData * fData_p);

void forward (Model * models, formatData * fData_p, const IloNumArray3 samplePaths, const IloNumArray3 coefSamplePaths,
		IloNumArray3 & candidateSol, IloNumArray & ub, IloNumArray & ub_l, IloNumArray & ub_r);

void backward (Model * models, formatData * fData_p, const IloNumArray3 candidateSol,
		IloNumArray & lb, std::unordered_set<std::string> & masterSol,
		const bool bendersFlag, const bool impvdBendersFlag,
		const bool integerFlag, const bool lagrangianFlag );

void LGupdate( IloModel & modelLGR, IloCplex & cplexLGR, IloObjective & objLGR,
		IloNumArray & multiplier, IloNumArray & scenObj,
		const Model model, const std::vector<int> state );
	
double avg ( std::vector<float> & v );

double std_dev ( std::vector<float> & v );

void usage (char *progname);


#endif // FUNCTIONS_H_INCLUDED
