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

void buildModel (model * models, formatData * fData_p);

void getSamplePaths (forwardPath & samplePaths, formatData * fData_p);

void forward (model * models, formatData * fData_p,
	const forwardPath samplePaths, IloNumArray3 & candidateSol,
	IloNumArray & ub_c, IloNumArray & ub_l, IloNumArray & ub_r);

void backward (model * models, formatData * fData_p,
	const IloNumArray3 candidateSol, IloNumArray & lb, bool cutFlag[4], const IloInt iter);

double LGsolve (model & LGmodel, IloNumArray & dualVar,
	IloNumVarArray z, const std::vector<int> state, const double ub );

double avg ( std::vector<float> & v );

double std_dev ( std::vector<float> & v );

void usage (char *progname);


#endif // FUNCTIONS_H_INCLUDED
