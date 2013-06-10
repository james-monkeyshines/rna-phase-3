#ifndef TREEMODELVALUE_H
#define TREEMODELVALUE_H

#include <vector>
#include <string>
using namespace std;

class TreeModelValue {
public :
    TreeModelValue() ;
    TreeModelValue( const string& stringTree,
                    const vector < double > & branchLengths,
                    const vector < double > & modelParameters,
                    double logLikelihood ) ;
    string stringTree;
    vector < double > branchLengths; // Branch lengths
    vector < double > modelParameters; // Model Parameters
    double logLikelihood; // loglikelihood of combination

    static int compareEntries( const void * v1 , const void * v2 );
};

#endif // TREEMODELVALUE_H




