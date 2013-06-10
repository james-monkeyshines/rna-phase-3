#include "Tree/HeterogeneousCluster.h"

#include "assert.h"

#include <iostream>


HeterogeneousCluster::HeterogeneousCluster( unsigned int totalNumberSpecies ) :
LengthCluster( totalNumberSpecies ){
}

HeterogeneousCluster::HeterogeneousCluster( unsigned int totalNumberSpecies,
                                            double initialLength ):
LengthCluster( totalNumberSpecies, initialLength ){
}

void HeterogeneousCluster::addParam( const Cluster& other ){
    LengthCluster::addParam(other);
    averageRate += ((HeterogeneousCluster*)&other)->averageRate;
    for (unsigned int i = 0; i < parameters.size(); ++i){
        parameters[i] += ((HeterogeneousCluster*)&other)->parameters[i];
    }
    for ( map<unsigned int, unsigned int>::const_iterator i =
               ((HeterogeneousCluster*)&other)->getModelCount().begin();
          i != ((HeterogeneousCluster*)&other)->getModelCount().end(); ++i){
        modelCount[(*i).first] += (*i).second;
    }
}    

void HeterogeneousCluster::addModelParameters( double averageSubstitutionRate,
                                               const vector<double> & param,
                                               unsigned int model ){
    assert( number == 1 );
    averageRate = averageSubstitutionRate;
    parameters = param;
    modelCount[model]=1;
}
