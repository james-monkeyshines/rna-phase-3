#ifndef HETEROGENEOUSCLUSTER_H
#define HETEROGENEOUSCLUSTER_H

#include "Tree/LengthCluster.h"


#include <map>

/**
 * a HeterogeneousCluster is a cluster
 * with 1 branch length information, an average substitution rate and
 * a set of substitution model parameters
 */
class HeterogeneousCluster : public LengthCluster {
protected:
    vector<double> parameters;
    double averageRate;
    
    map<unsigned int, unsigned int> modelCount;
    
public:    
    /** ************************************************************************
     * LengthCluster
     * @semantics  default constructor
     ************************************************************************ */
    HeterogeneousCluster( unsigned int totalNumberSpecies );
    
    /** ************************************************************************
     * LengthCluster
     * @semantics  constructor with 1 initial value
     ************************************************************************ */
    HeterogeneousCluster( unsigned int totalNumberSpecies, double initialLength );
    
    
    /** ************************************************************************
     * dup
     * @semantics  duplique the current cluster to store it in a cluster set
     ************************************************************************ */    
    virtual Cluster* dup() const{
        return new HeterogeneousCluster(*this);
    }
    
     /** ************************************************************************
     * addParam
     * @input      another cluster
     * @semantics  during a consensus we have to 'merge' the parameters of
     *             similar cluster, merge to this cluster the other one
     ************************************************************************ */    
    virtual void addParam( const Cluster& other );
    
     /** ************************************************************************
     * addModelParameters
     * @input      initial averageRate and parameters
     * @semantics  after the cluster is created we want to give it initial
     *             values from the tree that created it
     ************************************************************************ */    
    void addModelParameters( double averageSubstitutionRate,
                             const vector<double> & param,
                             unsigned int model );
    
    /** ************************************************************************
     * getAverageSubstitutionRate
     * @semantics  return the average substitution rate
     ************************************************************************ */
    inline double getAverageSubstitutionRate(){
        return averageRate;
    }

    /** ************************************************************************
     * getAverageSubstitutionRate
     * @semantics  return the average substitution rate
     ************************************************************************ */
    inline const vector< double > & getParameters(){
        return parameters;
    }
    
    /** ************************************************************************
     * getModelCount
     * @semantics  return the total count for each model
     ************************************************************************ */
    inline const map<unsigned int, unsigned int> & getModelCount(){
        return modelCount;
    }
    
};

#endif //HETEROGENEOUSCLUSTER_H
