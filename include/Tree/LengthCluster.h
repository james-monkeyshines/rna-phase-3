#ifndef LENGTHCLUSTER_H
#define LENGTHCLUSTER_H

#include "Tree/Cluster.h"

/**
 * a LengthCluster is a cluster
 * with 1 branch length information 
 */
class LengthCluster : public Cluster {
protected:
    double length;
    unsigned int number;

public:    
    /** ************************************************************************
     * LengthCluster
     * @semantics  default constructor
     ************************************************************************ */
    LengthCluster( unsigned int totalNumberSpecies );
    
    /** ************************************************************************
     * LengthCluster
     * @semantics  constructor with 1 initial value
     ************************************************************************ */
    LengthCluster( unsigned int totalNumberSpecies, double initialLength );
    
    /** ************************************************************************
     * getLength
     * @semantics  return the length
     ************************************************************************ */
    inline double getLength(){
        return length;
    }
    
    /** ************************************************************************
     * getNumber
     * @semantics  return the number of calls to addParam 
     ************************************************************************ */
    inline unsigned int getNumber(){
        return number;
    }

    
    /** ************************************************************************
     * dup
     * @semantics  duplique the current cluster to store it in a cluster set
     ************************************************************************ */    
    virtual Cluster* dup() const{
        return new LengthCluster(*this);
    }
    
     /** ************************************************************************
     * addParam
     * @input      another cluster
     * @semantics  during a consensus we have to 'merge' the parameters of
     *             similar cluster, merge to this cluster the other one
     ************************************************************************ */    
    virtual void addParam( const Cluster& other );
};

#endif //LENGTHCLUSTER_H
