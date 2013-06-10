#ifndef HETEROGENEOUSMOLCLKCONSENSUSTREE_H
#define HETEROGENEOUSMOLCLKCONSENSUSTREE_H

#include "Tree/HeterogeneousConsensusTree.h"


class HeterogeneousMolClkConsensusTree : public HeterogeneousConsensusTree{

public:
    /** ************************************************************************
     * MolClkConsensusTree
     * @input      parameters, a ParametersSet containing information about how
     *             the sampling was done
     * @semantics  initialise the consensus tree with the same parameters used
     *             to initialise the sampling
     ************************************************************************ */
    HeterogeneousMolClkConsensusTree( ParametersSet& parameters );
    
    /** ************************************************************************
     * clone
     * @input      a ParametersSet used to create a new ConsensusTree
     ************************************************************************ */
    virtual ConsensusTree* clone( ParametersSet& parameters ) const;
    
    /** ************************************************************************
     * ~ConsensusTree
     * @semantics  virtual destructor
     ************************************************************************ */
    virtual ~HeterogeneousMolClkConsensusTree(){};
    
protected:    
    /** ************************************************************************
     * HeterogeneousMolClkConsensusTree
     * @input      a registration name to save prototype in the ConsensusTree
     *             factory
     * @semantics  save prototypes of this class and of its descendant in the
     *             factory
     ************************************************************************ */
    HeterogeneousMolClkConsensusTree( const string & registrationName );

    /** ************************************************************************
     * createCluster
     * @input      a node of a sampled tree
     * @return     create an empty cluster for that node
     * @semantics  create the type of cluster specific for the consensus tree
     *             and fill its parameters part.
     ************************************************************************ */
    virtual Cluster* createCluster(BasicNode* node);
    
    /** ************************************************************************
     * getDist
     * @input      two clusters
     * @return     the distance to put between father and child
     * @semantics  each kind of cluster contains information to set the length
     *             in the final consensus tree. getDist is virtual too allow
     *             modification by descendant class
     ************************************************************************ */
    virtual double getDist( Cluster* father, Cluster* child );
    
private:
    static HeterogeneousMolClkConsensusTree prototype;    
};

#endif //HETEROGENEOUSMOLCLKCONSENSUSTREE_H
