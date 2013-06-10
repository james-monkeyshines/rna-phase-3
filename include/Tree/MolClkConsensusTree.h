#ifndef MOLCLKCONSENSUSTREE_H
#define MOLCLKCONSENSUSTREE_H

#include "Tree/ConsensusTree.h"

class MolClkConsensusTree : public ConsensusTree{

public:

    /** ************************************************************************
     * MolClkConsensusTree
     * @input      parameters, a ParametersSet containing information about how
     *             the sampling was done
     * @semantics  initialise the consensus tree with the same parameters used
     *             to initialise the sampling
     ************************************************************************ */
    MolClkConsensusTree( ParametersSet& parameters );
    
    /** ************************************************************************
     * clone
     * @input      a ParametersSet used to create a new ConsensusTree
     ************************************************************************ */
    virtual ConsensusTree* clone( ParametersSet& parameters ) const;
    
    /** ************************************************************************
     * ~ConsensusTree
     * @semantics  virtual destructor
     ************************************************************************ */
    virtual ~MolClkConsensusTree(){};
    
protected:    
    /** ************************************************************************
     * MolClkConsensusTree
     * @input      a registration name to save prototype in the ConsensusTree
     *             factory
     * @semantics  save prototypes of this class and of its descendant in the
     *             factory
     ************************************************************************ */
    MolClkConsensusTree( const string & registrationName );

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
    static MolClkConsensusTree prototype;
};

#endif //MOLCLKCONSENSUSTREE_H
