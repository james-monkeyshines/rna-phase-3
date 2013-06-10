#ifndef OPTIMIZERNODE_H
#define OPTIMIZERNODE_H

#include "Tree/InferenceNode.h"

class OptimizerTree;


class OptimizerNode : public InferenceNode{

public:
    /** ************************************************************************
     * OptimizerNode
     * @semantic  the tree the node will be attached to
     * @semantic  constructor of a detached node
     ************************************************************************ */
    OptimizerNode( OptimizerTree* ptree );

    /** ************************************************************************
     * OptimizerNode
     * @input     leafLabel, the label we attach the node to
     * @semantic  constructor of a leaf node
     ************************************************************************ */
    OptimizerNode( const string& label, OptimizerTree* ptree );
    
    /** ************************************************************************
     * OptimizerNode
     * @input     the treemap and the tree
     * @semantic  constructor of an arborescence of node
     ************************************************************************ */
    OptimizerNode( const TreeMap& treeMap, OptimizerTree* ptree );

    /** ************************************************************************
     * OptimizerNode
     * @semantic  destructor of node
     ************************************************************************ */
    virtual ~OptimizerNode();

    /** ************************************************************************
     * initialisationOpti
     * @input         leafNode, a boolean to indicate if the node is a leaf
     * @semantic      initialisation of the two previous arrays once the tree
     *                has got a sequence table and a model
     * @precondition  the model is known, pnodeModel != NULL,
     *                (call InferenceNode::initialisation( bool leafNode ) first)
     ************************************************************************ */
    void initialisationOpti( bool leafNode );

 
    /** ************************************************************************
     * computeBeliefPrior
     * @semantic       fill the belief and prior arrays the call starts from
     *                 the root go down to the leaves computing the belief
     *                 arrays for all internal nodes then return to the root
     *                 filling the prior arrays
     * @preconditions  partialLikelihood(s) were properly initialised thanks to
     *                 a call to likelihood(), do not forget to invalidate the
     *                 whole tree before in case some nodes are in a trusted
     *                 saved state which shouldn't be trusted.
     ************************************************************************ */
    void computeBeliefPrior();
       
    /** ***********************************************************************
     * belief and prior storage for each sequences
     * (ie. each model in this case)
     * dim vector = number of models
     * dim array : number of symbol * number of rate * number of site
     ************************************************************************ */
    vector< array3D< double > > belief;
    vector< array3D< double > > prior;
            
  
};

#endif //OPTIMIZERNODE_H
