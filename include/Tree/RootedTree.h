#ifndef ROOTEDTREE_H
#define ROOTEDTREE_H

#include "Tree/InferenceTree.h"

class ParametersSet;

class RootedTree : virtual public InferenceTree {

protected:
    /** ************************************************************************
     * RootedTree
     * @semantics  default constructor of an empty rooted tree
     *             (used by prototypes)
     ************************************************************************ */
    RootedTree();

    /** ************************************************************************
     * RootedTree
     * @semantics  common constructor of an empty rooted tree
     *             (used by descendants)
     ************************************************************************ */
    RootedTree(ParametersSet& parameters);

    /** ************************************************************************
     * ~RootedTree
     * @semantics  destructor of an rooted inference tree
     ************************************************************************ */
    virtual ~RootedTree();

    //optionnal name of a cluster or species to root the tree.
    //In effect, we just have to add the opposite cluster of rooting
    //to enforce this constraint.
    string rooting;

public:
    /** ************************************************************************
     * checkBinary/isBinary
     * @semantics  check that the tree is a strictly bifurcating binary tree
     *             after construction. checkBinary exit with an error message
     *             if this is not the case
     ************************************************************************ */
    void checkBinary(BasicNode* node = NULL);
    BasicNode* isBinary(BasicNode* node);
        
     /** ************************************************************************
     * constructFromString
     * @input          the string representation of the tree
     * @semantics      initialise the tree with a user-specified string
     ************************************************************************ */
    virtual void constructFromString( string stringTree );

    /** ************************************************************************
     * constructRandomly
     * @input      the maximum length for the branches during the creation
     * @semantics  create a random initial tree for the species in table
     ************************************************************************ */
    virtual void constructRandomly( double maxBranchLength );

    /** ************************************************************************
     * findBranches
     * @input      maxdistance, the distance from the root or from the ancestor
     *             of the given cluster
     * @input      excluded, a node to exclude during the descent (can be NULL)
     * @input      crossOnly, a boolean
     * @input      orig, a cluster to be used for the root of the search
     *             and to define its limit (can/must be NULL if no cluster
     *             to be used)
     * @return     sum, the sum of the distance in the vector
     * @return     a vector of node (equivalent to branches) with distance
     *             information
     * @semantics  Return the branches within the specified range from the
     *             (cluster) root. If equalOnly = false all the branches whose origin
     *             is at a distance lower than the specified one from
     *             the root are returned. The returned vector contains
     *             1. the destination of the branches
     *             2. the length of the branch in the specified radius
     *             (if destination is within the specified range then the
     *             branch length is returned)
     *             if crossOnly = false then the returned vector contains
     *             only the branches that crosses the distance value.
     *             this call is a recursive call
     ************************************************************************ */
    double findBranches( vector< pair<InferenceNode*,double> > & branchVector,
                       double distanceMax, InferenceNode* excluded = NULL, bool crossOnly = false, ClustersTreeNode* orig = NULL ) const;
    double findBranchesAux( vector< pair<InferenceNode*,double> > & branchVector,
                       double distanceLeft, InferenceNode* excluded, bool crossOnly,
                       InferenceNode* node ) const;


    /** ************************************************************************
     * loadDataAndModel
     * @input      sequences and model
     * @semantics  cf InferenceTree::loadDataAndModel
     *             redefined in case a root was specified (to add a cluster
     *             in the clusterTree)
     ************************************************************************ */
    virtual void loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel );
};

#endif
