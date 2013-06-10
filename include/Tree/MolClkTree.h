#ifndef MOLCLKTREE_H
#define MOLCLKTREE_H

#include "Tree/RootedTree.h"

#include <assert.h>
#include <math.h>

class ParametersSet;

class MolClkTree : public RootedTree {

private:
    double height;

protected:
    /** ************************************************************************
     * MolClkTree
     * @semantics  default constructor of an empty rooted tree
     *             (used by prototypes)
     ************************************************************************ */
    MolClkTree();

    /** ************************************************************************
     * MolClkTree
     * @semantics  common constructor of an empty rooted tree
     *             (used by descendants)
     ************************************************************************ */
    MolClkTree(ParametersSet& parameters);

    /** ************************************************************************
     * ~MolClkTree
     * @semantics  destructor of an rooted inference tree
     ************************************************************************ */
    virtual ~MolClkTree();

public:
    /** ************************************************************************
     * clockTree
     * @return     true if the subtree was equilibrated
     * @return     length, the length from startNode to the leaves
     * @semantics  method used to give equal length to a tree, return the
     *             length from a node to leaves,...
     ************************************************************************ */
     bool clockTree(double& length, InferenceNode* startNode = NULL);

    /** ************************************************************************
     * rescale
     * @input         scaling factor
     * @input         first node of the subtree to rescale (the length to its
     *                parent is rescaled too)
     * @post          invalidate recursively FROM startNode after the call
     * @semantics     method used to rescale the branch lengths of a subtree
     ************************************************************************ */
    void rescale(double scalingFactor, InferenceNode* startNode);

    /** ************************************************************************
     * rescale
     * @input         scaling factor
     * @semantics     method used to change the total length of the tree
     *                The new length parameter is length*scalingFactor
     *                and branches are multiplied by scaz\lingFactor to stay
     *                consistent
     ************************************************************************ */
    void rescale(double scalingFactor);

     /** ************************************************************************
     * constructFromString
     * @input          the string representation of the tree
     * @semantics      initialise the tree with a user-specified string
     * @precondition   one leaf label in stringTree is the outgroup
     ************************************************************************ */
    virtual void constructFromString( string stringTree );

    /** ************************************************************************
     * constructRandomly
     * @input      the sequence table of the species to include in the new tree
     * @input      the maximum length for the branches during the creation
     * @semantics  create a random initial tree for the species in table
     ************************************************************************ */
    virtual void constructRandomly( double maxBranchLength );

    /** ************************************************************************
     * setBranchVector
     * @input              new branch vector
     * @preconditions      branches compatible with a rooted tree
     ************************************************************************ */
    virtual void setBranchVector( const vector<double> & branchesLength );

    /** ************************************************************************
     * getHeight
     * @return    the distance from the root to each leaves
     ************************************************************************ */
    inline double getHeight() const {
        assert(!isnan(height));
        assert(!isinf(height));
        assert(height);
        return height;
    }

};

#endif
