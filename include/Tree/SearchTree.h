#ifndef SEARCHTREE_H
#define SEARCHTREE_H

#include "Tree/OptimizerTree.h"

class SearchTree : virtual public InferenceTree, virtual public OptimizerTree{

public:
    /** ************************************************************************
     * ~MLSearchTree
     * @semantics  destructor of a ML search inference tree
     ************************************************************************ */
    virtual ~SearchTree(){};

    /** ************************************************************************
     * processNext
     * @semantics  go to the next step of the search
     * @return     bool false if the search is over.
     ************************************************************************ */
    virtual bool processNext() = 0;

    /** ************************************************************************
     * optimizeFlag
     * @semantics  set the model optimization flag
     ************************************************************************ */
    virtual void optimizeFlag( bool flag ) = 0;

    /** ************************************************************************
     * empiricalFreqsFlag
     * @semantics  set the empirical frequencies flag
     ************************************************************************ */
    virtual void empiricalFreqsFlag( bool flag ) = 0;

    /** ************************************************************************
     * printResults
     * @semantics  print the results
     ************************************************************************ */
    virtual void printResults( ostream& outputStream ) = 0;

    /** ************************************************************************
     * printTree
     * @semantics  print only the tree
     ************************************************************************ */
    virtual void printTree( ostream& outputStream ) = 0;
};

#endif //SEARCHTREE_H
