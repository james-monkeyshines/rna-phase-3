#ifndef SEARCHTREEINSERTIONBASIC_H
#define SEARCHTREEINSERTIONBASIC_H

#include "Tree/SearchTreeBasic.h"

#include <deque>
#include <stack>

class OptimizerNode;


/*
 * we define in this class, useful primitives for the search heuristic
 * that works using partial trees (Step-wise addition, Branch-and-Bound,
 * Exhaustive search).
 * this class helps to handle the issues that arise when a tree is
 * constructed step-by-step (ie problems with clusters and rooting point)
 */

class SearchTreeInsertionBasic : virtual public OptimizerTree, public SearchTreeBasic{
protected:
    /** ************************************************************************
     * SearchTreeInsertionBasic
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< SearchTree >
     ************************************************************************ */
    SearchTreeInsertionBasic( const string & registrationName );

    /** ************************************************************************
     * SearchTreeInsertionBasic
     * @semantics  constructor of an empty SearchTreeStepWiseAdd according to the
     *             parameters, this constructor should be called by all
     *             descendant
     ************************************************************************ */
    SearchTreeInsertionBasic( ParametersSet& treeParameters );

    /** ************************************************************************
     * insertion/removal
     * @input          the node to insert
     * @input          its futur parent
     * @input          the branch to perform the insertion
     * @preconditions  node is a leaf
     * @preconditions  insertPoint is a valid branch to isert the node
     * @semantics      split the branch designed by insertPoint, put newParent
     *                 in the middle and add node to parent
     *
     ************************************************************************ */
    void insertion( OptimizerNode * node, OptimizerNode * newParent, OptimizerNode * insertNode );
    void removal( OptimizerNode * node, OptimizerNode * newParent, OptimizerNode * insertNode );

    /** ************************************************************************
     * insertionValidate/removalValidate
     * @input          the node to insert
     * @input          its futur parent
     * @input          the branch to perform the insertion
     * @preconditions  node is a leaf
     * @preconditions  insertion/removal was called just before with the same parameters
     * @semantics      split the branch designed by insertPoint, put newParent
     *                 in the middle and add node to parent
     ************************************************************************ */
    void insertionValidate( OptimizerNode * node, OptimizerNode * newParent, OptimizerNode * insertNode );
    void removalValidate( OptimizerNode * node, OptimizerNode * newParent, OptimizerNode * insertNode );

public:
    /** ************************************************************************
     * initialisation
     * @input          the sequences and the model to use with this tree
     * @preconditions  no data or model already loaded
     * @semantics      call SearchTreeBasic::initialise
     *                 prepare the first tree
     ************************************************************************ */
    virtual bool initialisation( SequenceTable * ptable, Model * pmodel );


    /** ************************************************************************
     * createInitialThreeLeavesTree
     * @semantics  use the three species provided to create an initial tree
     *             3-species unrooted tree
     ************************************************************************ */
    void createInitialThreeLeavesTree( unsigned int spec1, unsigned int spec2, unsigned int spec3 );

    /** ************************************************************************
     * retrievePossibleBranches
     * @input          the index of the species we are trying to add
     * @output         a vector of branches
     * @semantics      return the possible branches to add this new species
     *                 without breaking the clusters constraints.
     *                 NB:if the insertion is done, some clusters might
     *                 have to be moved around and the tree rerooted.
     ************************************************************************ */
    void retrievePossibleBranches( vector< InferenceNode* >& possibleBranches,
                                   unsigned int speciesId );



protected:

    //clusters who are not yet clearly established into the tree
    //will be parked at their first sibling
    //waitingClusters tells at which node the cluster is waiting
    //clustersQueue gives a queue of waiting clusters at the node
    //(convenient redundancy?)
    map< ClustersTreeNode*, OptimizerNode* > waitingClusters;
    map< OptimizerNode*, deque<ClustersTreeNode*> > clustersQueue;

    //the root has three children and some clusters will apply only to
    //two of them
    deque<ClustersTreeNode*> rootClusters;
    //tempOutgroup is the child of root which is not in rootClusters
    BasicNode* tempOutgroup;

    stack< pair < OptimizerNode*, OptimizerNode* > > rootSave;

private:

    /** ************************************************************************
     * reroot/cancelReroot
     * @semantics  update rootClusters and clusters waiting at the oldRoot
     *             when a rerooting is performed. outgroup is used
     *             as the new tempOutgroup if insertionPoint was not
     *             in rootClusters.
     ************************************************************************ */
    void reroot( OptimizerNode* oldRoot, OptimizerNode* newRoot, ClustersTreeNode* insertionPoint, OptimizerNode* outgroup );
    void cancelReroot( OptimizerNode* cancelledNewRoot, OptimizerNode* restoredRoot, ClustersTreeNode* oldInsertionPoint, OptimizerNode* outgroup );
};

#endif
