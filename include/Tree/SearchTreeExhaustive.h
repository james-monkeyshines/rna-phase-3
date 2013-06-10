#ifndef SEARCHTREEEXHAUSTIVE_H
#define SEARCHTREEEXHAUSTIVE_H

#include "Tree/SearchTreeInsertionBasic.h"

#include <deque>
#include <stack>


using namespace std;

/** ***************************************************************************
 * Exhaustive ML search, species are added sequentially in a random order.
 * All the possible topologies are tried. (but not the subtopologies)
 * Algorithms are complicated by the presence of clusters. To make it
 * easier, algorithms were adapted from the branch-bound search. Indeed, the
 * exhaustive search is a kind of special case of the branch-bound search
 * without optimisations and check of the upperbound for subtopologies
 * it makes the code look a bit more complicated than it should be but
 * it avoids trouble with clusters and rerooting. The computational cost
 * associated with this code overhead is not an issue.
 **************************************************************************** */
class SearchTreeExhaustive: virtual public OptimizerTree, public SearchTreeInsertionBasic{
protected:
    /** ************************************************************************
     * SearchTreeExhaustive
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< SearchTree >
     ************************************************************************ */
    SearchTreeExhaustive( const string & registrationName );

    /** ************************************************************************
     * SearchTreeExhaustive
     * @semantics  constructor of an empty SearchTreeExhaustive according to the
     *             parameters, this constructor should be called by all
     *             descendant
     ************************************************************************ */
    SearchTreeExhaustive( ParametersSet& treeParameters );


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
     * ~SearchTreeExhaustive
     * @semantics  destructor of a the step wise addition heuristic ML search
     ************************************************************************ */
    ~SearchTreeExhaustive();

    /** ************************************************************************
     * clone
     * @semantics  the clone method for the factory<SearchTree>
     ************************************************************************ */
    virtual OptimizerTree* clone(ParametersSet& parameters) const;

    /** ************************************************************************
     * processNext
     * @semantics  go to the next step of the search
     * @return     bool true if the search is over.
     ************************************************************************ */
    virtual bool processNext();

    /** ************************************************************************
     * printResults
     * @semantics  print the results
     ************************************************************************ */
    virtual void printResults( ostream& outputStream );

    /** ************************************************************************
     * printTree
     * @semantics  print only the tree
     ************************************************************************ */
    virtual void printTree( ostream& outputStream );

private:
    /** ************************************************************************
     * backTrack
     * return      true if the search is finished
     * @semantics  method called by processNext and buildRec when a backtrack
     *             is needed during the search. It is the case when the
     *             likelihood (+penalty) of a subtree constructed in buildRec
     *             is lower than the upperBound or when a full tree was
     *             constructed and found to be the best one.
     *             if return==false then the next possible tree has been built
     ************************************************************************ */
    bool backTrack();

    /** ************************************************************************
     * buildRec
     * return      true if the search is finished
     *             false if a new tree is found
     * @semantics  recursive method:
     *             1.the current tree is optimised and we check the bound
     *             2.if the bound was crossed we backtrack and we continue the
     *               search by another call to buildRec
     *             3.otherwise we prepare the next layer to add the next
     *               species (and we call buildRec again)
     *             4.if a better candidate tree is found (ie, better upperbound)
     *               we save the new best state and we backtrack
     *             5.true is returned when the backtrack is one of the above
     *               case is not possible anymore
     ************************************************************************ */
    bool buildRec();


    static SearchTreeExhaustive prototype;


    vector<unsigned int> additionOrder;

    //since it is a backtracking algorithm we use a stack...
    //at each step we have a set of possible branches to try
    stack< vector< OptimizerNode* > > branchesStack;
    //this iterator records what should be the next try
    stack< vector< OptimizerNode* >::const_iterator > branchIter;


    //to save the best state found during the search
    string bestTree;
    double bestLnLik, bestLnPenalty;
    vector< double > bestModelParameters;
    vector< double > bestModelPenalty;

    vector< double > initialModelParameters;
    vector< double > initialModelPenalty;

    //nodes to insert are stored here (created at initialisation)
    vector< pair< OptimizerNode*, OptimizerNode* > > insertNodes;

};

#endif //SEARCHTREEBRANCHBOUND_H
