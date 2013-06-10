#ifndef SEARCHTREESTEPWISEADD_H
#define SEARCHTREESTEPWISEADD_H

#include "Tree/SearchTreeInsertionBasic.h"

#include <deque>

/** ***************************************************************************
 * Heuristic ML search, species are added sequentially in a random order
 * All the possible branches are tried when adding a new one and the
 * tree with the best likelihood is kept for the next step.
 * Algorithms are complicated by the presence of clusters
 **************************************************************************** */
class SearchTreeStepwiseAdd : virtual public OptimizerTree, public SearchTreeInsertionBasic{
protected:
    /** ************************************************************************
     * SearchTreeStepwiseAdd
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< SearchTree >
     ************************************************************************ */
    SearchTreeStepwiseAdd( const string & registrationName );

    /** ************************************************************************
     * SearchTreeStepwiseAdd
     * @semantics  constructor of an empty SearchTreeStepWiseAdd according to the
     *             parameters, this constructor should be called by all
     *             descendant
     ************************************************************************ */
    SearchTreeStepwiseAdd( ParametersSet& treeParameters );

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
     * ~SearchTreeStepwiseAdd
     * @semantics  destructor of a the step wise addition heuristic ML search
     ************************************************************************ */
    ~SearchTreeStepwiseAdd();

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

    static SearchTreeStepwiseAdd prototype;


    vector<unsigned int> additionOrder;
    unsigned int index;
};

#endif //SEARCHTREESTEPWISEADD_H
