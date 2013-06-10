#ifndef SEARCHTREESTARSEARCH_H
#define SEARCHTREESTARSEARCH_H

#include "Tree/SearchTreeBasic.h"

#include <deque>

/** ***************************************************************************
 * Heuristic ML search, species are added sequentially in a random order
 * All the possible branches are tried when adding a new one and the
 * tree with the best likelihood is kept for the next step.
 * Algorithms are complicated by the presence of clusters
 **************************************************************************** */
class SearchTreeStarSearch : virtual public OptimizerTree, public SearchTreeBasic{
protected:
    /** ************************************************************************
     * SearchTreeStarSearch
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< SearchTree >
     ************************************************************************ */
    SearchTreeStarSearch( const string & registrationName );

    /** ************************************************************************
     * SearchTreeStarSearch
     * @semantics  constructor of an empty SearchTreeStarSearch according to the
     *             parameters, this constructor should be called by all
     *             descendant
     ************************************************************************ */
    SearchTreeStarSearch( ParametersSet& treeParameters );

    
    /** ************************************************************************
     * recursiveBuilt
     * @input      ptable, the sequences to know tips label
     * @input      clustersTreeNode, the node from where to build
     * @input      an initialized node to be used as an ancestor for the
     *             cluster
     * @semantics  initial tree building method (recursive) when clusters are
     *             used
     ************************************************************************ */
    void recursiveBuilt(SequenceTable * ptable, ClustersTreeNode* clustersTreeNode, OptimizerNode* clusterAncestor);
    
    /** ************************************************************************
     * link/unlink
     * @input          node1, node2 the two nodes to link/unlink
     * @semantics      add/remove a link
     * @preconditions  some...
     ************************************************************************ */
    void link( OptimizerNode* node1, OptimizerNode* node2, OptimizerNode* newParent );
    void unlink( OptimizerNode* node1, OptimizerNode* node2, OptimizerNode* newParent );
            
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
     * ~SearchTreeStarSearch
     * @semantics  destructor of a the step wise addition heuristic ML search
     ************************************************************************ */
    ~SearchTreeStarSearch();

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

    static SearchTreeStarSearch prototype;
    //store the nodes with unresolved bifurcations
    list< OptimizerNode* > centers;
};

#endif //SEARCHTREESTEPWISEADD_H
