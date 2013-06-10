#ifndef INFERENCETREE_H
#define INFERENCETREE_H

#include <string>

#include "Util/statlib.h"
#include "Util/array3D.h"

#include "Tree/Tree.h"
#include "Tree/InferenceNode.h"

class BasicNode;
class array_bag;
class randombox;

class SequenceTable;
class Model;


#define NUMBER_LOOKUP_TABLE 3
//if a node has more than 3 children expect some trouble
//with lookup. In most cases a node can have at most three children
//(and 2 would be enough if we were sure it would never become root).
//ML trees should handle the problem of multifurcating nodes by
//themselves


class InferenceTree : public Tree {
protected:

    /** ************************************************************************
     * pmodel
     * a pointer to the substitution model used
     ************************************************************************ */
    Model * pmodel;

    /** ************************************************************************
     * ptable
     * a pointer to the sequence table
     ************************************************************************ */
    SequenceTable * ptable;

    /** ************************************************************************
     * clustersStream
     * a pointer to the file stream to read clusters
     ************************************************************************ */
    ifstream* clustersStream;

    /** ************************************************************************
     * clustersTree
     * a pointer to the clusters
     ************************************************************************ */
    ClustersTree* clustersTree;

    /** ************************************************************************
     * lookupResize
     * @semantics  For ML optimizer Tree there is no guarantee that the number
     *             of children will remain below 3.
     *             this function resize lookup table to appropriate values
     *             once the tree is built
     ************************************************************************ */
    void lookupResize();

public:

    inline Model* getModel() const{
        return pmodel;
    }
    inline SequenceTable* getTable() const{
        return ptable;
    }


    /** ************************************************************************
     * InferenceTree
     * @semantics  default constructor of an empty inference tree
     ************************************************************************ */
    InferenceTree();

    /** ************************************************************************
     * InferenceTree
     * @input      the string representation of the tree
     * @semantics  constructor from a string
     ************************************************************************ */
    InferenceTree( string stree );

    /** ************************************************************************
     * InferenceTree
     * @input      another tree
     * @semantics  copy constructor
     ************************************************************************ */
    InferenceTree( InferenceTree & );

    /** ************************************************************************
     * ~InferenceTree
     * @semantics  destructor of an inference tree
     ************************************************************************ */
    virtual ~InferenceTree();

    /** ************************************************************************
     * invalidateNodes
     * @semantics  switch all the nodes to the working state for the given
     *             category (-1 = all)
     ************************************************************************ */
    virtual void invalidateNodes( int cat );

    /** ************************************************************************
     * retrieveNodes
     * @semantics  switch all the nodes back to their saved state for the given
     *             category (-1 = all)
     ************************************************************************ */
    virtual void retrieveNodes( int cat );

    /** ************************************************************************
     * saveNodesAux
     * @semantics  the current partialLikehood configuration becomes the new
     *             saved state (-1 = all)
     ************************************************************************ */
    virtual void saveNodes( int cat );

    /** ************************************************************************
     * operator =
     * @input      another tree
     * @semantics  the src tree is copied over the tree, optimisation are used
     *             if the two trees have the same model and sequence table
     ************************************************************************ */
    InferenceTree & operator = ( const InferenceTree & src );

    /** ************************************************************************
     * operator =
     * @input      another tree using same model and sequence table
     * @input      a boolean to decide whether the computation must be copied
     *             too
     * @semantics  the src tree is copied over the tree, optimisation are used
     *             if the two trees have the same model and sequence table
     ************************************************************************ */
    InferenceTree & quickCopy( const InferenceTree & src,
                               bool withPartialLikelihood );

    /** ************************************************************************
     * toString
     * @input      outputLengths, set to true if the representation must include
     *             lengths information
     * @semantics  return the string representation of this tree. The tree
     *             should be sorted before.
     ************************************************************************ */
    string toStringNumbered( bool topologyOnly = true );

    /** ************************************************************************
     * move
     * @input         the BasicNode to move, the new father
     * @semantics     move a node to another part of the tree
     *                the distance to the parent remains the same
     * @precondition  no cycle allowed, src and dest in the tree, ...
     ************************************************************************ */
    void move( InferenceNode* src, InferenceNode* dest );


    /** ************************************************************************
     * swap
     * @input         the two nodes to swap
     * @semantics     swap 2 nodes, nodes take the distance to their parent
     *                with them
     * @precondition  no cycle allowed, src and dest in the tree, ...
     ************************************************************************ */
    void swap( InferenceNode* node1, InferenceNode* node2 );


    /** ************************************************************************
     * remove
     * @input         the InferenceNode to remove
     * @semantics     remove a node from the tree
     ************************************************************************ */
    void remove( InferenceNode* src );


    /** ************************************************************************
     * add
     * @input         the InferenceNode to add, the new father
     * @semantics     add a node to the tree
     * @precondition  no cycle allowed, dest in the tree, ...
     * @precondition  TEMPORARILY src MUST HAVE BEEN INITIALISED PROPERLY
     *                (ie, must have been a node of the tree when
     *                initialisation was called)
     ************************************************************************ */
    void add( InferenceNode* src, InferenceNode* dest, double distance );


    /** ************************************************************************
     * setBranchLength
     * @input         the ID of the branch
     * @input         the new length
     * @semantics     change the length of the specified branch,
     *                invalidate the likelihood computation
     * @precondition  0 <= branchId < numberBranches
     ************************************************************************ */
    void setBranchLength( unsigned int branchId, double newLength );

    /** ************************************************************************
     * makeRoot
     * @input         the new root for the tree
     * @semantics     change the root of the tree
     * @precondition  newRoot is an internal node
     ************************************************************************ */
    bool makeRoot( InferenceNode* newRoot );

    /** ************************************************************************
     * loglikelihood
     * @semantics      return the loglikelihood of the tree
     * @preconditions  model and sequence loaded
     ************************************************************************ */
    double loglikelihood();


    /** ************************************************************************
     * loadDataAndModel
     * @input           the sequences and the model to use with this tree
     * @preconditions   no data or model already loaded
     ************************************************************************ */
    virtual void loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel );

    /** ***********************************************************************
     * a memory space allocated to avoid the allocation at each call to
     * likelihood for the nodes
     ************************************************************************ */
    vector< array3D< double >* > lookup;


    /** ************************************************************************
     * constructRandomly
     * @input      the sequence table of the species to include in the new tree
     * @input      the maximum length for the branches during the creation
     * @semantics  create a random initial tree for the species in table
     ************************************************************************ */
    virtual void constructRandomly( double maxBranchLength );


    /** ************************************************************************
     * constructFromString
     * @input      the string representation of the tree
     * @semantics  initialise the tree with a user-specified string
     ************************************************************************ */
    virtual void constructFromString( string stringTree );

#ifdef DEBUG1
    /** ************************************************************************
     * checkIdentical
     * @semantics  Function provided for debugging purpose, check the identical
     *             array for all the nodes
     ************************************************************************ */
    void checkIdentical();
#endif
    

    /** ************************************************************************
     * loadClusters
     * @input           filename
     * @semantics       Declare a cluster file
     *                  The reading of the file is delayed until the sequence
     *                  has been loaded.
     ************************************************************************ */
    void loadClusters( const string& fileName );

    /** ************************************************************************
     * processClusters
     * @semantics       load clusters from their fileStream, check consistency
     *                  with the tree
     ************************************************************************ */
    void processClusters();

#ifdef DEBUG1
    /** ************************************************************************
     * checkClusters
     * @semantics       check whether the tree is conform to the clusters
     *                  (debug method)
     ************************************************************************ */
    void checkClusters();
#endif
    
    /** ************************************************************************
     * assignClusters
     * @semantics       check whether the user's tree is conform to the
     *                  clusters. Method called to check user's input.
     *                  ancestor nodes of clusters are initialised with the
     *                  right node in the tree
     ************************************************************************ */
    void assignClusters();


    /** ************************************************************************
     * findClusterBranches
     * @input      a cluster
     * @return     sum, the sum of the distance in the vector
     * @return     a vector of node (equivalent to branches) with distance
     *             information
     * @semantics  Return all the branches (children nodes) which belong to the
     *             cluster (with distance information)
     ************************************************************************ */
    double findClusterBranches( vector< pair<InferenceNode*,double> > & branchVector,
                                ClustersTreeNode* clustersTreeNode ) const;

    /** ************************************************************************
     * findClusterAncestor
     * @input      a node
     * @return     the MRCA of the cluster this nodes belongs too
     ************************************************************************ */
    InferenceNode* findClusterAncestor(InferenceNode* node) const;

    /** ************************************************************************
     * findInsertion
     * @input      a node
     * @return    the ClustersTreeNode this node belongs to.
     ************************************************************************ */
    ClustersTreeNode* findInsertion( InferenceNode* node ) const;
    ClustersTreeNode* findInsertionRec( InferenceNode* node, ClustersTreeNode* ancestorCluster ) const;

    /** ************************************************************************
     * getClusterDistance
     * @input     a node
     * @return    orig, if clusters were defined, orig is the clustersTreeNode
     *            which contains all the descendants of node (STRICT inclusion),
     *            ie if node is already the ancestor of a cluster then orig
     *            will be the next bigger cluster.
     * @return    the distance to the ancestor node of the cluster this node
     *            belongs to. The distance to the root is returned if no
     *            cluster are defined.
     ************************************************************************ */
    double getClusterDistance( InferenceNode* node, ClustersTreeNode*& orig ) const;

    /** ************************************************************************
     * getRootDistance
     * @input     a node
     * @return    distance to the root
     ************************************************************************ */
    double getRootDistance( BasicNode* node ) const;

};

#endif //INFERENCETREE_H




