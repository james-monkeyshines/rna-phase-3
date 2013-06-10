#ifndef INFERENCENODE_H
#define INFERENCENODE_H

#include "Tree/BasicNode.h"

#include "Sequence/SequenceTable.h"

#include "Util/array3D.h"

#include <vector>

class InferenceTree;
class Model;

class Cluster;
class ClustersTree;
class ClustersTreeNode;

using namespace std;

class InferenceNode : public BasicNode {

protected:
    /** ************************************************************************
     * pnodeModel
     * a pointer to the substitution model used by this node
     ************************************************************************ */
    Model* pnodeModel;

    /** ************************************************************************
     * InferenceNode
     * @semantic  default constructor
     ************************************************************************ */
    InferenceNode();

public:

    /** ************************************************************************
     * lookup
     * a pointer to some reserved memory space in the inference tree
     ************************************************************************ */
    vector<array3D<double> *> lookup;

    /** ************************************************************************
     * lookupResize
     * @return     the biggest number of children in the subtree
     * @semantics  For ML optimizer Tree there is no guarantee that the number
     *             of children will remain below 3.
     *             this function resize lookup table to appropriate values
     *             once the tree is built, the call propagates recursively
     *             to the leaves (tip are excluded since no lookup table)
     ************************************************************************ */
    unsigned int lookupResize();
    
    /** ************************************************************************
     * getModel
     * @return   the substitution model used by this node
     ************************************************************************ */
    inline Model* getModel() const {
       return pnodeModel;
    }
    /** ************************************************************************
     * setModel
     * @semantics   change the substitution model used by this node
     ************************************************************************ */
    inline void setModel(Model* newModel) {
       pnodeModel = newModel;
    }


    /** ***********************************************************************
     * partial likelihood storage for each sequences
     * (ie. each model in this case)
     * dim vector = number of models
     * dim array : number of site * number of symbol * number of rate
     * the working space is where new values are computed, the save space
     * is for last validated value.
     * when a call to invalidate is performed, partialLikelihood will direct
     * to the work space and it means the computation has to be done every
     * time, when a call to validate is performed, save and work are swapped,
     * therefore partialLikelihood point to Save and the computation is
     * trusted
     ************************************************************************ */
    vector< array3D< double >* > partialLikelihood;
    vector< array3D< double >* > partialLikelihoodSave;
    vector< array3D< double >* > partialLikelihoodWork;

    /** ***********************************************************************
     * to speed up the likelihood computation :
     * since for a given ancestor lots of sites have not changed among its
     * children,
     * for each site, we store the common state in identical
     * (and we store -1 if there is a variation)
     * size(identical[i]) == t_model->getSequencesLength()
     * invariantSites are identical too but they are treated separatly
     * There is a vector for each model (size(lastCalc)==numberModelsUsed)
     ************************************************************************ */
    vector< vector<int> > identical;

    /** ***********************************************************************
     * lastCalc
     * in lastCalc we store the site where the computation has been done for
     * the first and the last time, with a given identical state for all
     * descendants size(lastCalc[i]) == t_model->getNumberStates()
     * There is a vector for each model (size(lastCalc)==numberModelsUsed)
     ************************************************************************ */
    vector< vector<int> > lastCalc;

    /** ***********************************************************************
     * ptree
     * a pointer to the tree
     ************************************************************************ */
     InferenceTree* ptree;


    /** ************************************************************************
     * InferenceNode
     * @semantic  the tree the node will be attached to
     * @semantic  constructor of a detached node
     ************************************************************************ */
    InferenceNode( InferenceTree* ptree );

    /** ************************************************************************
     * InferenceNode
     * @input     leafLabel, the label we attach the node to
     * @semantic  constructor of a leaf node
     ************************************************************************ */
    InferenceNode( const string& label, InferenceTree* ptree );

    /** ************************************************************************
     * InferenceNode
     * @input     the treemap and the tree
     * @semantic  constructor of an arborescence of node
     ************************************************************************ */
    InferenceNode( const TreeMap& treeMap, InferenceTree* ptree );

    /** ************************************************************************
     * InferenceNode
     * @semantic  destructor of a leaf node
     ************************************************************************ */
    virtual ~InferenceNode();

private:
    /** ************************************************************************
     * initInferenceNode
     * @semantics   constructor primitive
     ************************************************************************ */
  //  void initInferenceNode( InferenceTree* ptree, bool leafNode );

public:
    /** ************************************************************************
     * initialisation
     * @input     leafNode, a boolean to indicate if the node is a leaf
     * @semantic  initialisation of a node once the tree has got a sequence
     *            table and a model
     ************************************************************************ */
    void initialisation( bool leafNode );

    /** ************************************************************************
     * quickCopy
     * @input         the source node
     * @input         withPartialLikelihood, to copy the partialLikelihood array
     *                too
     * @semantic      copy of a node, the children are not taken into account
     *                the tree neither. See recursiveCopy for a copy of the
     *                tree structure. This function is USED ONLY FOR EFFICIENT
     *                COPY of an inference tree (InferenceTree::operator=),
     *                and during a recursiveCopy. Do not use carelessly.
     *                if src.calculated and withPartialLikelihood the partial
     *                likelihood array is copied (its quicker than doing the
     *                computation again but it is still slow)
     * @precondition  model and sequenceTable are the same for the two nodes,
     *                the node is an internal node
     ************************************************************************ */
    InferenceNode& quickCopy( const InferenceNode& src,
                               bool withPartialLikelihood);

    /** ************************************************************************
     * getChild
     * @input      the id of the children (warning, this can change)
     * @semantics  return the child
     ************************************************************************ */
    InferenceNode* getChild(unsigned int childId);

    /** ************************************************************************
     * toStringNumbered
     * @input      outputLengths, set to true if the representation must include
     *             lengths information
     * @semantics  This method outputs the string representation of the tree
     *             from this node, sequence number are used instead of labels
     ************************************************************************ */
    string toStringNumbered( bool topologyOnly ) const;

    /** ************************************************************************
     * likelihood
     * @semantic  fill the partialLikelihood array according to its children
     *            and a substitution model
     ************************************************************************ */
    void likelihood();

    /** ************************************************************************
     * computeIdentical
     * @preconditions  the tree that contains the node has a sequence table
     *                 loaded and THE CHILDREN IDENTICAL IS UP TO DATE
     *                 the identical of the parent should be modified too
     * @semantics      fill the identical array of the node, this array is
     *                 used to speed up the likelihood computation
     ************************************************************************ */
    void computeIdentical();

    /** ************************************************************************
     * updateIdenticalRecursively
     * @preconditions  the tree that contains the node has a sequence table
     *                 loaded and THE CHILDREN IDENTICAL IS UP TO DATE
     * @input          stoppingPoint, an ancestor of the actual node
     * @semantics      fill the identical array of the node and all its
     *                 ancestors BEFORE stoppingPoint
     ************************************************************************ */
    void updateIdenticalRecursively( InferenceNode* stoppingPoint = NULL );

    /** ************************************************************************
     * invalidateRecursively
     * @input          stoppingPoint, an ancestor of the actual node
     * @input          cat, the symbol category concerned (-1=all)
     * @semantics      switch the partialLikelihood array to the working state
     *                 from this node and all its ancestor BEFORE
     *                 stoppingPoint
     ************************************************************************ */
    void invalidateRecursively( int cat, InferenceNode* stoppingPoint = NULL );

    /** ************************************************************************
     * retrieveRecursively
     * @input          stoppingPoint, an ancestor of the actual node
     * @input          cat, the symbol category concerned (-1=all)
     * @semantics      switch the partialLikelihood array to the saved state
     *                 from this node and all its ancestor BEFORE
     *                 stoppingPoint
     ************************************************************************ */
    void retrieveRecursively( int cat, InferenceNode* stoppingPoint = NULL );

    /** ************************************************************************
     * saveRecursively
     * @input          stoppingPoint, an ancestor of the actual node
     * @input          cat, the symbol category concerned (-1=all)
     * @semantics      all the working array IN USE from this node to stopping
     *                 point becomes the saved state, and partial likelihood
     *                 is switched to the saved state
     ************************************************************************ */
    void saveRecursively( int cat, InferenceNode* stoppingPoint = NULL );


    /** ************************************************************************
     * recursiveCopy
     * @semantics      recursive function used to copy an inference Node
     *                 and its children the reference to the tree (ptree,
     *                 ,... are omitted), the actual ptree is kept,
     *                 new children are created.
     ************************************************************************ */
    void recursiveCopy( InferenceNode * srcNode, InferenceNode * parent );

#ifdef DEBUG1
    /** ************************************************************************
     * checkIdentical
     * @semantics  Function provided for debugging purpose, check the identical
     *             array for all the nodes
     ************************************************************************ */
    void checkIdentical();
#endif

    /** ************************************************************************
     * assignClustersRec
     * @return          the clusters in the tree for this node
     * @semantics       recursive helper for the function
     *                  InferenceTree::assignClusters
     *                  Check conformance of the user's tree with the defined
     *                  clades, initialise ancestor nodes of each cluster
     ************************************************************************ */
     Cluster* assignClustersRec( ClustersTree* clustersTree );


    /** ************************************************************************
     * getCluster
     * @return          the cluster in the tree for this node
     * @semantics       recursive helper for the function
     *                  InferenceTree::checkClusters
     * @preconditions   sequence table loaded
     ************************************************************************ */
    Cluster* getCluster() const;


    /** ************************************************************************
     * findClusterBranches
     * @input      the including cluster
     * @return     sum, the sum of the distance in the branchVector
     * @return     a vector of node (equivalent to branches) with distance
     *             information
     * @semantics  InferenceTree::findClusterBranches helper
     *             recursive call, this function should be called with the
     *             ancestor of the cluster.
     ************************************************************************ */
    double findClusterBranches( vector< pair<InferenceNode*,double> > & branchVector,
                                              ClustersTreeNode* clustersTreeNode );


    /** ************************************************************************
     * findIncludedBranches
     * @return     sum, the sum of the distance in the branchVector
     * @return     a vector of node (equivalent to branches) with distance
     *             information
     * @semantics  This function is similar to the previous one, given a
     *             branch in clustersTreeNode (ie the branch from this node
     *             to its parent), this function will return all
     *             the branches that are in clustersTreeNode and child of
     *             this node. (The branch from this node to its parent is included
     *             as well). In practice, this function just call
     *             the private function findClusterBranchesRec
     ************************************************************************ */
    double findIncludedBranches( vector< pair<InferenceNode*,double> > & branchVector,
                                              ClustersTreeNode* clustersTreeNode );


private:

    /** ************************************************************************
     * getClusterRec
     * @return          the cluster in the tree for this node
     * @semantics       recursive helper for the function
     *                  InferenceTree::checkClusters
     * @preconditions   sequence table loaded
     ************************************************************************ */
    void getClusterRec( Cluster* cluster ) const;

    /** ************************************************************************
     * findClusterBranchesRec
     * @input      the including cluster
     * @return     sum, the sum of the distance in the branchVector
     * @return     a vector of node (equivalent to branches) with distance
     *             information
     * @semantics  InferenceTree::findClusterBranches helper
     *             recursive call
     ************************************************************************ */
    double findClusterBranchesRec( vector< pair<InferenceNode*,double> > & branchVector,
                                   ClustersTreeNode* clustersTreeNode );

protected:


    /** ************************************************************************
     * optimizedLikelihoodLeafParent
     * @semantics  optimized version of the likelihood function for internal
     *            nodes with 2 leaves children
     ************************************************************************ */
    void optimizedLikelihoodLeafParent();


    /** ************************************************************************
     * loafLeafIdentical
     * @semantics  fill the identical array of the leaf nodes.
     *             fill the partialLikelihood array of the leaf nodes.
     *             this function is called ONCE during initialisation(true)
     *             Symbols in sequences are checked during the process
     ************************************************************************ */
    void loafLeaf();

};

#endif //INFERENCENODE_H
