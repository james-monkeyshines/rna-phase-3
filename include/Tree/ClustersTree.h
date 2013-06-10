#ifndef CLUSTERSTREE_H
#define CLUSTERSTREE_H

#include "Tree/ClustersSet.h"
#include "Tree/InferenceNode.h"

#include <list>
#include <vector>

class Cluster;
class ClustersTreeNode;

/** ***************************************************************************
 * ClustersTree is used to manage a set of coherent clusters/clades
 * They correspond to the user definition of restriction on the topology space
 **************************************************************************** */

class ClustersTree : public ClustersSet{

public:
    /** ************************************************************************
     * ClustersTree
     * @semantics  constructor
     *             create an empty clusters tree with just a root cluster
     *             cluster can be added with the method addcluster
     ************************************************************************ */
    ClustersTree(SequenceTable* table);

    /** ************************************************************************
     * ClustersTree
     * @semantics  constructor from a file all the clusters in the file are
     *             added
     ************************************************************************ */
    ClustersTree( ifstream& fileStream, SequenceTable* table );

    /** ************************************************************************
     * ClustersTree
     * @semantics  destructor of the class, all the ClustersTreeNode created
     *             should be destroyed by the deletion of the element in top
     *             level (cf ClustersTreeNode destructor). All the instances
     *             of Cluster should have been assigned to a ClustersTreeNode
     *             and will be destroyed as well
     ************************************************************************ */
    ~ClustersTree();


    /** ************************************************************************
     * list
     * @return     a vector with all the clustersTreeNode ordered in
     *             inclusive order
     * @semantics  useful to browse the constraints on topology quickly
     ************************************************************************ */
    const vector< ClustersTreeNode* > & getList() const{
        return listNode;
    }

    /** ************************************************************************
     * addCluster
     * @return     the result from insert
     * @semantics  add a cluster to the tree
     ************************************************************************ */
    inline bool addCluster( Cluster* cluster ){
        bool res = insert(cluster);
        if( res ){
            makeList();
        }
        return res;
    }

    /** ************************************************************************
     * findInsertion
     * @return     ClustersTreeNode
     * @semantics  return the smallest cluster(TreeNode) that contains
     *             the given cluster, the given species.
     *             return root if no cluster found
     ************************************************************************ */
    ClustersTreeNode* findInsertion( unsigned int speciesId ) const;
    ClustersTreeNode* findInsertion( Cluster* cluster ) const;


    /** ************************************************************************
     * findByAncestor
     * @return     ClustersTreeNode
     * @semantics  return a cluster tree node whose ancestor matches the
     *             given node (NULL if nothing found)
     ************************************************************************ */
    ClustersTreeNode* findByAncestor( InferenceNode* node ) const;

    inline ClustersTreeNode* getRoot(){
        return root;
    }

    /** ************************************************************************
     * unrootCheck
     * @return     a boolean, false means the tree must be rooted according
     *             to the clusters
     * @semantics  we cannot have two complementary clusters specified
     *             simultaneously with unrooted tree, this function does the
     *             check. 
     ************************************************************************ */
    bool unrootCheck();

protected:
    /** ************************************************************************
     * addTreeClusters
     * @return     the top cluster
     * @semantics  add to the clusters the set of clusters defined by the
     *             topology, (recursive call)
     ************************************************************************ */
    Cluster* addTreeClusters( const string& stringTree, const string& clusterName,
                                 unsigned int& clusterId, SequenceTable* table );


    /** ************************************************************************
     * insert
     * @return     always true, (abort if conflict)
     * @semantics  add a cluster to the set and to the tree structure
     *             this function should not be called outside the constructor
     *             since the list is created at construction
     ************************************************************************ */
    bool insert( Cluster* cluster );


    /** ************************************************************************
     * makeList
     * @semantics  create the list of clusters (if A in B then indexA<indexB)
     *             (recursive call)
     ************************************************************************ */
    void makeList();

    /** ************************************************************************
     * the top level node of the clusters hierarchy, clusters of all
     * species
     ************************************************************************ */
    ClustersTreeNode* root;

    /** ************************************************************************
     * a vector of all the clusters in the set (if A in B then indexA<indexB)
     ************************************************************************ */
    vector< ClustersTreeNode* > listNode;

    /** ************************************************************************
     * keep a pointer to the sequence table
     ************************************************************************ */
    SequenceTable* table;
};

#endif //CLUSTERSTREE_H
