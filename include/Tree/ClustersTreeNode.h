#ifndef CLUSTERSTREENODE_H
#define CLUSTERSTREENODE_H

#include <list>
#include <vector>

class Cluster;

using namespace std;

/** ****************************************************************************
 * the node used to build a hierarchy of Clusters in ClustersTreeSet
 **************************************************************************** */
class ClustersTreeNode : public list<ClustersTreeNode*>{
public:
    /** ************************************************************************
     * ClustersTreeNode
     * constructor of the class
     * @input           the cluster
     * @preconditions   cluster != null
     * @semantics       cluster is be deleted during destruction
     *                  of the instance, children are deleted too
     ************************************************************************ */
    ClustersTreeNode( Cluster* cluster );

    /** ************************************************************************
     * ClustersTreeNode
     * constructor of the class
     * @input           the cluster
     * @input           a pointer to the parent node
     * @preconditions   cluster != null
     * @semantics       cluster is be deleted during destruction
     *                  of the instance, children are deleted too
     ************************************************************************ */
    ClustersTreeNode( Cluster* cluster, ClustersTreeNode* parent );

    /** ************************************************************************
     * ~ClustersTreeNode
     * destructor of the class
     * @semantics       delete the instance of Cluster and destroy the children
     ************************************************************************ */
    ~ClustersTreeNode();

    /** ************************************************************************
     * makeListRec
     * @semantics  create the list of underneath clusters (if A in B then
     *             indexA<indexB)
     *             recursive call, helper for ClustersTree::makeList()
     ************************************************************************ */
    void makeListRec( vector< ClustersTreeNode* >& listNode );

    /** ************************************************************************
     * add
     * @semantics  add clustersTreeNode to the list of children
     ************************************************************************ */
    void add( ClustersTreeNode* clustersTreeNode );

    /** ************************************************************************
     * remove
     * @semantics  remove clustersTreeNode from the list of children
     ************************************************************************ */
    void remove( ClustersTreeNode* clustersTreeNode );

     /** ************************************************************************
     * remove
     * @semantics  remove clustersTreeNode from the list of children
     ************************************************************************ */
    inline ClustersTreeNode* getParent() const{
        return parent;
    }


    /** ************************************************************************
     * the contained cluster
     ************************************************************************ */
    Cluster* cluster;

    /** ************************************************************************
     * findInsertion
     * @semantics  return the smallest cluster containing the species/ cluster
     *             WARNING: No check performed (for cluster insertion), return
     *             an incompatible cluster if the cluster cannot be inserted
     ************************************************************************ */
    ClustersTreeNode* findInsertion( unsigned int speciesId );
    ClustersTreeNode* findInsertion( Cluster* cluster );

    
    
   /** ************************************************************************
     * a pointer to the parent
     ************************************************************************ */
    ClustersTreeNode* parent;
    
protected:

};
#endif //CLUSTERSTREENODE_H
