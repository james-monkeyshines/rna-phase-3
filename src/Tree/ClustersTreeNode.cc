#include "Tree/ClustersTreeNode.h"

#include <algorithm>

#include "Tree/Cluster.h"

#include <assert.h>

using namespace std;

ClustersTreeNode::ClustersTreeNode( Cluster* cluster ){
    this->cluster = cluster;
    this->parent = NULL;
}

ClustersTreeNode::ClustersTreeNode( Cluster* cluster, ClustersTreeNode* parent ){
    this->cluster = cluster;
    this->parent = parent;
}

ClustersTreeNode::~ClustersTreeNode(){
    delete cluster;
    for ( iterator iter = begin(); iter != end(); ++iter ){
        delete *iter;
    }
}


void ClustersTreeNode::makeListRec( vector< ClustersTreeNode* >& listNode ){
    for ( list<ClustersTreeNode*>::const_iterator iter = begin();
          iter != end(); ++iter ){
        (*iter)->makeListRec(listNode);
    }
    listNode.push_back( this );
}

ClustersTreeNode* ClustersTreeNode::findInsertion( unsigned int speciesId ){
    for ( list<ClustersTreeNode*>::const_iterator iter = begin();
          iter != end(); ++iter ){
        if ((*iter)->cluster->contains( speciesId )){
            return (*iter)->findInsertion( speciesId );
        }
    }
    //nothing found? then this cluster is the one
    return this;
}

ClustersTreeNode* ClustersTreeNode::findInsertion( Cluster* cluster ){
    for ( list<ClustersTreeNode*>::const_iterator iter = begin();
          iter != end(); ++iter ){
        if ((*iter)->cluster->contains( *cluster )){
            return (*iter)->findInsertion( cluster );
        }
    }
    //nothing found? then this cluster is the one
    //no compatibility check is done
    return this;
}

void ClustersTreeNode::add( ClustersTreeNode* clustersTreeNode ){
    clustersTreeNode->parent = this;
    push_back(clustersTreeNode);
}

void ClustersTreeNode::remove( ClustersTreeNode* clustersTreeNode ){
    clustersTreeNode->parent = NULL;
    list<ClustersTreeNode*>::iterator iter = find(begin(),end(),clustersTreeNode);
    assert(iter!=end());
    erase(iter);
}
