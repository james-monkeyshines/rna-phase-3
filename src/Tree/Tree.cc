#include "Tree/Tree.h"

#include <assert.h>

#include <math.h>

#include <algorithm>
#include <functional>
#include <iostream>

/*tolower function*/
#include <cctype>
#include <cstdlib>

using namespace std;

#include "Tree/TreeMap.h"

#include "Tree/BasicNode.h"


Tree::Tree() {
    root = NULL;
}

Tree::Tree( string stree ) {
    // Create a depth map of the tree
    TreeMap map( stree );
    // Create the nodes
    root = new BasicNode( map );
    root->setParentDistance(-1.0);
    createIndex();
}

Tree::~Tree(){
    if (root){
        delete root;
    }
}

bool Tree::makeRoot( BasicNode* newRoot ){
    if ( newRoot->makeRoot() ){
        root->updateNodesCountRecursively();
        root = newRoot;
        createIndex();
        return true;
    }
    return false;
}


string Tree::toString( bool topologyOnly ) const{
    return ( root->toString( topologyOnly ) );
}
string Tree::toStringCheck( bool topologyOnly ){
    return ( root->toStringCheck( topologyOnly ) );
}


void Tree::createIndex() {
    unsigned int nodeToExpand;

    // construct the nodeRefVect ensure that for all (i,j) / i < j,
    // nodeRefVector[j] is NOT an ancestor of nodeRefVector[i]
    // this will help future computation when the state of an ancestor
    // requires the state of its childs to be known
    nodeRefVector.clear();
    sort();

    if ( root != NULL ) {
        nodeRefVector.push_back( root );
    }
    else {
        return;
    }

    nodeToExpand = 0;
    while ( nodeToExpand != nodeRefVector.size() ) {
        // Expand out the node
        for ( list<BasicNode*>::iterator iter =
                  nodeRefVector[nodeToExpand]->getChildrenList().begin();
              iter != nodeRefVector[nodeToExpand]->getChildrenList().end(); ++iter ) {
            nodeRefVector.push_back( *iter );
        }
        //next node
        ++nodeToExpand;
    }
}

void Tree::move( BasicNode* src, BasicNode* dest ){
    assert(src!=dest);
    assert( !dest->isDescendant(src) );
    BasicNode* parent = src->getParent();
    assert( parent != dest );
    if (parent != dest){
        double distance = src->getParentDistance();
        BasicNode* ancestor = getMostRecentAncestor(parent, dest);
        parent->removeChild( src );
        dest->addChild( src, distance );
        parent->updateNodesCountRecursively(ancestor);
        src->updateNodesCountRecursively();
        createIndex();
    }
}

void Tree::remove( BasicNode* src ){
    BasicNode* parent = src->getParent();
    //cannot remove the root
    assert( parent );
    parent->removeChild( src );
    parent->updateNodesCountRecursively();
    createIndex();
}

void Tree::add( BasicNode* src, BasicNode* dest, double distance ){
    //src should not be a node of the tree
    assert( !dest->isDescendant(src) );
    dest->addChild( src, distance );
    dest->updateNodesCountRecursively();
    createIndex();
}

void Tree::getBranchVector( vector < double > & branchesLength ) const {

    branchesLength.clear();

    vector < BasicNode * >::const_iterator iter = nodeRefVector.begin();
    ++iter;
    while ( iter != nodeRefVector.end() ){
        branchesLength.push_back( (*iter)->getParentDistance() );
        ++iter;
    }
}

void Tree::setBranchVector( const vector < double > & branchesLength ) {
    vector< double >::const_iterator iterLength = branchesLength.begin();
    vector < BasicNode * >::iterator iter = nodeRefVector.begin();
    ++iter;
    while ( iter != nodeRefVector.end() ){
        (*iter)->setParentDistance(*iterLength);
        ++iterLength;
        ++iter;
    }
    assert( iterLength == branchesLength.end() );
}


double Tree::getBranchLength( unsigned int branchId ) const {
    assert( branchId < getNumberBranches() );
    return nodeRefVector[branchId+1]->getParentDistance();
}

BasicNode* Tree::getBranchOrig( unsigned int branchId ){
    return nodeRefVector[branchId+1]->getParent();
}

BasicNode* Tree::getBranchDest( unsigned int branchId ){
    return nodeRefVector[branchId+1];
}

void Tree::setBranchLength( unsigned int branchId, double newDist ) {
    assert( branchId < getNumberBranches() );
    assert( newDist >= 0 );
    nodeRefVector[branchId+1]->setParentDistance( newDist );
}

BasicNode* Tree::getMostRecentAncestor( BasicNode* node1, BasicNode* node2 ) const{
    unsigned int nbLeaves1, nbInternals1;
    unsigned int nbLeaves2, nbInternals2;
    node1->getNodesCount( nbInternals1, nbLeaves1 );
    node2->getNodesCount( nbInternals2, nbLeaves2 );
    while( node1 != node2 ){
        bool change = false;
        //climb from 1 as much as possible without bypassing the most common recent ancestor
        while ( (nbInternals1 < nbInternals2) || (nbLeaves1 < nbLeaves2) ){
            node1 = node1->getParent();
            assert(node1);
            node1->getNodesCount( nbInternals1, nbLeaves1 );
            change = true;
        }
        //climb from 2 as much as possible without bypassing the most common recent ancestor
        while ( (nbInternals2 < nbInternals1) || (nbLeaves2 < nbLeaves1) ){
            node2 = node2->getParent();
            assert(node2);
            node2->getNodesCount( nbInternals2, nbLeaves2 );
            change = true;
        }
        //if it has not been possible to do anything (same number of internal, ...) both can climb one step
        if (!change){
            node1 = node1->getParent();
            assert(node1);
            node1->getNodesCount( nbInternals1, nbLeaves1 );
            node2 = node2->getParent();
            assert(node2);
            node2->getNodesCount( nbInternals2, nbLeaves2 );
        }
    }
    assert(node1);
    return node1;
}

void Tree::checkIndex(){
    cout << "WARNING: checking nodeRefVector" << endl;
    vector<BasicNode*> oldVector = nodeRefVector;
    createIndex();
    assert( oldVector == nodeRefVector );
}

void Tree::checkNodesCount(){
    cout << "WARNING: checking nodes count" << endl;
    for ( vector< BasicNode* >::reverse_iterator iter = nodeRefVector.rbegin();
          iter != nodeRefVector.rend(); ++iter ){
        if ( (*iter)->isLeaf() ){
            assert((*iter)->nbInternalNodes==0);
            assert((*iter)->nbLeaves==1);
        }
        else{
            unsigned int nbInternalNodes = 1;
            unsigned int nbLeaves = 0;
            for ( list<BasicNode*>::iterator childIter = (*iter)->getChildrenList().begin();
                  childIter != (*iter)->getChildrenList().end(); ++childIter ){
                nbInternalNodes += (*childIter)->nbInternalNodes;
                nbLeaves += (*childIter)->nbLeaves;
            }
            assert((*iter)->nbInternalNodes==nbInternalNodes);
            assert((*iter)->nbLeaves==nbLeaves);
        }
    }
}

