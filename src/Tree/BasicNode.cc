#include "Tree/BasicNode.h"

#include <stdio.h>
#include <assert.h>

#include <iostream>

// STL algorithms class library
#include <algorithm>

#include "Tree/TreeMap.h"

#include "Util/array2D.h"

#include "Tree/Tree.h"

class SequenceTable;


// Constructors
BasicNode::BasicNode() {
    parent = NULL;
    pdistance = -1.0;
    nbInternalNodes = 1;
    nbLeaves = 0;
}

BasicNode::BasicNode( const string& label ) {
    parent = NULL;
    pdistance = -1.0;
    this->label = label;
    nbInternalNodes = 0;
    nbLeaves = 1;
}

BasicNode::BasicNode( const TreeMap& map ){    
    if ( map.getNumberChildren() ){
        parent = NULL;
        for ( unsigned int i = 0; i < map.getNumberChildren() ; ++i ) {
            addChild( new BasicNode( map.getChildMap(i) ), map.getDistance(i) );
        }
        updateNodesCount();
    }
    else{
        parent = NULL;
        nbInternalNodes = 0;
        nbLeaves = 1;
        this->label = map.getLabel();
    }
}

BasicNode::BasicNode( const BasicNode & src ) {
    assert( 0 );
    if(src.nbLeaves){};
}

BasicNode & BasicNode::operator = ( const BasicNode & src ) {
    assert(!src.isLeaf());    
    pdistance = -1.0;
    nbLeaves = src.nbLeaves;
    nbInternalNodes = src.nbInternalNodes;
    childrenList.clear();
    return ( * this );
}


// Destructor
BasicNode::~BasicNode() {
    for ( list< BasicNode* >::iterator iter = childrenList.begin();
        iter != childrenList.end(); ++iter ){
        delete (*iter);
    }
}

bool BasicNode::isDescendant(BasicNode* hypotheticalAncestor) const {
    assert (hypotheticalAncestor!=this);
    BasicNode* temp = this->getParent();
    while( temp != NULL ){
        if (temp == hypotheticalAncestor){
            return true;
        }
        temp = temp->getParent();
    }
    return false;
}


BasicNode* BasicNode::getChild(unsigned int childId){
    assert( childId < childrenList.size() );
    list< BasicNode* >::iterator childIter = childrenList.begin();
    while (childId>0){
        ++childIter;
        --childId;
    }
    return *childIter;
}


const BasicNode* BasicNode::getChild(unsigned int childId) const{
    assert( childId < childrenList.size() );
    list< BasicNode* >::const_iterator childIter = childrenList.begin();
    while (childId>0){
        ++childIter;
        --childId;
    }
    return *childIter;
}

string BasicNode::toString( bool topologyOnly ) const{

    //if leaf node
    if ( isLeaf() ) {
        return label;
    }
    //else
    string tree= "(";
    
    for ( list< BasicNode* >::const_iterator iter = childrenList.begin();
          iter != childrenList.end(); ++iter ){
        //comma separator
        if ( iter != childrenList.begin() ){
            tree += ',';
        }
        tree += (*iter)->toString( topologyOnly );
        if( !topologyOnly ){
            char d[30];
            sprintf( d, ":%.8f", (*iter)->pdistance );
            tree += d;
        }        
    }
    tree += ")";
    return tree;
}


bool BasicNode::operator < ( const BasicNode & other) const{
    unsigned int nbInternalNodes2;
    unsigned int nbLeaves2;
    
    other.getNodesCount(nbInternalNodes2, nbLeaves2);
    if (nbLeaves != nbLeaves2){
        return (nbLeaves<nbLeaves2);
    }
    //if equal compare with the number of internal nodes
    if (nbInternalNodes != nbInternalNodes2){
        return (nbInternalNodes<nbInternalNodes2);
    }
    //if still equal compare the number of children
    if ( getNumberChildren() != other.getNumberChildren() ){
        return ( getNumberChildren() < other.getNumberChildren() );
    }
    //then compare each child
    list<BasicNode*>::const_iterator iter1 = childrenList.begin();
    list<BasicNode*>::const_iterator iter2 = other.getChildrenList().begin();
    while( iter1 != childrenList.end() ){      
        if( *(*iter1) < *(*iter2) ){
            return true;
        }
        if( *(*iter2) < *(*iter1) ){
            return false;
        }
        ++iter1;
        ++iter2;
    }
    // at last use the topology, if labels are equal then the nodes
    // are really equal and the answer can be true or false
    return ( toString() < other.toString() );
}

bool BasicNode::makeRoot(){
    assert( !isLeaf() );
    return makeParent( NULL );
}

bool BasicNode::makeParent( BasicNode * newParent ) {
    //add the old parent (if it exists) to the list of child
    //and iterate the process while the root is not reached
    if ( parent != NULL ){
        parent->makeParent( this );
        addChild( parent, pdistance );
    }
    //the new parent must have been one of the children so remove it
    if (newParent != NULL){
        list< BasicNode* >::iterator childIter =
            find( childrenList.begin(), childrenList.end(), newParent );
        childrenList.erase( childIter );
    }
    else{
        //end of the recusive call, the old root is reached
        pdistance = -1.0;
    }
    parent = newParent;
    return true;
}              

void BasicNode::getNodesCount( unsigned int & nbInternalNodes,
                                    unsigned int & nbLeaves ) const{
    nbLeaves = this->nbLeaves;
    nbInternalNodes = this->nbInternalNodes;
}

void BasicNode::updateNodesCount(){
    unsigned int nbInternalNodesChild;
    unsigned int nbLeavesChild;
    
    if ( isLeaf() ){
        nbInternalNodes = 0;
        nbLeaves = 1;
    }
    else{
        nbInternalNodes = 1;
        nbLeaves = 0;
        for ( list<BasicNode*>::iterator iter = childrenList.begin();
              iter != childrenList.end(); ++iter ){
            (*iter)->getNodesCount(nbInternalNodesChild, nbLeavesChild);
            nbInternalNodes += nbInternalNodesChild;
            nbLeaves += nbLeavesChild;
        }
    }
}

void BasicNode::updateNodesCountRecursively(BasicNode* stoppingNode){
    BasicNode* tempNode = this;
    while (tempNode != stoppingNode){
        assert(tempNode);
        tempNode->updateNodesCount();
        tempNode = tempNode->getParent();
    }
}


void BasicNode::sort(){
    //sort each child first
    for ( list<BasicNode*>::iterator childIter = childrenList.begin();
          childIter != childrenList.end(); ++childIter ){
        (*childIter)->sort();
    }
    //then sort the children of this node (the sorting algorithm is not
    //sophisticated but we expect 2 children only)
    if( !isLeaf() ){
        BasicNode* temp;
        for ( list<BasicNode*>::iterator iter1 = childrenList.begin();
              iter1 != childrenList.end(); ++iter1 ){
            for ( list<BasicNode*>::iterator iter2 = --childrenList.end();
                iter2 != iter1; --iter2 ){
                //use BasicNode::operator <
                if ( *(*iter2) < *(*iter1) ){
                    temp = *iter2;
                    *iter2 = *iter1;
                    *iter1 = temp;
                }
            }
        }
    }
}

void BasicNode::removeChild( BasicNode* childNode ){
    list< BasicNode* >::iterator childIter =
            find( childrenList.begin(), childrenList.end(), childNode );
    childrenList.erase( childIter );
    childNode->parent = NULL;
}

void BasicNode::addChild( BasicNode* childNode, double distance ){    
    childrenList.push_back( childNode );
    childNode->setParent(this);
    childNode->setParentDistance( distance );
}


string BasicNode::toStringCheck( bool topologyOnly ) const{

    //if leaf node
    if ( isLeaf() ) {
        return label;
    }
    //else
    string tree= "(";
    char d[30];
    for ( list< BasicNode* >::const_iterator iter = childrenList.begin();
          iter != childrenList.end(); ++iter ){
        //comma separator
        if ( iter != childrenList.begin() ){
            tree += ",";
        }
        tree += (*iter)->toStringCheck( topologyOnly );
        if( !topologyOnly ){
            sprintf( d, ":%.8f", (*iter)->pdistance );
            tree += d;
        }
        unsigned int a,b;
        (*iter)->getNodesCount( a, b );
        sprintf( d, ":(%d,%d)", a, b );
        tree += d;
        
    }
    tree += ")";
    return tree;
}
