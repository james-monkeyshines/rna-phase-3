#ifndef TREE_H
#define TREE_H

#include "configfix.h"

#include <string>
#include <vector>

#include "Tree/BasicNode.h"

using namespace std;

class Tree {
protected:
    /** the root of the tree                                                  */
    BasicNode * root;

    /** ************************************************************************
     * sort
     * @semantics   sort the tree to have a common string representation
     *              for equivalent tree
     ************************************************************************ */
    inline void sort(){
        root->sort();
    }

public:

    /** an array pointing at all the nodes (see createIndex)                  */
    vector < BasicNode * > nodeRefVector;

    /** ************************************************************************
     * Tree
     * @semantics  default constructor of an empty tree
     ************************************************************************ */
    Tree();

    /** ************************************************************************
     * Tree
     * @input      the string representation of the tree
     * @semantics  constructor from a string
     ************************************************************************ */
    Tree( string stree );

    /** ************************************************************************
     * ~Tree
     * @semantics  destructor of a tree
     ************************************************************************ */
    virtual ~Tree();


    /** ************************************************************************
     * getRoot
     * @semantics  return the root of the tree (const version)
     ************************************************************************ */
    inline BasicNode* getRoot() const{
        return root;
    }

    /** ************************************************************************
     * setRoot
     * @semantics  set a root for the tree
     ************************************************************************ */
    inline void setRoot( BasicNode* root ){
        this->root = root;
    }

    /** ************************************************************************
     * toString
     * @input      topologyOnly, set to false if the representation must include
     *             lengths information
     * @semantics  return the string representation of this tree. The tree
     *             should be sorted before.
     ************************************************************************ */
    string toString( bool topologyOnly = false ) const;

    /** ************************************************************************
     * toStringCheck
     * @input      topologyOnly, set to false if the representation must include
     *             lengths information
     * @semantics  return the string representation of this tree.
     ************************************************************************ */
    string toStringCheck( bool topologyOnly = false );

    /** ************************************************************************
     * getBranchVector
     * @input      branchesLength, the vector to fill
     * @semantics  return the lengths of the branches of the tree (preorder
     *             traversal), BECAREFUL, the tree is not sorted before.
     *             Use createIndex
     ************************************************************************ */
    void getBranchVector( vector < double > & branchesLength ) const;

    /** ************************************************************************
     * setBranchVector
     * @input      branchesLength, the vector with the new lengths
     * @semantics  change the lengths of all the branches
     ************************************************************************ */
    virtual void setBranchVector( const vector < double > & branchesLength );


    /** ************************************************************************
     * move
     * @input         the BasicNode to move, the new father
     * @semantics     move a node to another part of the tree
     *                the distance to the parent remains the same
     * @precondition  no cycle allowed, src and dest in the tree, ...
     ************************************************************************ */
    void move( BasicNode* src, BasicNode* dest );


    /** ************************************************************************
     * remove
     * @input         the BasicNode to remove
     * @semantics     remove a node from the tree
     ************************************************************************ */
    void remove( BasicNode* src );


    /** ************************************************************************
     * add
     * @input         the BasicNode to add, the new father
     * @semantics     add a node to the tree
     * @precondition  no cycle allowed, dest in the tree, ...
     ************************************************************************ */
    void add( BasicNode* src, BasicNode* dest, double distance );


    /** ************************************************************************
     * makeRoot
     * @input         the new root for the tree
     * @semantics     change the root of the tree
     * @precondition  newRoot is an internal node
     ************************************************************************ */
    bool makeRoot( BasicNode* newRoot );


    /** ************************************************************************
     * getBranchLength
     * @input         the ID of the branch
     * @semantics     return the length of the specified branch
     * @precondition  0 <= branchId < numberBranches
     ************************************************************************ */
    double getBranchLength( unsigned int branchId ) const;

    /** ************************************************************************
     * getBranchOrig
     * @input         the ID of the branch
     * @semantics     return a node incident to the edge (the father)
     * @precondition  0 <= branchId < numberBranches
     ************************************************************************ */
     BasicNode* getBranchOrig( unsigned int branchId );

    /** ************************************************************************
     * getBranchDest
     * @input         the ID of the branch
     * @semantics     return a node incident to the edge (the child)
     * @precondition  0 <= branchId < numberBranches
     ************************************************************************ */
    BasicNode* getBranchDest( unsigned int branchId );

    /** ************************************************************************
     * setBranchLength
     * @input         the ID of the branch
     * @input         the new length
     * @semantics     change the length of the specified branch
     * @precondition  0 <= branchId < numberBranches
     ************************************************************************ */
    void setBranchLength( unsigned int branchId, double newLength );


    /** ************************************************************************
     * getMostRecentAncestor
     * @input          node1, node2 the two nodes
     * @semantics      return the most common recent ancestor of the two nodes
     * @preconditions  nodes are in the same consistent tree
     ************************************************************************ */
     BasicNode* getMostRecentAncestor( BasicNode* node1, BasicNode* node2 ) const;

    inline unsigned int getNumberNodes() const {
        return nodeRefVector.size();
    }

    inline unsigned int getNumberBranches() const {
        return ( nodeRefVector.size() - 1);
    }

    inline unsigned int getNumberTips() const {
        unsigned int numberNodes, numberLeaves;
        root->getNodesCount( numberNodes, numberLeaves );
        return ( numberLeaves );
    }

    /** ************************************************************************
     * checkIndex//checkNodesCount
     * @semantics   debug function, can be used to check that the tree is left
     *              in a correct state while performing unprotected but
     *              efficient low-level node operation
     ************************************************************************ */
    void checkIndex();
    void checkNodesCount();

    void createIndex();

};

#endif //TREE_H




