#ifndef BASICNODE_H
#define BASICNODE_H

#include <string>
#include <list>

using namespace std;

class Tree;

class TreeMap;

class Model;


class BasicNode {
private:

    /** ************************************************************************
     * children of the node
     ************************************************************************ */
    list< BasicNode * > childrenList;

    /** ************************************************************************
     * BasicNode
     * @semantic  copy constructor, access barred
     ************************************************************************ */
    BasicNode ( const BasicNode & );

public:
    /** ************************************************************************
     * the number of internal nodes below that node (this one included)
     ************************************************************************ */
    unsigned int nbInternalNodes;

    /** ************************************************************************
     * the number of leaves below that node (1 if the node is a leaf)
     ************************************************************************ */
    unsigned int nbLeaves;


protected:

    /** ************************************************************************
     * operator=
     * @semantic  = operator, specific use for InferenceNode::operator=
     *            to be used with an internal node only,
     ************************************************************************ */
    BasicNode & operator = ( const BasicNode & );

    /** ************************************************************************
     * the label of the node if any
     ************************************************************************ */
    string label;

    /** ************************************************************************
     * the parent of the node
     ************************************************************************ */
    BasicNode * parent;

    /** ************************************************************************
     * the distance to the parent
     ************************************************************************ */
    double pdistance;


public:

    /** ************************************************************************
     * cnumber
     * an index for the leaves node (used in InferenceTree::toStringNumbered)
     * @semantics    we define this field here even though there is no
     *               SequenceTable to make a correspondance between species
     *               name and id. It is useful for the consensus method
     *               which read tree in number format
     ************************************************************************ */
    int cnumber;


    /** ************************************************************************
     * BasicNode
     * @semantic  constructor of a detached node
     ************************************************************************ */
    BasicNode();

    /** ************************************************************************
     * BasicNode
     * @semantic  constructor of a leaf node
     ************************************************************************ */
    BasicNode( const string& label );

    /** ************************************************************************
     * BasicNode
     * @semantic  constructor of an arborescence of node
     ************************************************************************ */
    BasicNode( const TreeMap& map );

private:
    /** ************************************************************************
     * initBasicNode
     * @semantics   constructor primitive
     ************************************************************************ */
    void initBasicNode();

public:
    /** ************************************************************************
     * ~BasicNode
     * @semantic  destructor
     ************************************************************************ */
    virtual ~BasicNode();

    /** ************************************************************************
     * toString
     * @input      topologyOnly, set to true to remove branch lengths in the
     *             outputs
     * @semantics  This method outputs the string representation of the tree
     *             from this node
     ************************************************************************ */
    string toString( bool topologyOnly = false ) const;



    /** ************************************************************************
     * toStringCheck
     * @input      topologyOnly, set to true to remove branch lengths in the
     *             outputs
     * @semantics  This method outputs the string representation of the tree
     *             from this node, nbInternalNodes and nbLeaves are printed
     *             in order to check the output.
     ************************************************************************ */
    string toStringCheck( bool topologyOnly = false ) const;


    /** ************************************************************************
     * operator <
     * @semantics  ordering over the topology space : the aim is to have a
     *             unique string representation for equivalent tree
     ************************************************************************ */
    bool operator < ( const BasicNode & ) const;

    /** ************************************************************************
     * isLeaf
     * @semantics  is the node a leaf (ie. no child) ?
     ************************************************************************ */
    inline bool isLeaf() const {
        return ( childrenList.empty() );
    }

    /** ************************************************************************
     * isDescendant
     * @return     true if the node is a descendant of src
     * @semantics  search if src is one of the ancestors
     ************************************************************************ */
    bool isDescendant(BasicNode* hypotheticalAncestor) const;

    /** ************************************************************************
     * getLabel
     * @semantics  the label of the node (should be a leaf node)
     ************************************************************************ */
    inline const string& getLabel() const {
        return ( label );
    }

    /** ************************************************************************
     * changeLabel
     * @semantics  change the label of the node (should be a leaf node)
     ************************************************************************ */
    inline void changeLabel( const string& newLabel ) {
        label = newLabel;
    }

    /** ************************************************************************
     * getChildrenList
     * @semantics  return the vector of children
     ************************************************************************ */
    inline list< BasicNode* >& getChildrenList(){
        return childrenList;
    }
    inline const list< BasicNode* >& getChildrenList() const{
        return childrenList;
    }

    /** ************************************************************************
     * getChild
     * @input      the Id of the child (warning, this can change)
     * @semantics  return the pointer to the child
     ************************************************************************ */
    BasicNode* getChild(unsigned int childId);
    const BasicNode* getChild(unsigned int childId) const;

    /** ************************************************************************
     * getChildDistance
     * @input      the Id of the child (warning, this can change)
     * @semantics  return the distance to the given childId
     ************************************************************************ */
    inline double getChildDistance(unsigned int childId) const{
        return getChild( childId )->getParentDistance();
    }

    /** ************************************************************************
     * getNumberChildren
     * @semantics  return the number of children
     ************************************************************************ */
    inline unsigned int getNumberChildren() const{
        return(childrenList.size());
    }

    /** ************************************************************************
     * getParent
     * @semantics  return a pointer to the parent node
     ************************************************************************ */
    inline BasicNode* getParent() const{
        return parent;
    }

    inline void setParent( BasicNode* newParent) {
        parent = newParent;
    }


    /** ************************************************************************
     * getParentDistance
     * @semantics  return the distance to the parent
     ************************************************************************ */
    inline double getParentDistance() const{
        return pdistance;
    }


     /** ************************************************************************
     * setParentDistance
     * @input      the new distance
     * @semantics  change the distance to the parent
     ************************************************************************ */
    inline void setParentDistance( double newDistance ){
        pdistance = newDistance;
    }

    /** ************************************************************************
     * getNodesCount
     * @return     the number of internal nodes beyond this node, the number
     *             of leaves beyond this node (the node is included in the
     *             count)
     * @semantics  return the number of internal nodes and leaves in the tree
     *             from this node to the input/output parameters.
     *             This is a recursive call. Use this function with the root
     *             of the tree to retrieve the number of internal nodes and
     *             the number of leaves in the tree
     ************************************************************************ */
    void getNodesCount( unsigned int & nbInternalNodes,
                             unsigned int & nbLeaves ) const;


    /** ************************************************************************
     * sort
     * @semantics   sort the nodes according to the size of the inner trees
     ************************************************************************ */
    void sort();


    /** ************************************************************************
     * makeRoot
     * @preconditions  The node is not a leaf node
     * @semantics      The internal node becomes the new root of the tree.
     *                 nbInternalNodesCounter and nbLeaves are not updated,
     *                 use Tree::makeRoot(BasicNode*) instead
     ************************************************************************ */
    bool makeRoot();

    /** ************************************************************************
     * removeChild
     * @preconditions  child is a child !
     * @semantics      detach the child and its successor from this node.
     *                 Use updateNodesCountRecursively after this call.
     ************************************************************************ */
    void removeChild( BasicNode* childNode );

    /** ************************************************************************
     * addChild
     * @input          child, the new child
     * @input          distance, the distance to the new child
     * @preconditions  The node is not attached and no cycle
     * @semantics      attach a newChild. Use updateNodesCountRecursively
     *                 after this call
     ************************************************************************ */
    void addChild( BasicNode* childNode, double distance );

    /** ************************************************************************
     * makeparent
     * @semantics  Recursive function used to change the root of the tree in
     *             makeRoot. While the old root is not reached, the current
     *             parent become a child and the new parent is the one given in
     *             parameter. This function is then called for the old parent
     *             so as to :
     *             1) let it know its child (newParent) is its parent now.
     *             2) make its own parent its child
     ************************************************************************ */
    bool makeParent( BasicNode * newParent );

    /** ************************************************************************
     * updateNodesCount
     * @semantics  Update the value of the parameters nbInternalNodes and
     *             nbLeavesNodes
     ************************************************************************ */
    void updateNodesCount();

    /** ************************************************************************
     * updateNodesCountRecursively
     * @semantics  Update the value of the parameters nbInternalNodes and
     *             nbLeavesNodes from this node to the stoppingNode.
     *             (root by default)
     ************************************************************************ */
    void updateNodesCountRecursively(BasicNode* stoppingNode = NULL);

};


#endif
