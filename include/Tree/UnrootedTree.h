#ifndef UNROOTEDTREE_H
#define UNROOTEDTREE_H

#include "Tree/InferenceTree.h"

class ParametersSet;

class UnrootedTree : virtual public InferenceTree {
private:
    //outgroup can be the name of a specie or of a cluster
    string outgroupName;

protected:
    BasicNode* outgroupNode;
    Cluster* outgroupCluster;

    /** ************************************************************************
     * UnrootedTree
     * @semantics  default constructor of an empty unrooted tree
     *             (used by prototypes)
     ************************************************************************ */
    UnrootedTree();

    /** ************************************************************************
     * UnrootedTree
     * @semantics  common constructor of an empty unrooted tree
     *             (used by descendants)
     ************************************************************************ */
    UnrootedTree(ParametersSet& parameters);

    /** ************************************************************************
     * ~UnrootedTree
     * @semantics  destructor of an unrooted inference tree
     ************************************************************************ */
    virtual ~UnrootedTree();

public:
    /** ************************************************************************
     * unroot
     * @return     true if the tree was rooted
     * @semantics  method used to unroot a rooted tree
     ************************************************************************ */
     bool unroot();

    /** ************************************************************************
     * checkBinary/isBinary
     * @semantics  check that the tree is a strictly bifurcating binary tree
     *             after construction (the root must have 3 children).
     *             checkBinary exit with an error message if this is not the
     *             case
     ************************************************************************ */
    void checkBinary(BasicNode* node = NULL);
    BasicNode* isBinary(BasicNode* node);
    
    /** ************************************************************************
     * installOutgroup
     * @semantics  use the outgroup to root the tree
     ************************************************************************ */
    bool installOutgroup();

    /** ************************************************************************
     * createOutgroup
     * @semantics  initialise outgroupNode/outgroupCluster (with findOutgroup)
     *             and install it (findOutgroup must return true)
     * @return     installOutgroup call result
     ************************************************************************ */
    bool createOutgroup();

    /** ************************************************************************
     * findOutgroup
     * @semantics  look for the outgroupNode/outgroupCluster in the tree
     * @return     a "outgroup" node if cluster or outgroup species in the tree
     *             (ie, a pointer to the leaf node ou outgroupCluster->ancestor)
     ************************************************************************ */
    BasicNode* findOutgroup();

    /** ************************************************************************
     * isOutgroup
     * @return     true if the parameter is the outgroup
     * @semantics  lower case comparison between outgroup and the query
     ************************************************************************ */
    bool isOutgroup( const string& queryOutgroup);

    /** ************************************************************************
     * changeOutgroup
     * @input the outgroup to use to normalize the tree
     * @semantics  change the outgroup specie
     ************************************************************************ */
    void changeOutgroup( const string& outgroup );


     /** ************************************************************************
     * constructFromString
     * @input          the string representation of the tree
     * @semantics      initialise the tree with a user-specified string
     ************************************************************************ */
    virtual void constructFromString( string stringTree );


    /** ************************************************************************
     * constructRandomly
     * @input      the maximum length for the branches during the creation
     * @semantics  create a random initial tree for the species in table
     ************************************************************************ */
    virtual void constructRandomly( double maxBranchLength );

    /** ************************************************************************
     * loadDataAndModel
     * @input           he sequences and the model to use with this tree
     * @preconditions   no data or model already loaded
     * @semantics       this method is redefined in the unrooted tree
     *                  to check whether the outgroup provided by the user
     *                  is sensible
     ************************************************************************ */
    virtual void loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel );

};

#endif
