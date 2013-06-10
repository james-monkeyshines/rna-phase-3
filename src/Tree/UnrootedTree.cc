#include "Tree/UnrootedTree.h"

#include "Util/ParametersSet.h"

#include <iostream>
#include <algorithm>
#include <assert.h>

#include "Tree/Cluster.h"
#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

UnrootedTree::UnrootedTree() : InferenceTree(){
}

UnrootedTree::UnrootedTree(ParametersSet& parameters) : InferenceTree(){
    if (parameters.findParameter("Clusters file")){
        loadClusters(parameters.stringParameter("Clusters file"));
    }
    changeOutgroup( parameters.stringParameter("Outgroup") );
    assert( outgroupName != "" );
}

UnrootedTree::~UnrootedTree(){
}

void UnrootedTree::checkBinary(BasicNode* node){
    BasicNode* errorNode = isBinary(node);
    if (errorNode){
        if (errorNode==getRoot()){
            cerr << "error while creating a tree, the root should have 3 children... aborting" << endl;
            cerr << root->toString(true) << endl;
	}
        else{
            cerr << "error while creating a tree, the tree is not binary..." << endl
                 << errorNode->toString(true) << endl
                 << "aborting" << endl;
        }
        exit(EXIT_FAILURE);
    }
}

BasicNode* UnrootedTree::isBinary(BasicNode* node){
    if (node == NULL){
        node = getRoot();
        if (node->getNumberChildren() != 3){
        //bug mandrake gcc 3.3.2... do not remove
            unsigned int numberChildren = root->getNumberChildren();
            if (root->getNumberChildren()==3){
                cerr << "Compilation issue with PHASE and linux (mandrake?)... please report the error and ask for assistance" << endl;
            }
            return node;
        }
    }
    else{
        //if it is a leaf exit
        if (!node->getNumberChildren()) return NULL;
        if (node->getNumberChildren() != 2){
            return node;
        }
    }
    for(list<BasicNode*>::iterator iter = (node->getChildrenList()).begin();
        iter != (node->getChildrenList()).end(); ++iter){
        BasicNode* errorNode = isBinary(*iter);
        if (errorNode) return errorNode;        
    }
    return NULL;
}

bool UnrootedTree::unroot(){
    //if relevant, make sure clusters have been assigned before calling unroot
    
    //if ( root->getNumberChildren() != 2 ) return false;
    //WARNING: strange issue with the previous line and gcc 3.3.2 on mandrake distrib
    //function used to return false even though the root had exactly two children.
    //the command was tweaked and it seems to work fine now.
    BasicNode* node = getRoot();
    unsigned int numberChildren = 0;
    unsigned int nbc = node->getNumberChildren();
    if (nbc == 2){
        numberChildren=root->getNumberChildren();
        if (numberChildren!=2){
		    cerr << "Compilation issue with PHASE and linux (mandrake?)... please report the error and ask for assistance" << endl;
		    exit(EXIT_FAILURE);
        }
    }
    else{
        numberChildren=root->getNumberChildren();
        if (numberChildren==2){
		    cerr << "Compilation issue with PHASE and linux (mandrake?)... please report the error and ask for assistance" << endl;
		    exit(EXIT_FAILURE);
        }
        return false; //three or more children, do not unroot
    }
    
    InferenceNode* removedNode = (InferenceNode*)(root->getChild(0));
    InferenceNode* keepedNode = (InferenceNode*)(root->getChild(1));
    if ( removedNode->isLeaf() ){
        //this is the STL swap algorithm, not the InferenceTree::swap
        std::swap(removedNode, keepedNode);
    }
    else{
        if ( clustersTree && clustersTree->findByAncestor(removedNode) ){
            //this is the STL swap algorithm, not the InferenceTree::swap
            std::swap(removedNode, keepedNode);
            assert(!clustersTree->findByAncestor(removedNode));
        }
    }
    if( removedNode->isLeaf() || (clustersTree && clustersTree->findByAncestor(removedNode)) ){
        cerr << "Sorry, cannot unroot your tree as expected." << endl;
        exit(EXIT_FAILURE);
    }


    double newDist = removedNode->getParentDistance()+
                     keepedNode->getParentDistance();
    keepedNode->setParentDistance(newDist);
    while ( removedNode->getNumberChildren() ){
        move( removedNode->getChild(0), (InferenceNode*)root);
    }
    remove(removedNode);
    delete(removedNode);
    createIndex();
    return true;
}

void UnrootedTree::changeOutgroup( const string& outgroupName ){
    this->outgroupName = outgroupName;
    transform<string::iterator,string::iterator,int(*)(int)>
        ( (this->outgroupName).begin(), (this->outgroupName).end(),
          (this->outgroupName).begin(), tolower );
    outgroupNode = NULL;
    outgroupCluster = NULL;
}

bool UnrootedTree::createOutgroup() {
    findOutgroup();
    assert( outgroupNode );
    return installOutgroup();
}



BasicNode* UnrootedTree::findOutgroup(){
    outgroupNode = NULL;
    outgroupCluster = NULL;
    int speciesId = ptable->findSequence(outgroupName);
    if (speciesId!=-1){
        vector<BasicNode*>::reverse_iterator iter = nodeRefVector.rbegin();
        while ( !outgroupNode && iter!=nodeRefVector.rend() ){
            if ( (*iter)->isLeaf() ){
                string leafLabel( (*iter)->getLabel() );
                transform<string::iterator,string::iterator,int(*)(int)>
                ( leafLabel.begin(), leafLabel.end(), leafLabel.begin(), tolower );
                if ( outgroupName == leafLabel ) {
                    outgroupNode = (*iter);
                }
            }
            ++iter;
        }
    }
    else{
        outgroupCluster = clustersTree->getCluster(outgroupName);
        assert(outgroupCluster);
        outgroupNode = outgroupCluster->getAncestor();
    }
    return outgroupNode;
}

bool UnrootedTree::installOutgroup() {
    //if the outgroup is a cluster, the outgroup node
    //might have changed
    if (outgroupCluster){
        outgroupNode = outgroupCluster->getAncestor();
    }
    assert(outgroupNode);
    assert(outgroupNode->getParent());

    if (clustersTree){
        clustersTree->getRoot()->cluster->setAncestor( outgroupNode->getParent() );
    }
    if ( outgroupNode->getParent() == root ) {
        return false;
    }
    else {
        makeRoot( (InferenceNode*)(outgroupNode->getParent()) );
        createIndex();
    }
    return true;
}

bool UnrootedTree::isOutgroup(const string & queryOutgroup){
    string test(queryOutgroup);
    transform<string::iterator,string::iterator,int(*)(int)>
        ( test.begin(), test.end(), test.begin(), tolower );
    return (test == outgroupName);
}

void UnrootedTree::constructFromString( string stringTree ){
    InferenceTree::constructFromString( stringTree );
    if ( unroot() ){
        cerr << "WARNING: your tree has been unrooted" << endl;
    }
    //unrooted trees do not have to be strictly bifurcating
    //checkBinary();
    createOutgroup();
    //branch lengths must be positive
    vector<BasicNode*>::iterator iter = nodeRefVector.begin();
    while ( (++iter) != nodeRefVector.end() ){
        if ( (*iter)->getParentDistance() < 0.0 ){
            (*iter)->setParentDistance( 0.1 );
        }
        if ( (*iter)->getParentDistance() == 0.0 ){
            (*iter)->setParentDistance( 1e-8 );
        }
    }
}

void UnrootedTree::constructRandomly( double maxBranchLength ){
    InferenceTree::constructRandomly( maxBranchLength );
    unroot();
    createOutgroup();
}


void UnrootedTree::loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel ){
    InferenceTree::loadDataAndModel( psequenceTable, pmodel );
    //check that the clusters list is not incompatible with an unrooted tree
    if (clustersTree){
        if ( !clustersTree->unrootCheck() ){
            cerr << "Sorry, when unrooted tree are used you cannot specify complementary clusters" << endl
                 << "Remove \"" << clustersTree->getRoot()->front()->cluster->getName()
                 << "\" or \"" << clustersTree->getRoot()->back()->cluster->getName() << '"' << endl;
            exit(EXIT_FAILURE);
        }
    }
    //check the outgroup now that we have the sequence table
    int speciesId = ptable->findSequence(outgroupName);
    if (speciesId==-1){
        //if the outgroup is not a species name, it might be a cluster name
        if( !clustersTree ){
            cerr << "Error : outgroup taxa \"" << outgroupName << "\" unknown" << endl;
            exit(EXIT_FAILURE);
        }
        outgroupCluster = clustersTree->getCluster(outgroupName);
        if( !outgroupCluster ){
            cerr << "Error : outgroup taxa (or cluster) \"" << outgroupName << "\" unknown" << endl;
            exit(EXIT_FAILURE);
        }
    }
    //check that outgroup does not contradict the clusters
    if (clustersTree){
        Cluster* testCluster;
        if (speciesId==-1){
            testCluster = outgroupCluster->dup();
        }
        else{
            testCluster = new Cluster(ptable->getNumberSpecies());
            testCluster->addSpecies(speciesId);
        }
        testCluster->invert();
        if (!clustersTree->compatible( *testCluster )){
            cerr << "Error : outgroup taxa \"" << outgroupName << "\" not compatible with your list of clusters" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

