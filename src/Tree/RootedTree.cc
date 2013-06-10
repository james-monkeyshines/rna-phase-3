#include "Tree/RootedTree.h"

#include "Util/ParametersSet.h"

#include <iostream>
#include <assert.h>

#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

RootedTree::RootedTree() : InferenceTree(){
}

RootedTree::RootedTree(ParametersSet& parameters) : InferenceTree(){
    if (parameters.findParameter("Clusters file")){
        loadClusters(parameters.stringParameter("Clusters file"));
    }
    //if an outgroup is defined with rooted tree, use it to define
    //a cluster
    if (parameters.findParameter("Outgroup")){
        rooting = parameters.stringParameter("Outgroup");
    }
}

RootedTree::~RootedTree(){
}

void RootedTree::constructFromString( string stringTree ){
    InferenceTree::constructFromString( stringTree );
    //rooted trees do not have to be strictly bifurcating
    //checkBinary();
    //but root must have two children
    if (root->getNumberChildren()!=2){
        cerr << "Error, two children are expected for the root, the given tree is invalid:" << endl;
        cerr << stringTree << endl;
        exit(EXIT_FAILURE);
    }
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

void RootedTree::constructRandomly( double maxBranchLength ){
    InferenceTree::constructRandomly( maxBranchLength );
}

void RootedTree::checkBinary(BasicNode* node){
    BasicNode* errorNode = isBinary(node);
    if (errorNode){
        cerr << "error while creating a tree, the tree is not binary..." << endl
             << errorNode->toString(true) << endl
             << "aborting" << endl;
        exit(EXIT_FAILURE);
    }
}


BasicNode* RootedTree::isBinary(BasicNode* node){
    if (node == NULL){
        node = root;
    }
    //if it is a leaf exit
    if (!node->getNumberChildren()) return NULL;
    //otherwise, two children expected
    if (node->getNumberChildren() != 2){
        return node;
    }
    //propagate recursive call and return the answer
    for(list<BasicNode*>::iterator iter = (node->getChildrenList()).begin();
        iter != (node->getChildrenList()).end(); ++iter){
        BasicNode* errorNode = isBinary(*iter);
        if (errorNode) return errorNode;
    }
    return NULL;
}


double RootedTree::findBranches( vector< pair<InferenceNode*,double> > & branchVector,
                       double distanceMax, InferenceNode* excluded, bool crossOnly, ClustersTreeNode* orig ) const{
    double sum = 0.0;
    BasicNode* startingNode = root;
    if (orig!=NULL){
        assert(clustersTree);
        startingNode = orig->cluster->getAncestor();
        assert( excluded != startingNode );
    }
    if ( startingNode != root ){
        if (!crossOnly){
            branchVector.push_back( pair<InferenceNode*,double>((InferenceNode*)startingNode, startingNode->getParentDistance() ) );
            sum += startingNode->getParentDistance();
        }
    }
    for ( list<BasicNode*>::iterator childIter = (startingNode->getChildrenList()).begin();
          childIter != (startingNode->getChildrenList()).end(); ++childIter ){
        if ( *childIter != excluded ){
            sum += findBranchesAux( branchVector, distanceMax,
                                excluded, crossOnly, (InferenceNode*)*childIter );
        }
    }
    return sum;
}




double RootedTree::findBranchesAux( vector< pair<InferenceNode*,double> > & branchVector,
                       double distanceLeft, InferenceNode* excluded, bool crossOnly,
                       InferenceNode* node ) const{
    double sum;
    double parentDistance = node->getParentDistance();
    //termination condition
    if ( parentDistance >= distanceLeft) {
        sum = distanceLeft;
        branchVector.push_back( pair<InferenceNode*,double>((InferenceNode*)node, distanceLeft) );
    }
    else{
        sum = 0.0;
        if (!crossOnly){
            branchVector.push_back( pair<InferenceNode*,double>((InferenceNode*)node, parentDistance) );
            sum += parentDistance;
        }
        if ( !clustersTree || !clustersTree->findByAncestor(node) ){
            for ( list<BasicNode*>::iterator childIter = (node->getChildrenList()).begin();
                  childIter != (node->getChildrenList()).end(); ++childIter ){
                if ( *childIter != excluded ){
                    sum += findBranchesAux( branchVector, distanceLeft-parentDistance, excluded, crossOnly, (InferenceNode*)*childIter );
                }
            }
        }
    }
    return sum;
}

void RootedTree::loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel ){
    InferenceTree::loadDataAndModel( psequenceTable, pmodel);
    //if a root was clearly defined
    if (rooting!=""){
        int speciesId = psequenceTable->findSequence( rooting );
        Cluster* rootCluster = new Cluster(psequenceTable->getNumberSpecies(), "compl_"+rooting);
        //if there is just one species for the outgroup:
        if (speciesId != -1){
            //if no clusters were defined, initialise a cluster tree from scratch
            if (!clustersTree){
                clustersTree = new ClustersTree( psequenceTable );
            }
            //prepare the (maybe unique) cluster
            rootCluster->addSpecies( (unsigned int)speciesId );
            rootCluster->invert();
        }
        else{
            if (!clustersTree){
                cerr << "Sorry, species " << rooting << " not found in your datafile" << endl;
                exit(EXIT_FAILURE);
            }
            else{
                Cluster* outgroupCluster = clustersTree->getCluster( rooting );
                if (!outgroupCluster){
                    cerr << "Sorry, species or cluster " << rooting << " not found in your datafile/clusters file" << endl;
                    exit(EXIT_FAILURE);
                }
                rootCluster->invert(*outgroupCluster) ;
            }
        }
        clustersTree->addCluster( rootCluster );
    }
}
