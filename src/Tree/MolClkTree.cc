
#include "Tree/MolClkTree.h"

#include "Util/ParametersSet.h"

#include <iostream>
#include <assert.h>


MolClkTree::MolClkTree() : RootedTree(){
}

MolClkTree::MolClkTree(ParametersSet& parameters) : RootedTree(parameters){
}

MolClkTree::~MolClkTree(){
}


bool MolClkTree::clockTree(double& returnHeight, InferenceNode* startNode){
    //if no node provided start from the root
    if (!startNode){
        startNode = (InferenceNode*)root;
    }

    //if it is a leaf return true (recursive termination)
    //and the length to the parent, if this length is 0 set it to 0.05
    if ( startNode->isLeaf() ){
        returnHeight = startNode->getParentDistance();
        if (returnHeight<=0.0){
            returnHeight = 0.05;
            startNode->setParentDistance(0.05);
            return false;
        }
        return true;
    }

    bool res = true;
    bool subtreeWasRooted;
    double subtreeHeight;
    double maxSubtreeHeight = 0;
    vector<double> subtreeHeightsVector;

    for ( list<BasicNode*>::iterator iter = startNode->getChildrenList().begin();
          iter != startNode->getChildrenList().end(); ++iter ){
        subtreeWasRooted = clockTree( subtreeHeight,  (InferenceNode*)*iter );
        subtreeHeightsVector.push_back(subtreeHeight);
        if (maxSubtreeHeight < subtreeHeight){
            maxSubtreeHeight = subtreeHeight;
        }
        res = res && subtreeWasRooted;
    }
    int i = 0;
    list<BasicNode*>::iterator iter = startNode->getChildrenList().begin();
    while( iter != startNode->getChildrenList().end() ){
        rescale(maxSubtreeHeight/subtreeHeightsVector[i], (InferenceNode*)*iter);
        ++i;
        ++iter;
    }
    if ( startNode == (InferenceNode*)root ){
        returnHeight = maxSubtreeHeight;
    }
    else{
        returnHeight = startNode->getParentDistance() + maxSubtreeHeight;
    }
    if (!res){
        startNode->partialLikelihood = startNode->partialLikelihoodWork;
    }
    
    return res;
}

void MolClkTree::rescale(double scalingFactor){
    assert(scalingFactor>0.0);
    height = height * scalingFactor;
    rescale(scalingFactor, (InferenceNode*)root);
}

void MolClkTree::rescale(double scalingFactor, InferenceNode* startNode){
    if (startNode!=(InferenceNode*)root){
        startNode->setParentDistance(scalingFactor*(startNode->getParentDistance()));
    }
    if (!startNode->isLeaf()){
        startNode->partialLikelihood = startNode->partialLikelihoodWork;
        for ( list<BasicNode*>::iterator iter = startNode->getChildrenList().begin();
              iter != startNode->getChildrenList().end(); ++iter ){
            rescale(scalingFactor, (InferenceNode*)*iter);
        }
    }
}

void MolClkTree::constructFromString( string stringTree ){
    RootedTree::constructFromString( stringTree );
    bool clock = clockTree(height);
    //if stringTree was just a topology, avoid the warning
    if ( clock && (stringTree!=(toString(true)+';')) ){
        cerr << "WARNING: your tree has been modified to follow the molecular "
             << "clock assumption." << endl
             << "New length from root to leaves = " << getHeight() << endl;
    }
}

void MolClkTree::constructRandomly( double maxBranchLength ){
    RootedTree::constructRandomly( maxBranchLength/10.0 );
    clockTree(this->height);
}

void MolClkTree::setBranchVector( const vector<double> & branchVector ){
    RootedTree::setBranchVector( branchVector );
    //recompute height
    bool clock = clockTree( height, NULL );
    assert(clock);
}
