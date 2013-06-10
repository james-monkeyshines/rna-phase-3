#include "Tree/SearchTreeStarSearch.h"

#include "Util/ParametersSet.h"
#include "Util/Optimise.h"

#include "Models/Model.h"

#include "Tree/OptimizerNode.h"
#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

#include <float.h>
#include <iostream>
#include <iomanip>

using namespace std;


//register the tree to the tree factory with the name Star decomposition
SearchTreeStarSearch SearchTreeStarSearch::prototype("Star decomposition");


SearchTreeStarSearch::SearchTreeStarSearch( const string & registrationName ) :
SearchTreeBasic(registrationName){
}

SearchTreeStarSearch::SearchTreeStarSearch(ParametersSet& parameters) :
InferenceTree(),OptimizerTree(),SearchTreeBasic(parameters){
    root = new OptimizerNode(this);
    createIndex();
}

SearchTreeStarSearch::~SearchTreeStarSearch(){
}

OptimizerTree* SearchTreeStarSearch::clone(ParametersSet& parameters) const{
    return new SearchTreeStarSearch(parameters);
}

void SearchTreeStarSearch::recursiveBuilt(SequenceTable * ptable, ClustersTreeNode* clustersTreeNode, OptimizerNode* clusterAncestor){

    Singleton< randombox > & randBox = Singleton< randombox >::instance();
    clustersTreeNode->cluster->setAncestor( clusterAncestor );
    //get the cluster
    Cluster* cluster = clustersTreeNode->cluster->dup();
    for ( ClustersTreeNode::iterator iter = clustersTreeNode->begin();
          iter!=clustersTreeNode->end(); ++iter ){
        cluster->remove( *((*iter)->cluster) );
        OptimizerNode* child = new OptimizerNode( this );
        recursiveBuilt( ptable, *iter, child );
        add( child, clusterAncestor, randBox.ran()*.6 );
        if (child->getNumberChildren()>2){
            centers.push_back( child );
        }
    }
    //elements left in cluster are species assigned to this cluster and not
    //its children
    vector<unsigned int> species;
    cluster->getSpecies( species );
    for( vector< unsigned int >::iterator iter = species.begin();
         iter != species.end(); ++iter ){
        OptimizerNode* leaf = new OptimizerNode( ptable->species[*iter], this );
        add( leaf, clusterAncestor, randBox.ran()*.6 );
    }
    delete cluster;
}

bool SearchTreeStarSearch::initialisation( SequenceTable * ptable, Model * pmodel ){
    SearchTreeBasic::initialisation( ptable, pmodel );
    
    Singleton< randombox > & randBox = Singleton< randombox >::instance();

    centers.clear();
    if (clustersTree){
        recursiveBuilt( ptable, clustersTree->getRoot(), (OptimizerNode*)root );
    }
    else{
        for( unsigned int i = 0; i < ptable->getNumberSpecies(); ++i ){
            OptimizerNode* leaf = new OptimizerNode( ptable->species[i], this );
            add( (OptimizerNode*)leaf, (OptimizerNode*)root, randBox.ran()*.6);
        }
    }
    if (root->getNumberChildren()>3){
        centers.push_back( (OptimizerNode*)root );
    }
    //tree ready, create the index
    createIndex();
    //initialise outgroupNode
    findOutgroup();
    //resize lookup to accomodate the multifurcations
    lookupResize();
    Optimise::optimiseQuasiNewton( *this, *pmodel, optimizeModel, empiricalFreqs, optimizeModelPenalty ? 1e-9 : 1e-8, false, false );
    //return true since we have an unresolved initial state
    return true;
}

bool SearchTreeStarSearch::processNext(){

    if ( centers.empty() ){
       return false;
    }
    //count the number of optimizations to perform
    unsigned int nbTry = 0;
    for ( list< OptimizerNode* >::iterator iter = centers.begin(); iter!=centers.end(); ++iter ){
        //count the number of neighbours at each multifurcations
        unsigned int numberLink = (*iter)->getNumberChildren();
        if ( *iter != root ) ++numberLink;
        //a center should have at least 4 links
        assert(numberLink>3);
        if (numberLink>4){
            nbTry += numberLink*(numberLink-1) / 2;
        }
        //special case with four links since linking A,B or linking C,D is equivalent  
        else{
            nbTry += 3;
        }
    }

    cout << "Trying: " << nbTry << " combinations" << endl;
    double bestLikelihood = -DBL_MAX/3.0;
    OptimizerNode* bestNode1 = NULL;
    OptimizerNode* bestNode2 = NULL;
    list< OptimizerNode* >::iterator bestCenterIter = centers.begin();
    vector< double > oldBranches;
    getBranchVector( oldBranches );
    vector< double > bestBranches;
    vector< double > oldModelParameters;
    if (optimizeModel){
        pmodel->getAllParameters( oldModelParameters );
    }
    vector< double > bestModelParameters;
    vector< double > oldModelPenalty;
    double bestPenalty = 0.0;
    if (optimizeModelPenalty){
        bestPenalty = -DBL_MAX/3.0;
        pmodel->getAllPenaltyParameters( oldModelPenalty );
    }
    else{
        bestPenalty = pmodel->getLnPenalty();
    }
    vector< double > bestModelPenalty;
    unsigned int nb = 0;
    //each iteration of processNext, a multifurcation lose one child and a new node is added
    //initialise this new node
    OptimizerNode* newParent = new OptimizerNode(this);
    
    for ( list< OptimizerNode* >::iterator centersIter = centers.begin(); centersIter!=centers.end(); ++centersIter ){

        //find all the pairs for the given center
        /** possible pairs? *******************************************************************
        //                  4 links                           |       5links or more
        //root         | 3 pairs (between the three children  | all pairs of children possible
        //             |          which are not outgroup)     |       
        //other node   | 3 pairs (between the three children) | pairs with the parent also possible
        ************************************************************************************** */
        unsigned int numberChildren = (*centersIter)->getNumberChildren();
        vector< pair <OptimizerNode*, OptimizerNode*> > allowedLinks;
        allowedLinks.clear();
        //link within all children are allowed unless this is the root
        for ( list< BasicNode* >::iterator iterChild1 = ((*centersIter)->getChildrenList()).begin();
              iterChild1 !=  ((*centersIter)->getChildrenList()).end(); ++iterChild1 ){
            //if centre is the root and there are 4 children we must be sure not to pair the outgroup
            if ( (*iterChild1==outgroupNode) && (numberChildren == 4) ){
                //nothing, wait next iteration
                assert(*centersIter==root);
            }
            else{
                //if centre is not the root and more than 3 children, prepare pairs with the parent
                if ( (*centersIter != root) && ( numberChildren>3 ) ){
                    allowedLinks.push_back( pair <OptimizerNode*, OptimizerNode*>
                        ((OptimizerNode*)(*iterChild1), (OptimizerNode*)((*centersIter)->getParent())) );
                }
                list< BasicNode* >::iterator iterChild2 = iterChild1;
                ++iterChild2;
                while ( iterChild2 !=  ((*centersIter)->getChildrenList()).end() ){
                    //if centre is the root and there are 4 children we must be sure not to pair the outgroup
                    if ( (*iterChild2==outgroupNode) && (numberChildren == 4) ){
                        //nothing, wait next iteration
                        assert(*centersIter==root);
                    }
                    else{
                        allowedLinks.push_back( pair <OptimizerNode*, OptimizerNode*>
                        ((OptimizerNode*)(*iterChild1), (OptimizerNode*)(*iterChild2)) );
                    }
                    ++iterChild2;
                }
            }
        }
        //now try them...
        for ( vector< pair <OptimizerNode*, OptimizerNode*> >::iterator iter = allowedLinks.begin();
              iter != allowedLinks.end(); ++iter ){
            OptimizerNode* node1 = iter->first;
            OptimizerNode* node2 = iter->second;
            
            // prepare the tree
            link( node1, node2, newParent );
            
            // try optimization
            
            //cosmetic separation
            if (nb) cout << "--------" << endl;
            cout << ++nb << '.' << toString(true) << endl;
            Optimise::optimiseQuasiNewton( *this, *pmodel, optimizeModel, empiricalFreqs, optimizeModelPenalty ? 1e-9 : 1e-8, false, false  );
            double lik = loglikelihood();
            double pen = pmodel->getLnPenalty();
            if(optimizeModelPenalty){
                cout << "optimized value=" << lik + pen << "     (";
            }
            cout << "maxLikelihood=" << lik;
            if(optimizeModelPenalty){
                cout << ", maxPenalty=" << pen << ')';
            }
            else{
                if (bestPenalty){
                    cout << "      (penalty=" << bestPenalty << ')';
                }
            }
            cout << endl;
            if (lik+pen>bestLikelihood+bestPenalty){
                bestLikelihood = lik;
                bestPenalty = pen;
                bestNode1 = node1;
                bestNode2 = node2;
                bestCenterIter = centersIter;
                getBranchVector( bestBranches );
                if (optimizeModel){
                    pmodel->getAllParameters( bestModelParameters );
                }
                if (optimizeModelPenalty){
                    pmodel->getAllPenaltyParameters( bestModelPenalty );
                }
            }
                    
            /* ********** roll back for next try *********** */
            unlink( node1, node2, newParent );
            setBranchVector( oldBranches );
            if (optimizeModel){
                pmodel->setAllParameters( oldModelParameters );
            }
            if (optimizeModelPenalty){
                pmodel->setAllPenaltyParameters ( oldModelPenalty );
            }
            if (optimizeModel || optimizeModelPenalty){
                pmodel->validChange();
            }
        }
    }
    assert(nb==nbTry);
    link( bestNode1, bestNode2, newParent );
    setBranchVector( bestBranches );
    if (optimizeModel){
        pmodel->setAllParameters( bestModelParameters );
    }
    if (optimizeModelPenalty){
        pmodel->setAllPenaltyParameters ( bestModelPenalty );
    }
    if (optimizeModel || optimizeModelPenalty){
        pmodel->validChange();
    }
    //update centers
    if ( ( (*bestCenterIter)->getNumberChildren() == 2 ) || 
         ( ( (*bestCenterIter)->getNumberChildren() == 3 ) && ( *bestCenterIter == root ) ) ){
         centers.erase( bestCenterIter );
    }
    return true;
}


void SearchTreeStarSearch::link( OptimizerNode* node1, OptimizerNode* node2, OptimizerNode* newParent ){

    Singleton< randombox > & randBox = Singleton< randombox >::instance();
    OptimizerNode* oldParent = (OptimizerNode*)node1->getParent();
    
    //standard case
    if ( oldParent == node2->getParent() ){
        double branch1 = node1->getParentDistance();
        double branch2 = node2->getParentDistance();
        //link node1 and node2 with newParent
        remove(node1);
        remove(node2);
        double linkDist = randBox.ran()*min(branch1, branch2);
        add( node1, newParent, branch1 - linkDist );
        add( node2, newParent, branch2 - linkDist );
        add( newParent, oldParent, linkDist );
        //reroot if necessary
        if ( (node1==outgroupNode) || (node2==outgroupNode) ){
            installOutgroup();
        }
        else{
            createIndex();
        }
    }
    else{
        assert(oldParent->getParent() == node2 );
        double branch1 = node1->getParentDistance();
        double branch2 = oldParent->getParentDistance();
        remove(node1);
        remove(oldParent);
        double linkDist = randBox.ran()*min(branch1, branch2);
        add( node1, newParent, branch1 - linkDist );
        add( oldParent, newParent, linkDist );
        add( newParent, node2, branch2 - linkDist );
        createIndex();
    }
}

void SearchTreeStarSearch::unlink( OptimizerNode* node1, OptimizerNode* node2, OptimizerNode* newParent ){
    //standard case
    if ( newParent == node2->getParent() ){
        if (newParent==root){
            assert(root->getNumberChildren()==3);
            OptimizerNode* thirdChild = NULL;
            for ( list< BasicNode* >::iterator iterChild = (root->getChildrenList()).begin();
                iterChild !=  (root->getChildrenList()).end(); ++iterChild ){
                if ( (*iterChild!=node1) && (*iterChild!=node2) ){
                    assert(thirdChild==NULL);
                    thirdChild = (OptimizerNode*)(*iterChild);
                }
            }
            makeRoot( thirdChild );
            assert(thirdChild!=NULL);
            double branch1 = node1->getParentDistance() + newParent->getParentDistance();
            double branch2 = node2->getParentDistance() + newParent->getParentDistance();
            //unlink node1, node2 and newParent from root==thirdChild
            remove(newParent);
            remove(node1);
            remove(node2);
            add( node1, thirdChild, branch1);
            add( node2, thirdChild, branch2);
            createIndex();
        }
        else{
            OptimizerNode* oldParent = (OptimizerNode*)newParent->getParent();
            createIndex();            
            double branch1 = node1->getParentDistance() + newParent->getParentDistance();
            double branch2 = node2->getParentDistance() + newParent->getParentDistance();
            remove(node1);
            remove(node2);
            remove(newParent);
            add( node1, oldParent, branch1);
            add( node2, oldParent, branch2);
            createIndex();
        }
    }
    else{
        assert(newParent->getParent() == node2 );
        assert(newParent->getNumberChildren()==2);
        OptimizerNode* secondChild = NULL;
        for ( list< BasicNode* >::iterator iterChild = (newParent->getChildrenList()).begin();
            iterChild !=  (newParent->getChildrenList()).end(); ++iterChild ){
            if (*iterChild!=node1){
                assert(secondChild==NULL);
                secondChild = (OptimizerNode*)(*iterChild);
            }
        }
        assert(secondChild!=NULL);
        double branch1 = node1->getParentDistance() + secondChild->getParentDistance();
        double branch2 = newParent->getParentDistance() + secondChild->getParentDistance();
        remove(newParent);
        remove(node1);
        remove(secondChild);
        add( node1, secondChild, branch1 );
        add( secondChild, node2, branch2 );
        createIndex();
    }
}


void SearchTreeStarSearch::printResults( ostream& outputStream ){
    outputStream << "Model found:" << endl;
    pmodel->printParameters( outputStream );
    
    outputStream << "Tree found:" << endl;
    outputStream << toString() << ';' << endl;

    outputStream.setf(ios::fixed);
    outputStream << setprecision(4);
    double lik = loglikelihood();
    double pen = 0.0;
    if (optimizeModelPenalty){
        pen = pmodel->getLnPenalty();
    }
    if(optimizeModelPenalty){
        outputStream << "optimized value=" << lik + pen << "     (";
    }
    outputStream << "maxLikelihood=" << lik;
    if(optimizeModelPenalty){
        outputStream << ", maxPenalty=" << pen << ')';
    }
    else{
        if (pen){
            outputStream << "      (penalty=" << pen << ')';
        }
    }
    outputStream << endl;
}

void SearchTreeStarSearch::printTree( ostream& outputStream ){
    // Save a tree to a file on its own
    outputStream << toString(false) << ';' << endl;
}
