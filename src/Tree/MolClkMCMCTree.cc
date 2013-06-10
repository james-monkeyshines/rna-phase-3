#include "Tree/MolClkMCMCTree.h"

#include "Util/ParametersSet.h"
#include "Models/PerturbatorGaussParameter.h"

#include <iostream>
#include <iomanip>
#include <assert.h>

#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

//register the tree to the tree factory with the name Rooted MCMC tree with molecular clock
MolClkMCMCTree MolClkMCMCTree::prototype("Rooted MCMC tree with molecular clock");

void MolClkMCMCTree::constructFromString( string stringTree ){
    MolClkTree::constructFromString( stringTree );
    checkBinary();
}


MolClkMCMCTree::MolClkMCMCTree( const string & registrationName ) :
InferenceTree(), MCMCTreeBasic(registrationName), MolClkTree(){
}

MolClkMCMCTree::MolClkMCMCTree(ParametersSet& parameters) :
InferenceTree(), MCMCTreeBasic(parameters), MolClkTree(parameters){
}

MCMCTree* MolClkMCMCTree::clone(ParametersSet& parameters) const{
    return new MolClkMCMCTree(parameters);
}

MolClkMCMCTree::~MolClkMCMCTree(){
    if (treeHeightPerturbator){
        delete treeHeightPerturbator;
    }
    if (branchLengthPerturbator){
        delete branchLengthPerturbator;
    }
}

void MolClkMCMCTree::getCloseFarChildren( InferenceNode* internalNode,
                         InferenceNode*& closerChild,
                         InferenceNode*& farerChild ){
    list<BasicNode*>::const_iterator iter =
            (internalNode->getChildrenList()).begin();
    closerChild = (InferenceNode*)(*iter);
    ++iter;
    if ( closerChild->getParentDistance() > (*iter)->getParentDistance() ){
        farerChild = closerChild;
        closerChild = (InferenceNode*)(*iter);
    }
    else{
        farerChild = (InferenceNode*)(*iter);
    }
}


double MolClkMCMCTree::changeBranchLength(){
    assert( lastLocalChange == INVALID_LOCAL );
    assert( lastGlobalChange == INVALID_GLOBAL );

    //choose a random node (not a leaf)
    int modifiedBranchId = (int)(getNumberBranches() * randBox.ran());
    internalNode = (InferenceNode*)getBranchOrig(modifiedBranchId);
    //if the node is the root expand or shrink the length from root to leaves
    if (internalNode==root){
        //use treeLengthPerturbator to modify treeLengthVariable
        double actualTreeHeight = getHeight();
        double newTreeHeight = actualTreeHeight;
        treeHeightPerturbator->perturbValue( &newTreeHeight );
        //rescale the tree
        rescale(newTreeHeight/actualTreeHeight);
        lastLocalChange = TREE_HEIGHT;        
        return 0.0;
    }
    //if it is not the root the node is moved between its children and parent
    //retrieve its child (and find the closest)
    InferenceNode* closerChild;
    InferenceNode* farerChild;
    getCloseFarChildren( internalNode, closerChild, farerChild );
    //store the distance to the parent and to the children
    double oldParentDist = internalNode->getParentDistance();
    double oldCloserChildDist = closerChild->getParentDistance();
    double oldFarerChildDist = farerChild->getParentDistance();
    double delta = oldParentDist + oldCloserChildDist;
    double var = oldParentDist/delta;
    branchLengthPerturbator->perturbValue( &var );
    //var is in [0;+1]
    double newParentDist = var*delta;
    /**
     * if rebound against the father (low bound 0.0) then nothing
     * if rebound against the child (high bound 1.0) initiate a NNI
     */
    internalNode->setParentDistance(newParentDist);
    closerChild->setParentDistance(delta-newParentDist);
    farerChild->setParentDistance(oldFarerChildDist+oldParentDist-newParentDist);
    if ( branchLengthPerturbator->getHighShock() && (treePriority!=0.0) ){
        //the node must not be a leaf
        bool NNIpossible = !closerChild->isLeaf();
        //and must not define a cluster too
        if ( clustersTree && NNIpossible ){
            if (clustersTree->findByAncestor(closerChild) != NULL){
                NNIpossible = false;
            }
        }
        if ( NNIpossible ){
            //the rebound initiates a NNI
            InferenceNode* swappedChild;
            if ( randBox.ran() < .5 ){
                swappedChild = (InferenceNode*)closerChild->getChild(0);
            }
            else{
                swappedChild = (InferenceNode*)closerChild->getChild(1);
            }
            pertSave.clear();
            pertSave.push_back(swappedChild);
            pertSave.push_back(farerChild);
            NNI(swappedChild,farerChild);
            lastLocalChange = NNI_LOCAL;
            ++numberLocalNNI;
            return 0.0;
        }
    }
    internalNode->invalidateRecursively(-1);
    pertSave.clear();
    pertSave.push_back(closerChild);
    pertSave.push_back(farerChild);
    lastLocalChange = LENGTH;
    return 0.0;
}


bool MolClkMCMCTree::validateBranchLength( bool validation ){
    assert(internalNode);
    assert( lastLocalChange != INVALID_LOCAL );
    double oldParentFrac;

    //after a NNI induced by a negative branch length
    switch(lastLocalChange){
        case TREE_HEIGHT:
          double oldTreeHeight;
          treeHeightPerturbator->validatePerturbationValue( validation, &oldTreeHeight );
          if (validation){
              saveNodes(-1);
          }
          else{
              rescale(oldTreeHeight/getHeight());
              retrieveNodes(-1);
          }
        break;
        case LENGTH:
          assert(pertSave.size()==2);
          branchLengthPerturbator->validatePerturbationValue( validation, &oldParentFrac );
          if (!validation){
              //store temporarily the refused parent dist in distDelta
              double distDelta = internalNode->getParentDistance();
              //compute the CONSTANT length delta = distance(father,closer_child)
              double delta = pertSave[0]->getParentDistance()+distDelta;
              internalNode->setParentDistance(oldParentFrac * delta);
              //store refusedValue - oldValue in distDelta
              distDelta -= oldParentFrac * delta;
              pertSave[0]->setParentDistance( (1.0-oldParentFrac) * delta );
              pertSave[1]->setParentDistance(pertSave[1]->getParentDistance()+distDelta);
              internalNode->retrieveRecursively(-1);
          }
          else{
              internalNode->saveRecursively(-1);
          }
        break;
        case NNI_LOCAL:
          assert(pertSave.size()==2);
          branchLengthPerturbator->validatePerturbationValue( validation, &oldParentFrac );
          if (validation){
              ++acceptedLocalNNI;
              //save from the closer child
              ((InferenceNode*)pertSave[1]->getParent())->saveRecursively(-1);
          }
          else{
              //revert the NNI
              NNI(pertSave[1],pertSave[0]);
              //store temporarily the refused parent dist in distDelta
              double distDelta = internalNode->getParentDistance();
              //closerChild = pertSave[0]->parent()
              pertSave[0] = (InferenceNode*)pertSave[0]->getParent();
              //compute the CONSTANT length delta = distance(father,closer_child)
              double delta = pertSave[0]->getParentDistance()+distDelta;
              internalNode->setParentDistance(oldParentFrac * delta);
              //store refusedValue - oldValue in distDelta
              distDelta -= oldParentFrac * delta;
              pertSave[0]->setParentDistance( (1.0-oldParentFrac) * delta );
              pertSave[1]->setParentDistance(pertSave[1]->getParentDistance()+distDelta);
              internalNode->retrieveRecursively(-1);
          }
        break;
        default:
          assert(0);
    }  //end switch
    lastLocalChange = INVALID_LOCAL;
    internalNode = NULL;
    return true;
}




double MolClkMCMCTree::globalTopologyChange(){
    assert( lastLocalChange == INVALID_LOCAL );
    assert( lastGlobalChange == INVALID_GLOBAL );

    double lnHastingRatio;
    
    /* Sweeping */
    if ( randBox.ran() < 0.5 ) {
        double branchDist;
        double dist;
        pertSave.resize(2);
        vector< pair<InferenceNode*,double> > attachPointVector;
        //choose a node (not the root nor its children)
        do{
            assert( root->getNumberChildren() == 2 );
            int sweepedBranchId = (int)((getNumberBranches()-2) * randBox.ran()) + 2;
            pertSave[0] = (InferenceNode*)getBranchDest(sweepedBranchId);


            //retrieve the distance from the node to the root
            //or its parent cluster
            ClustersTreeNode* orig = NULL;
            double distanceMax = getClusterDistance( pertSave[0], orig );

            //find the attach point (the new brother for node)
            attachPointVector.clear();
            branchDist = findBranches( attachPointVector, distanceMax, pertSave[0], false, orig );
//            if (attachPointVector.size()<=2){
//                cout << "impossible sweep:" << pertSave[0]->toString() << endl;
//            }
        } while (attachPointVector.size()<=2);


        InferenceNode* parent = (InferenceNode*)pertSave[0]->getParent();

        assert( parent!=NULL );
        assert( parent!=root );
        //store the distance from the parent to the gp
        oldDist = parent->getParentDistance();
        //store the old brother
        pertSave[1] = (InferenceNode*)parent->getChild(0);
        if ( pertSave[1] == pertSave[0] ){
            pertSave[1] = (InferenceNode*)parent->getChild(1);
        }

        //the current branch is in attachPointVector and we have to skip it
        //(the current branch is splitted in two by parent actually)
        branchDist -= min(pertSave[0]->getParentDistance(),pertSave[1]->getParentDistance());
        branchDist -= parent->getParentDistance();
        assert(branchDist>0.0);
        dist = randBox.ran() * branchDist;
        vector< pair<InferenceNode*,double> >::iterator iter = attachPointVector.begin();
        while ( (dist > iter->second) || (iter->first==pertSave[1]) || (iter->first==parent) ){
            if ( (iter->first!=pertSave[1]) && (iter->first!=parent) ){
                dist -= iter->second;
            }
            ++iter;
            assert(iter!=attachPointVector.end());
        }

        //perform the sweep
        sweeping( pertSave[0], iter->first, dist );
        ++numberSweepingPerturbation;
        lastGlobalChange = SWEEP_PROP;
        lnHastingRatio = 0.0;
    }
    /* Nearest Neighbor Interchange */
    else{
        //there are two versions of this algorithm, one standard and
        //one for trees with tip date and clusters
        if (!clustersTree){
            InferenceNode* closerChild = NULL;
            InferenceNode* farerChild;
            //choose a random node with at least one internal child
            //and cycle until a possible one is found
            while (!closerChild){
                int modifiedBranchId = (int)(getNumberBranches() * randBox.ran());
                //retrieve its child (and find the closest)
                getCloseFarChildren((InferenceNode*)getBranchOrig(modifiedBranchId), closerChild, farerChild);
                //NNi not allowed if the branch is a terminal branch
                if ( closerChild->isLeaf() ){
                    closerChild = NULL;
                }
            }
            InferenceNode* swappedChild;
            if ( randBox.ran() < .5 ){
                swappedChild = (InferenceNode*)closerChild->getChild(0);
            }
            else{
                swappedChild = (InferenceNode*)closerChild->getChild(1);
            }
            pertSave.clear();
            pertSave.push_back(swappedChild);
            pertSave.push_back(farerChild);
            NNI(swappedChild,farerChild);
            ++numberNNIPerturbation;
            lastGlobalChange = NNI_PROP;
            lnHastingRatio = 0.0;
        }
        //slide NNI
        else{
            pertSave.resize(2);
            double grandParentDist = 0.0;
            double child1Dist = 0.0;
            double child2Dist = 0.0;
            InferenceNode* parent = NULL;
            InferenceNode* grandParent = NULL;
            do {
               //choose a branch to slide (not linked to the root)
                int slidedBranchId = (int)((getNumberBranches()-2) * randBox.ran()) + 2;
                //recover parent and child of the slided branch
                parent = (InferenceNode*)getBranchOrig(slidedBranchId);
                pertSave[0]=(InferenceNode*)getBranchDest(slidedBranchId);
                assert(parent!=root);
                //retrieve the parent of parent and its second child
                grandParent = (InferenceNode*)parent->getParent();
                assert(grandParent);
                pertSave[1] = (InferenceNode*)parent->getChild(0);
                if ( pertSave[1] == pertSave[0] ){
                    pertSave[1] = (InferenceNode*)parent->getChild(1);
                }
                //the branch can slide above grandparent if it is not the root and if we are not
                //in a cluster
                if ( (grandParent!=root) && (!clustersTree || !clustersTree->findByAncestor(parent)) ){
                    grandParentDist = grandParent->getParentDistance();
                }
                //the branch can slide below secondChild if secondChild does not define a cluster and...
                if ( !pertSave[1]->isLeaf() && (!clustersTree || !clustersTree->findByAncestor(pertSave[1])) ){
                    double distChild = pertSave[0]->getParentDistance()-pertSave[1]->getParentDistance();
                    //... if secondChild is higher than pertSave[0], the child of the slided branch
                    if (distChild>0.0){
                        child1Dist = min(distChild, pertSave[1]->getChildDistance(0) );
                        child2Dist = min(distChild, pertSave[1]->getChildDistance(1) );
                    }
                }
            } while ( !grandParentDist && !child1Dist && !child2Dist ); //choose another branch if no slide possible
            //store the old distance
            oldDist = pertSave[1]->getParentDistance();
            //store temporarily the sum of the possible branches in lnHastingRatio
            lnHastingRatio = grandParentDist+child1Dist+child2Dist;
            //and decide where to slide
            double ran = (lnHastingRatio) * randBox.ran();
            ran -= grandParentDist;
            double denom = 0.0;
            //slide to the parent
            if (ran < 0.0){
                InferenceNode* newGrandParent = (InferenceNode*)grandParent->getParent();
                assert(newGrandParent);
                if ( (newGrandParent!=root) && (!clustersTree || !clustersTree->findByAncestor(grandParent)) ){
                    denom += grandParent->getParentDistance();
                }
                denom += pertSave[1]->getParentDistance() + grandParent->getChildDistance(0) + grandParent->getChildDistance(1);
                slideNNI( pertSave[0], grandParent, -ran );
            }
            else{
                InferenceNode* chosenChild = NULL;
                ran -= child1Dist;
                if (ran<0.0){
                    chosenChild = pertSave[1]->getChild(0);
                }
                else{
                    ran -= child2Dist;
                    chosenChild = pertSave[1]->getChild(1);
                }
                assert(ran<0.0);
                denom = parent->getParentDistance() + pertSave[1]->getParentDistance();
                if ( !chosenChild->isLeaf() && (!clustersTree || !clustersTree->findByAncestor(chosenChild)) ){
                    double distChild = pertSave[0]->getParentDistance()-pertSave[1]->getParentDistance()-chosenChild->getParentDistance();
                    if (distChild>0.0){
                        denom += min(distChild, chosenChild->getChildDistance(0) );
                        denom += min(distChild, chosenChild->getChildDistance(1) );
                    }
                }
                slideNNI( pertSave[0], chosenChild, chosenChild->getParentDistance()+ran );
            }
            ++numberNNIPerturbation;
            lastGlobalChange = NNI_SLIDE;
            lnHastingRatio = log(lnHastingRatio/denom);
        }
    }
    return lnHastingRatio;
}



bool MolClkMCMCTree::validateTopologyChange( bool validation ){
    assert( lastGlobalChange != INVALID_GLOBAL );
    InferenceNode* temp;
    //modified and modifiedStructure have been invalidated,
    //whatever the value of validation
    switch( lastGlobalChange ){
      case SWEEP_PROP:
        assert(pertSave.size()==2);
        if (validation){
            ++acceptedSweepingPerturbation;
            ((InferenceNode*)(pertSave[1]->getParent()))->saveRecursively(-1);
            ((InferenceNode*)(pertSave[0]->getParent()))->saveRecursively(-1);
        }
        else{
            InferenceNode* missedGrandParent = (InferenceNode*)(pertSave[0]->getParent()->getParent());
            sweeping( pertSave[0], pertSave[1], oldDist );
            missedGrandParent->retrieveRecursively(-1);
            ((InferenceNode*)(pertSave[0]->getParent()))->retrieveRecursively(-1);
        }
      break;
      case NNI_PROP:
        assert(pertSave.size()==2);
        assert(!clustersTree);
        //standard NNI
        if (validation){
            ++acceptedNNIPerturbation;
            ((InferenceNode*)(pertSave[1]->getParent()))->saveRecursively(-1);
        }
        else{
            NNI(pertSave[1],pertSave[0]);
                ((InferenceNode*)(pertSave[0]->getParent()))->retrieveRecursively(-1);
        }
      break;
      case NNI_SLIDE:
        assert(pertSave.size()==2);
        assert(clustersTree);
        //slideNNI
        if (validation){
            ++acceptedNNIPerturbation;
            temp = (InferenceNode*)(pertSave[0]->getParent());
            if (pertSave[1]->getParent()->getParent()==temp){
                //save from pertSave[1]->parent (the previous grandParent)
                temp = (InferenceNode*)(pertSave[1]->getParent());
            }
            else{
                if (temp->getParent()==pertSave[1]){
                    //save from temp (the parent of the slided branch)
                }
                else{
                    assert(0);
                }
            }
            temp->saveRecursively(-1);
        }
        else{
            temp = (InferenceNode*)pertSave[0]->getParent();
            if (pertSave[1]->getParent()->getParent()==temp){
                //retrieve from pertSave[0]->getParent once it is back in place
            }
            else{
                if (temp->getParent()==pertSave[1]){
                    //retrieve from pertSave[1] once pertSave[0] is back in place
                    temp = pertSave[1];
                }
                else{
                    assert(0);
                }
            }
            slideNNI(pertSave[0],pertSave[1],oldDist);
            temp->retrieveRecursively(-1);
        }
      break;
      default:
        assert(0);
      break;
    }
    lastGlobalChange = INVALID_GLOBAL;
    return true;
}


void MolClkMCMCTree::stopBurn(){
    MCMCTreeBasic::stopBurn();
    numberSweepingPerturbation = 0;
    acceptedSweepingPerturbation = 0;
    numberNNIPerturbation = 0;
    acceptedNNIPerturbation = 0;
    numberLocalNNI = 0;
    acceptedLocalNNI = 0;
    branchLengthPerturbator->stopBurn();
    treeHeightPerturbator->stopBurn();
}

void MolClkMCMCTree::initialiseMCMC( ParametersSet & parameters ){
    numberSweepingPerturbation = 0;
    acceptedSweepingPerturbation = 0;
    numberNNIPerturbation = 0;
    acceptedNNIPerturbation = 0;

    numberLocalNNI = 0;
    acceptedLocalNNI = 0;

    lastGlobalChange = INVALID_GLOBAL;
    lastLocalChange = INVALID_LOCAL;
    internalNode = NULL;

    if (parameters.findParameter("Branch lengths, prior")){
        cerr << "Priors on branch lengths not implemented with rooted trees"
             << "defaulting to a uniform prior" << endl
             << "Use the prior on the height of the tree" << endl;
    }
    //It is a node sliding, the parameter we modify is in fact the fraction
    //D(node->parent,node)/D(node->parent,node->closerchild)
    //bounded in [0.0,1.0]
    PriorField branchPriors("uniform(0.0,1.0)");
    branchLengthPerturbator = new PerturbatorHelper( hyperPerturbator, *this, UpdateMessage::BRANCH,
                                               "Branch lengths", branchPriors,
                                               parameters, "Branch lengths",
                                               .2, 0.0, 0.9, "initial step" );

    PriorField treeHeightPrior(parameters.stringParameter("Tree height, prior"));
    treeHeightPerturbator = new PerturbatorHelper( hyperPerturbator, *this, UpdateMessage::BRANCH,
                                             "Tree height", treeHeightPrior,
                                             parameters, "Tree height",
                                             .05, 0.0, 1.0, "initial step" );
    double min, max;
    treeHeightPerturbator->getExtrema( min, max );
    double height = getHeight();
    if ( (height<min) || (height>max) ){
        double h = treeHeightPerturbator->getMean();
        if (h>3.0){
            h = min + .1;
        }
        rescale(h/height);
        cerr << "WARNING: branch lengths modified to fit the prior on the height of the tree" << endl;
    }
    MCMCTreeBasic::initialiseMCMC(parameters);
}

void MolClkMCMCTree::printPerturbationParameters(ostream& outputStream){
    outputStream << setprecision(5);
    outputStream << "Sweeping acceptance rate                    =  "
                 << (double)acceptedSweepingPerturbation*100.0 / (double)numberSweepingPerturbation
                 << '%' << endl;
    outputStream << "Number of Sweeping proposals                =  "
                 << numberSweepingPerturbation << endl;
    outputStream << "Global NNI acceptance rate             =  "
                 << (double)acceptedNNIPerturbation*100.0 / (double)numberNNIPerturbation
                 << '%' << endl;
    outputStream << "Number of NNI proposals                =  "
                 << numberNNIPerturbation << endl;

    outputStream << setprecision(2);
    outputStream << "Local NNI (induced by a branch length proposal) acceptance rate =  "
                 << (double)acceptedLocalNNI*100.0 / (double)numberLocalNNI
                 << '%' << endl;
    outputStream << "Number of local NNI proposals                                   =  "
                 << numberLocalNNI << endl;

    branchLengthPerturbator->printPerturbationParameters( outputStream );
    treeHeightPerturbator->printPerturbationParameters( outputStream );

    MCMCTreeBasic::printPerturbationParameters( outputStream );
}


void MolClkMCMCTree::NNI( InferenceNode* swappedChild, InferenceNode* farerChild ){
    assert(swappedChild->getParent()->getParent() == farerChild->getParent());
    assert(swappedChild->getParent()->getParentDistance() < farerChild->getParentDistance());
    /**
     *       ------------------- farerChild
     *   * <                   ---------- swappedChild
     *       --- closerChild <
     *                        --------------------- X
     */
    double farerDist = farerChild->getParentDistance();
    double closerDist = swappedChild->getParent()->getParentDistance();
    double swappedDist = swappedChild->getParentDistance();
    //for both node the distance to the root should be the
    //same at the end of the process
    farerChild->setParentDistance(farerDist-closerDist);
    swappedChild->setParentDistance(swappedDist+closerDist);
    swap(farerChild, swappedChild);
}

void MolClkMCMCTree::slideNNI( InferenceNode* slideChild, InferenceNode* destChild, double dist ){

    InferenceNode* parent = (InferenceNode*)slideChild->getParent();
    InferenceNode* grandParent = (InferenceNode*)parent->getParent();
    InferenceNode* invalidate = NULL;
    InferenceNode* otherChild = (InferenceNode*)parent->getChild(0);
    assert(dist>0.0);
    double newDist;
    if (otherChild==slideChild){
        otherChild = (InferenceNode*)parent->getChild(1);
    }

    /*                                                                  otherChild
     *  grandgrandParent--->grandParent(possible destChild)---->parent<
     *                                                                  slideChild
     */
    //if we are sliding  to the grandParent
    if ( destChild == grandParent ){
        assert(dist<grandParent->getParentDistance());
        invalidate = grandParent;
        newDist = slideChild->getParentDistance() + parent->getParentDistance() + dist;
        //if grandparent was the ancestor of a cluster the parent will take this role
        if ( clustersTree ){
            ClustersTreeNode* clustersTreeNode = clustersTree->findByAncestor(grandParent);
            if (clustersTreeNode ){
                clustersTreeNode->cluster->setAncestor( parent );
            }
        }
    }
    /*              otherChild ---> destChild
     * ---->parent<
     *              slideChild
     */
    //sliding to a branch from otherChild
    else{
        assert(destChild->getParent() == otherChild);
        assert(dist<destChild->getParentDistance());
        invalidate = parent;
        newDist = slideChild->getParentDistance() - otherChild->getParentDistance() - destChild->getParentDistance() + dist;
        //if parent was the ancestor of a cluster, otherChild will take this role
        if ( clustersTree ){
            ClustersTreeNode* clustersTreeNode = clustersTree->findByAncestor(parent);
            if (clustersTreeNode ){
                clustersTreeNode->cluster->setAncestor( otherChild );
            }
        }
    }
    assert(newDist>0.0);


    //detach the branch from the tree.
    double branchDist = otherChild->getParentDistance() + parent->getParentDistance();
    parent->removeChild( otherChild );
    grandParent->removeChild( parent );
    grandParent->addChild( otherChild, branchDist );

    //reinsert parent
    branchDist = destChild->getParentDistance();
    assert( branchDist > dist );
    InferenceNode* newGrandParent = (InferenceNode*)destChild->getParent();
    newGrandParent->removeChild( destChild );
    newGrandParent->addChild( parent, branchDist - dist );
    parent->addChild( destChild, dist );
    slideChild->setParentDistance( newDist );

    invalidate->updateNodesCountRecursively();
    invalidate->updateIdenticalRecursively();
    invalidate->invalidateRecursively( -1 );
    createIndex();
}


void MolClkMCMCTree::sweeping( InferenceNode* node, InferenceNode* newBrother, double dist ){

    assert(dist>0.0);
    InferenceNode* parent = (InferenceNode*)(node->getParent());
    assert( parent );
    InferenceNode* oldGrandParent = (InferenceNode*)(parent->getParent());
    assert( oldGrandParent );
    InferenceNode* newGrandParent = (InferenceNode*)(newBrother->getParent());
    assert( newGrandParent );

    //save the old distance to the root
    double rootDistance = getRootDistance(node);
    //detach the node (and its parent) from the tree.
    InferenceNode* otherChild = (InferenceNode*)parent->getChild(0);
    if (otherChild==node){
        otherChild = (InferenceNode*)parent->getChild(1);
    }

    double branchDist = otherChild->getParentDistance() +
                       parent->getParentDistance();
    parent->removeChild( otherChild );
    oldGrandParent->removeChild( parent );
    oldGrandParent->addChild(otherChild, branchDist);

    //put parent at the attachPoint (between newGrandParent and newBrother)
    branchDist = newBrother->getParentDistance() - dist;
    newGrandParent->removeChild(newBrother);
    newGrandParent->addChild( parent, dist );
    parent->addChild( newBrother, branchDist );

    //resize the branch between node and parent to a correct value
    branchDist = rootDistance-getRootDistance(parent);
    assert( branchDist > 0.0 );
    node->setParentDistance( branchDist );

    if (clustersTree){
        //if the parent of the swapped branch was a cluster ancestor
        //newBrother becomes cluster ancestor
        ClustersTreeNode * clustersTreeNode = clustersTree->findByAncestor( parent );
        if ( clustersTreeNode ){
            clustersTreeNode->cluster->setAncestor( otherChild );
        }
        //if the new brother was a cluster ancestor
        //and if node came from this cluster
        clustersTreeNode = clustersTree->findByAncestor( newBrother );

        //we use otherChild to check whether node was a descendant of new
        //brother before it was moved
        if ( clustersTreeNode && otherChild->isDescendant(newBrother) ){
            clustersTreeNode->cluster->setAncestor( parent );
        }
    }
    //perform all the invalidation
    InferenceNode* ancestor = (InferenceNode*)getMostRecentAncestor(parent, oldGrandParent);
    parent->updateNodesCountRecursively(ancestor);
    parent->updateIdenticalRecursively(ancestor);
    parent->invalidateRecursively( -1, ancestor);
    oldGrandParent->updateNodesCountRecursively();
    oldGrandParent->updateIdenticalRecursively();
    oldGrandParent->invalidateRecursively( -1 );
    createIndex();
}


void MolClkMCMCTree::getAllPerturbationParameters( vector<double>& params ) const{
    MCMCTreeBasic::getAllPerturbationParameters(params);
    params.push_back((double)numberSweepingPerturbation);
    params.push_back((double)acceptedSweepingPerturbation);
    params.push_back((double)numberNNIPerturbation);
    params.push_back((double)acceptedNNIPerturbation);
    params.push_back((double)numberLocalNNI);
    params.push_back((double)acceptedLocalNNI);
    vector<double> p;
    treeHeightPerturbator->getAllPerturbationParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
    branchLengthPerturbator->getAllPerturbationParameters( p );
    params.insert( params.end(), p.begin(), p.end() );    
    assert( params.size() == MolClkMCMCTree::getNumberPerturbationParameters() );
}

void MolClkMCMCTree::setAllPerturbationParameters( const vector<double>& params){
    assert( params.size() == getNumberPerturbationParameters() );
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp = MCMCTreeBasic::getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    MCMCTreeBasic::setAllPerturbationParameters(p);
    numberSweepingPerturbation = (unsigned int)(*iter);
    ++iter;
    acceptedSweepingPerturbation = (unsigned int)(*iter);
    ++iter;
    numberNNIPerturbation = (unsigned int)(*iter);
    ++iter;
    acceptedNNIPerturbation = (unsigned int)(*iter);
    ++iter;
    numberLocalNNI = (unsigned int)(*iter);
    ++iter;
    acceptedLocalNNI = (unsigned int)(*iter);
    ++iter;
    nbp = treeHeightPerturbator->getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    treeHeightPerturbator->setAllPerturbationParameters ( p );
    nbp = branchLengthPerturbator->getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    branchLengthPerturbator->setAllPerturbationParameters ( p );
}

unsigned int MolClkMCMCTree::getNumberPerturbationParameters() const{
    unsigned int numberPerturbationParameters = 6 + MCMCTreeBasic::getNumberPerturbationParameters();
    numberPerturbationParameters += treeHeightPerturbator->getNumberPerturbationParameters();
    numberPerturbationParameters += branchLengthPerturbator->getNumberPerturbationParameters();
    return numberPerturbationParameters;
}

void MolClkMCMCTree::getAllPriorParameters( vector<double>& params ) const{
    MCMCTreeBasic::getAllPriorParameters(params);
    vector<double> p;
    treeHeightPerturbator->getAllPriorParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
    branchLengthPerturbator->getAllPriorParameters( p );
    assert( p.size() == 0 );
    //params.insert( params.end(), p.begin(), p.end() ); no prior at the moment
    assert( params.size() == getNumberPriorParameters() );
}

void MolClkMCMCTree::setAllPriorParameters( const vector<double>& params){
    assert( params.size() == getNumberPriorParameters() );
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp = MCMCTreeBasic::getNumberPriorParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    MCMCTreeBasic::setAllPriorParameters(p);
    nbp = treeHeightPerturbator->getNumberPriorParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    treeHeightPerturbator->setAllPriorParameters( p );
    assert( iter == params.end() );
}

unsigned int MolClkMCMCTree::getNumberPriorParameters() const{
    unsigned int numberPriorParameters = MCMCTreeBasic::getNumberPriorParameters();
    numberPriorParameters += treeHeightPerturbator->getNumberPriorParameters();
    return numberPriorParameters;
}

double MolClkMCMCTree::getLnPrior() const{
    double lnPrior = MCMCTreeBasic::getLnPrior();
    return treeHeightPerturbator->getLnPriorValue(getHeight()) + lnPrior;
}


