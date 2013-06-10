#include "Tree/HeterogeneousMCMCTree.h"

#include "Util/ParametersSet.h"
#include "Util/FileParser.h"

#include <iostream>
#include <iomanip>

#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

//register the tree to the tree factory with the name Heterogeneous MCMC tree
HeterogeneousMCMCTree HeterogeneousMCMCTree::prototype("Heterogeneous MCMC tree");


HeterogeneousMCMCTree::HeterogeneousMCMCTree( const string & registrationName ) :
InferenceTree(), MCMCTreeBasic(registrationName), RootedTree(){
}

HeterogeneousMCMCTree::HeterogeneousMCMCTree(ParametersSet& parameters) :
InferenceTree(), MCMCTreeBasic(parameters), RootedTree(parameters){
}

HeterogeneousMCMCTree::~HeterogeneousMCMCTree(){
    if (swapPerturbator){
        delete swapPerturbator;
    }
    if (branchPerturbator){
        delete branchPerturbator;
    }
}

MCMCTree* HeterogeneousMCMCTree::clone(ParametersSet& parameters) const{
    return new HeterogeneousMCMCTree(parameters);
}

void HeterogeneousMCMCTree::initSampling(ParametersSet& parameters, bool overwrite){
    MCMCTreeBasic::initSampling(parameters, overwrite );
    string outputFileBaseName( parameters.stringParameter( "Output file" ) );

    FileParser::confirmOpenFile( branchModelsFile,
                          outputFileBaseName + ".bm", overwrite );
}

void HeterogeneousMCMCTree::sample(){
    MCMCTreeBasic::sample();
    vector < BasicNode * >::const_iterator iter = nodeRefVector.begin();
    ++iter;
    //branch models sampling
    while ( iter != nodeRefVector.end() ){
        branchModelsFile << setw(4) << modelMap[((InferenceNode*)(*iter))->getModel()] << ' ';
        ++iter;
    }
    branchModelsFile << endl;
}

double HeterogeneousMCMCTree::changeBranchLength(){
    assert( lastLocalChange == INVALID_LOCAL );
    assert( lastGlobalChange == INVALID_GLOBAL );
    modifiedBranchId = ( unsigned int )( getNumberBranches() * randBox.ran() );
    double newValue = getBranchLength( modifiedBranchId );
    double lnHR =
        branchPerturbator->perturbValue( &newValue );
    assert( lnHR == 0.0 );
    setBranchLength( modifiedBranchId, newValue );
    if ( branchPerturbator->getLowShock() && (treePriority != 0.0) ){
        //the rebound might trigger a NNI
        InferenceNode* childNode = (InferenceNode*)getBranchDest(modifiedBranchId);
        //the node must not be a leaf
        bool NNIpossible = !childNode->isLeaf();
        //and must not define a cluster too
        if ( clustersTree && NNIpossible ){
            if (clustersTree->findByAncestor(childNode) != NULL){
                NNIpossible = false;
            }
        }
        if ( NNIpossible ){
            lastLocalChange = LOCAL_NNI;
            ++numberLocalNNI;
            pertSave.clear();
            pertSave.push_back( childNode );
            pertSave.push_back( (InferenceNode*)(childNode->getParent()) );
            chooseNNISwapNode( childNode );
            NNI( pertSave[0], pertSave[1], pertSave[2], pertSave[3] );
            assert(lnHR==0.0);
            return lnHR;
        }
    }
    lastLocalChange = LENGTH;
    return lnHR;
}


bool HeterogeneousMCMCTree::validateBranchLength( bool validation ){
    assert( modifiedBranchId != -1 );
    assert( lastLocalChange != INVALID_LOCAL );
    double oldValue;
    branchPerturbator->validatePerturbationValue( validation, &oldValue );
    //after a NNI induced by a negative branch length
    if (lastLocalChange == LOCAL_NNI){
        if (validation){
            ++acceptedLocalNNI;
            if (pertSave[0]->getParent()==pertSave[1]){
                pertSave[0]->saveRecursively(-1);
            }
            else{
                assert(pertSave[1]->getParent()==pertSave[0]);
                pertSave[1]->saveRecursively(-1);
            }
        }
        else{
            //come back
            NNI( pertSave[0], pertSave[1], pertSave[3], pertSave[2] );
            setBranchLength( modifiedBranchId, oldValue );
            if (pertSave[0]->getParent()==pertSave[1]){
                pertSave[0]->retrieveRecursively(-1);
            }
            else{
                assert(pertSave[1]->getParent()==pertSave[0]);
                pertSave[1]->retrieveRecursively(-1);
            }
        }
    }
    //after a normal branch modification
    else{
        assert(lastLocalChange == LENGTH);
        if (!validation){
            //we come back to the saved tree (the structure is still the same)
            setBranchLength( modifiedBranchId, oldValue );
            ((InferenceNode*)getBranchOrig(modifiedBranchId))->retrieveRecursively(-1);
        }
        else{
            ((InferenceNode*)getBranchOrig(modifiedBranchId))->saveRecursively(-1);
        }
    }
    lastLocalChange = INVALID_LOCAL;
    modifiedBranchId = -1;
    return true;
}




double HeterogeneousMCMCTree::globalTopologyChange(){
    assert( lastLocalChange == INVALID_LOCAL );
    assert( lastGlobalChange == INVALID_GLOBAL );

    double lnHastingRatio = 0.0;

    //SPR variables
    unsigned int branchId;
    unsigned int branchInsertId;
    InferenceNode* movedParent;
    InferenceNode* movedChild;
    InferenceNode* insertParent = NULL;
    InferenceNode* insertChild = NULL;
    bool found;
    InferenceNode* childNode = NULL;

    //model swap variables
    unsigned int oldModelId;
    unsigned int newModelId;
    vector< BasicNode * >::iterator iter;


    lastGlobalChange = (unsigned int)(randBox.ran() * 3);
    switch ( lastGlobalChange ){
      case SPR_PROP :            /* Subtree Pruning and Regrafting */
        ++numberSPRPerturbation;
        //choose the branch to move and the insertion point
        //(beware, special case if linked to the root)
        found = false;
        while (!found){
            branchId = (unsigned int)(getNumberBranches() * randBox.ran());
            movedParent = (InferenceNode*)getBranchOrig(branchId);
            movedChild = (InferenceNode*)getBranchDest(branchId);
            //find the brother of movedChild
            childNode = movedParent->getChild(0);
            if (childNode==movedChild){
                childNode = movedParent->getChild(1);
            }
            //if clusters are defined, collect a vector of possible insertion branches
            if (clustersTree){
                vector< pair<InferenceNode*, double> > candidateBranches;
                candidateBranches.clear();
                ClustersTreeNode * childCluster = findInsertion( movedChild );
                ClustersTreeNode * secCluster = NULL;
                //normal case if the branch is not linked to the root
                if (movedParent->getParent()){
                    secCluster = findInsertion( movedParent );
                }
                //if the branch is linked to the root we are interested in moving the branch movedChild--childNode
                else{
                    secCluster = findInsertion( childNode );
                }
                //special case where the child is a cluster ancestor,
                //the branch can "move" in both child and parent clusters
                //(ie the rooting point of the child cluster can change OR
                //the child cluster can move in the parent cluster
                if (childCluster!=secCluster){
                   //first add to branches all the possible rooting point
                    //for the childCluster
                    findClusterBranches( candidateBranches, childCluster );
                    /**                                            n1 --- > n2 ---> C
                     * example :    movedParent ---> movedChild  <
                     *                                             D
                     * we want to add the branches n1--->n2 and n2--->C to the possible
                     * moves (C defines another cluster and cannot be broken)
                     *                                                               */
                    //then add to branches all the possible branching point
                    //for the secCluster
                    findClusterBranches( candidateBranches, secCluster );
                    /**                                                             ...
                     * example :   C -->n1--->n2---> movedParent ---> movedChild  <
                     *                                                              ...
                     * we want to add the branches C--->n1 and n1--->n2 to the possible
                     * moves (C defines another cluster and cannot be broken)
                     *                                                               */
                }
                //if the branches is completely in a well defined cluster
                else{
                    findClusterBranches( candidateBranches, childCluster );
                }
                vector< InferenceNode* > allowedBranches;
                allowedBranches.reserve(candidateBranches.size());
                allowedBranches.clear();
                for ( vector< pair<InferenceNode*, double> >::iterator iter = candidateBranches.begin();
                      iter!=candidateBranches.end(); ++iter ){
                    if ( (iter->first != movedParent) &&
                         (iter->first != movedChild) &&
                         (iter->first->getParent() != movedParent) &&
                         (iter->first->getParent() != movedChild) &&
                         //we add a special impossible case for rooted tree here:
                         ( (movedParent!=root) || (iter->first->getParent()->getParent()!=root) ) ){
                        allowedBranches.push_back(iter->first);
                    }
                }
                 //if there is a place where we can move the branch then...
                found = (allowedBranches.size() != 0);
                if (found){
                    branchInsertId = ( int )( ( allowedBranches.size() * randBox.ran() ) );
                    insertChild = allowedBranches[branchInsertId];
                    insertParent = (InferenceNode*)insertChild->getParent();
                }
            }
            //usual case without clusters... there must be a possible move here (except 4 taxon case)
            else{
                do{
                    branchInsertId = ( branchId + 1 +
                    ( int )( ( getNumberBranches() - 1 ) * randBox.ran() ) ) %
                    getNumberBranches();
                    insertParent = (InferenceNode*)getBranchOrig(branchInsertId);
                    insertChild = (InferenceNode*)getBranchDest(branchInsertId);
                    //change the branch if the insertion is impossible, i.e
                    if ( (movedParent!=insertParent) && (movedParent!=insertChild) &&
                         (insertParent!=movedChild) &&
                         ( (movedParent!=root) || (insertParent->getParent()!=root) ) ){
                        found = true;
                    }
                } while ( (!found) && (getNumberBranches() > 6) );
            }
        }
        //save in pertSave the branch we are moving and where it comes from
        pertSave.clear();
        pertSave.push_back( movedParent );
        //if movedParent is the root, one might have to change moved child
        //we want the moved branch to be in the path from the insertion point to the root
        //--> swap the two child node
        if ( (!movedParent->getParent()) && (!insertParent->isDescendant(movedChild)) ){
            std::swap(childNode, movedChild);
        }
        pertSave.push_back( movedChild );

        //case 1, insertParent was below movedChild (or is after the first
        //manip we did in the case of movedParent==root)
        if ( insertParent->isDescendant(movedChild) ){
            pertSave.push_back( movedChild->getChild(0) );
            pertSave.push_back( movedChild->getChild(1) );
            //save the position of movedParent in the branch (movedParent->getParent()---->childOrig)
            oldFrac = pertSave[2]->getParentDistance() /
                  (pertSave[2]->getParentDistance()+pertSave[3]->getParentDistance());
        }
        //otherwise, it is guaranteed that movedParent != root and that
        //insertParent is not in a subTree started by movedChild
        else{
            pertSave.push_back( (InferenceNode*)(movedParent->getParent()) );
            pertSave.push_back( childNode );
            oldFrac = movedParent->getParentDistance() /
                    (movedParent->getParentDistance()+childNode->getParentDistance());
        }
        //Note that in both case we saved in oldFrac the ratio:
        //d(pertSave[2],branchMoved) / ( d(pertSave[2],branchMoved) + d(pertSave[3],branchMoved) )
        lnHastingRatio = log( insertBranch( insertParent, insertChild,
                                            movedParent, movedChild, randBox.ran() ) );
      break;
      case NNI_PROP :            /* Nearest Neighbor Interchange */
       ++numberNNIPerturbation;
        childNode = NULL;
        while ( childNode == NULL ){
            branchId = ( unsigned int )( getNumberBranches() * randBox.ran() );
            childNode = (InferenceNode*)(getBranchDest(branchId));
            //if childNode is not ok, put it back to NULL for the next try
            //NNi not allowed if the branch is a terminal branch
            if (childNode->isLeaf()){
                childNode = NULL;
            }
            //NNi not allowed if the branch is linked to a cluster
            else{
                if ( clustersTree && clustersTree->findByAncestor(childNode) ){
                    childNode = NULL;
                }
            }
        }
        pertSave.clear();
        pertSave.push_back( (InferenceNode*)(getBranchDest(branchId)) );
        pertSave.push_back( (InferenceNode*)(getBranchOrig(branchId)) );
        chooseNNISwapNode( (InferenceNode*)getBranchDest(branchId) );
        NNI( pertSave[0], pertSave[1], pertSave[2], pertSave[3] );
        lnHastingRatio = 0.0;
      break;
      case SWAP_PROP :            /* Model Swap */
        //add one to the counter of swap proposals
        swapPerturbator->addPerturbation();
        modelSave.clear();
        iter = nodeRefVector.begin();
        ++iter;
        while ( iter != nodeRefVector.end() ){
            //for each node, prob swap = step of swapPerturbator
            if ( randBox.ran() < swapPerturbator->getStep() ){
                //save the old model
               Model* oldModel = ((InferenceNode*)(*iter))->getModel();
               oldModelId = modelMap[oldModel];
               modelSave.push_back( pair<InferenceNode*, Model*>
                    ( (InferenceNode*)(*iter), oldModel ) );
               //choose a new model (different from the old one)
               newModelId = (int)(randBox.ran()*(numberModels-1));
               if ( newModelId >= oldModelId ){
                   ++newModelId;
               }
               ((InferenceNode*)(*iter))->setModel(
                       pheterogeneous->getModel( newModelId ) );
               ((InferenceNode*)(*iter)->getParent())->invalidateRecursively(-1);
            }
            ++iter;
        }
        lnHastingRatio = 0.0;
      break;
      default :
        assert(0);
      break;
    }
    return lnHastingRatio;
}

void HeterogeneousMCMCTree::chooseNNISwapNode( InferenceNode* nodeChild ){

    int swapId;

    assert(!nodeChild->isLeaf());
    InferenceNode* nodeParent = (InferenceNode*)(nodeChild->getParent());
    assert(nodeParent);

    //selection of the first node to swap (a children of nodeChild)
    swapId = (int)(nodeChild->getNumberChildren() * randBox.ran());
    pertSave.push_back( nodeChild->getChild( swapId ) ) ;

    //selection of the second node to swap (adjacent to nodeParent)
    //it cannot be nodeChild but it can be nodeParent's parent if it exists
    InferenceNode* swapNode = NULL;
    while (!swapNode){
        swapId = (int)((nodeParent->getNumberChildren()) * randBox.ran());
        swapNode = nodeParent->getChild( swapId );
        //if nodeChild is chosen, replace swapNode2 by the parent
        //if nodeParent is root, another node (there is only one choice for
        //binary tree) will be chosen
        if (swapNode == nodeChild){
            swapNode = (InferenceNode*)(nodeParent->getParent());
        }
    }
    pertSave.push_back( swapNode );
}


bool HeterogeneousMCMCTree::validateTopologyChange( bool validation ){
    assert( lastGlobalChange != INVALID_GLOBAL );
    //modified and modifiedStructure have been invalidated,
    //whatever the value of validation
    switch ( lastGlobalChange ) {
      case SPR_PROP:
        if (validation){
            ++acceptedSPRPerturbation;
            saveNodes(-1);
        }
        else{
            //come back
            insertBranch( pertSave[2], pertSave[3], pertSave[0], pertSave[1], oldFrac );
            retrieveNodes(-1);
        }
      break;
      case NNI_PROP:
        if (validation){
            ++acceptedNNIPerturbation;
            if (pertSave[0]->getParent()==pertSave[1]){
                pertSave[0]->saveRecursively(-1);
            }
            else{
                assert(pertSave[1]->getParent()==pertSave[0]);
                pertSave[1]->saveRecursively(-1);
            }
        }
        else{
            //come back
            NNI( pertSave[0], pertSave[1], pertSave[3], pertSave[2] );
            if (pertSave[0]->getParent()==pertSave[1]){
                pertSave[0]->retrieveRecursively(-1);
            }
            else{
                assert(pertSave[1]->getParent()==pertSave[0]);
                pertSave[1]->retrieveRecursively(-1);
            }
        }
      break;
      case SWAP_PROP:
        //the structure has not changed
        //warn the perturbator about the outcome (for adaptation of the step)
        swapPerturbator->validatePerturbation(validation);
        if (validation){
            for ( unsigned int i = 0; i < modelSave.size(); ++i ){
                ((InferenceNode*)((modelSave[i].first)->getParent()))->saveRecursively(-1);
            }
        }
        else{
            for ( unsigned int i = 0; i < modelSave.size(); ++i ){
                (modelSave[i].first)->setModel( modelSave[i].second );
                ((InferenceNode*)((modelSave[i].first)->getParent()))->retrieveRecursively(-1);
            }
        }
      break;
      default:
        assert(0);
      break;
    }
    lastGlobalChange = INVALID_GLOBAL;
    return true;
}


void HeterogeneousMCMCTree::stopBurn(){
    MCMCTreeBasic::stopBurn();
    numberSPRPerturbation = 0;
    acceptedSPRPerturbation = 0;
    numberNNIPerturbation = 0;
    acceptedNNIPerturbation = 0;
    numberLocalNNI = 0;
    acceptedLocalNNI = 0;
    swapPerturbator->stopBurn();
    branchPerturbator->stopBurn();
}

void HeterogeneousMCMCTree::initialiseMCMC( ParametersSet & parameters ){
    numberSPRPerturbation = 0;
    acceptedSPRPerturbation = 0;
    numberNNIPerturbation = 0;
    acceptedNNIPerturbation = 0;

    numberLocalNNI = 0;
    acceptedLocalNNI = 0;

    lastGlobalChange = INVALID_GLOBAL;
    lastLocalChange = INVALID_LOCAL;

    swapPerturbator = new PerturbatorBase("Swap", parameters,
                                            "Swap", .05, .01, .9,
                                            "initial probability");
    PriorField branchPriors(parameters.stringParameter("Branch lengths, prior"));
    branchPerturbator = new PerturbatorHelper( hyperPerturbator, *this, UpdateMessage::BRANCH,
                                               "Branch lengths", branchPriors,
                                               parameters, "Branch lengths",
                                               .05, 0.0, 0.5, "initial step" );
    MCMCTreeBasic::initialiseMCMC(parameters);
}

void HeterogeneousMCMCTree::printPerturbationParameters(ostream& outputStream){
    outputStream << setprecision(5);
    outputStream << "SPR acceptance rate                    =  "
                 << (double)acceptedSPRPerturbation*100.0 / (double)numberSPRPerturbation
                 << '%' << endl;
    outputStream << "Number of SPR proposals                =  "
                 << numberSPRPerturbation << endl;
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

    swapPerturbator->printPerturbationParameters( outputStream );
    branchPerturbator->printPerturbationParameters( outputStream );
    MCMCTreeBasic::printPerturbationParameters( outputStream );
}


void HeterogeneousMCMCTree::NNI( InferenceNode* node1, InferenceNode* node2,
                            InferenceNode* swapNode1, InferenceNode* swapNode2 ){
    InferenceNode* nodeParent = node2;
    InferenceNode* nodeChild = node1;
    InferenceNode* swapNodeParent = swapNode2;
    InferenceNode* swapNodeChild = swapNode1;
    if (node1->getParent() != node2){
        nodeParent = node1;
        nodeChild = node2;
        assert(nodeChild->getParent() == nodeParent);
        swapNodeParent = swapNode1;
        swapNodeChild = swapNode2;
    }


//    cout << "NNI branch" << nodeParent->toString(true) << " -->  " << nodeChild->toString(true) << endl;
//    cout << "Swap P " << swapNodeParent->toString(true) << endl;
//    cout << "Swap C " << swapNodeChild->toString(true) << endl;

    //swap
    if ( swapNodeParent != (InferenceNode*)(nodeParent->getParent()) ){
        swap( swapNodeParent, swapNodeChild );
        //parameters are automatically updated
    }
    else{
        //if node parent was the MRCA of a cluster we need to correct that
        //before the swap
        if (clustersTree){
            ClustersTreeNode* clusterParent = clustersTree->findByAncestor(nodeParent);
            if (clusterParent){
                clusterParent->cluster->setAncestor(nodeChild);
            }
        }
        //distance swapNodeParent->nodeParent
        double distance1 = nodeParent->getParentDistance();
        //distance nodeParent->nodeChild
        double distance2 = nodeChild->getParentDistance();
        //distance nodeChild->swapNodeChild
        double distance3 = swapNodeChild->getParentDistance();
        //break swapNodeParent->nodeParent
        swapNodeParent->removeChild(nodeParent);
        //break nodeParent->nodeChild
        nodeParent->removeChild( nodeChild );
        //break nodeChild->swapNodeChild
        nodeChild->removeChild( swapNodeChild );
        //attach swapNodeParent->nodeChild with distance1
        swapNodeParent->addChild(nodeChild, distance1);
        //attach nodeChild->nodeParent with distance2
        nodeChild->addChild( nodeParent, distance2 );
        //attach nodeParent->swapNodeChild with distance3
        nodeParent->addChild( swapNodeChild, distance3 );
        //swap the models of nodeChild and nodeParent
        Model* pmodel = nodeParent->getModel();
        nodeParent->setModel( nodeChild->getModel() );
        nodeChild->setModel( pmodel );
        //update parameters recursively in the root direction from node2
        nodeParent->updateNodesCountRecursively();
        nodeParent->updateIdenticalRecursively();
        nodeParent->invalidateRecursively(-1);
    }
    createIndex();
}

double HeterogeneousMCMCTree::insertBranch( InferenceNode* insert1, InferenceNode* insert2,
                                       InferenceNode* moved1, InferenceNode* moved2,
                                       double frac ){
    InferenceNode* insertParent = insert1;
    InferenceNode* insertChild = insert2;
    InferenceNode* movedParent = moved1;
    InferenceNode* movedChild = moved2;

    if (insert2->getParent()!=insert1){
        assert(insert1->getParent()==insert2);
        insertParent = insert2;
        insertChild = insert1;
        frac = 1.0 - frac;
    }
    if (moved2->getParent()!=moved1){
        assert(moved1->getParent()==moved2);
        movedParent = moved2;
        movedChild = moved1;
    }

    assert(insertParent!=movedParent);
    assert(insertChild!=movedParent);
    assert(insertParent!=movedChild);

//    cout << "Move" << movedParent->toString(true) << " -->  " << movedChild->toString(true) << endl;
//    cout << "Insert" << insertParent->toString(true) << " -->  " << insertChild->toString(true) << endl;

    double distance = insertChild->getParentDistance();
    double oldDistance;

    if (movedParent==root){
        assert(insertParent->isDescendant(movedChild));
    }

    //special case when the branch to move is between the root and the
    //insert branch, we have some trouble with the model assigned to
    //each branch there. Clusters are not affected by the move
    //there is no need to handle something for them here
    if ( insertParent->isDescendant(movedChild) ){
        /*                         subNodeY->..
         *root->....->movedChild <
         *                         subNodeX->...->insertParent->insertChild->....
         * ************************** CHANGED TO ****************************
         *                          insertParent->...->...->subNodeX->subNodeY->..
         *root->....->movedChild <
         *                          insertChild->....
         * we have to swap model for the nodes between subNodeX and
         * insertParent                                                       */
        //detach the two subtrees from the child of the branch to move
        assert(movedChild->getNumberChildren() == 2);
        InferenceNode* subNode1 = (InferenceNode*)(movedChild->getChild(0));
        InferenceNode* subNode2 = (InferenceNode*)(movedChild->getChild(1));
        oldDistance = subNode1->getParentDistance() +
                             subNode2->getParentDistance();
        movedChild->removeChild( subNode1 );
        movedChild->removeChild( subNode2 );

        //break the insert branch(subNodeX can be subNode1 or subNode2)
        //subNodeX---->...---->insertParent ---/ /---->insertChild
        insertParent->removeChild( insertChild );
        //revert the tree from insertOrig to subNodeX (X = 1 or 2)
        //subNodeX->...->insertParent  ==>  insertParent->...->subNodeX
        //and make sur each branch will keep its model
        insertParent->makeParent( NULL );
        //restore the insertBranch around the child of the moved branch
        movedChild->addChild( insertChild, (1.0-frac)*distance );
        movedChild->addChild( insertParent, frac*distance );
        //                          insertParent->...->...->subNodeX
        //root->....->movedChild <
        //                          insertChild->....

        //attach subNode1 to subNode2, the new parent subNodeX is the one
        //which was reverted in the previous process (the one in the subtree
        //which contains the insert branch)
        //therefore subNodeX is the node which has been given a parent
        //during the call to makeParent
        InferenceNode* parentNode;
        InferenceNode* childNode;
        if ( subNode1->getParent() ){
            parentNode = subNode1;
            childNode = subNode2;
        }
        else{
            parentNode = subNode2;
            childNode = subNode1;
        }
        parentNode->addChild( childNode, oldDistance );
        //                          insertParent->...->...->subNodeX->subNodeY->..
        //root->....->movedChild <
        //                          insertChild->....

        //restore to each branch its previous model
        //the branches movedChild->subNodeX and movedChild->subNodeY
        //have merged into subNodeX->subNodeY and we keep for them the model
        //which was for movedChild->subNodeY.
        //the branch insertParent->insertChild has split. Its model is used for
        //movedChild->insertChild. movedChild->insertParent has no model
        //and we use the model which was used for movedChild->subNodeX
        //so that the move is reversible and with a 'limited' amount of change
        Model* movedModel = parentNode->getModel();
        InferenceNode* temp = parentNode;
        while (temp!=insertParent){
            temp->setModel(((InferenceNode*)temp->getParent())->getModel());
            temp = (InferenceNode*)temp->getParent();
        }
        insertParent->setModel(movedModel);

        //update parameters :
        //subNodeY is ok (it was detached and it remains a valid subtree
        //with count, identical and invalidate OK).
        //insertDest is OK
        //nodes movedChild (and ancestors) are not valid anymore
        //parent node (ie subNodeX) is not valid too
        parentNode->updateNodesCountRecursively();
        parentNode->invalidateRecursively(-1);
        parentNode->updateIdenticalRecursively();
    }
    //otherwise normal and simple case
    else{
        /*             ... -> movedParent->movedChild->..
         *      root <
         *             ... ->insertParent->insertChild->..
         *  OR
         *      root->...->insertParent->insertChild->..-> movedParent->movedChild->..
         * ************************** CHANGED TO ****************************
         *                                         movedChild->..
         *root->....->insertParent -> movedParent<
         *                                         insertChild->....        */
        assert(movedParent->getNumberChildren()==2);
        InferenceNode* childNode = (InferenceNode*)movedParent->getChild(0);
        if ( childNode == movedChild ){
            childNode = (InferenceNode*)movedParent->getChild(1);
        }
        InferenceNode* parentNode = (InferenceNode*)movedParent->getParent();
        /*                           movedChild
         * parentNode->movedParent <
         *                           childNode                              */

        //movedParent can define a cluster in a special case where movedParent is insertParent ancestor
        //in such a case the cluster ancestor must be moved to childNode
        //insertChild can define a cluster too, if moved parent was a descendant then we have to translate
        //the cluster too...
        if( clustersTree ){
            ClustersTreeNode* clustersTreeNode = clustersTree->findByAncestor( movedParent );
            if (clustersTreeNode){
                assert(insertParent->isDescendant(movedParent));
                clustersTreeNode->cluster->setAncestor( childNode );
            }
            else{
                ClustersTreeNode* clustersTreeNode = clustersTree->findByAncestor( insertChild );
                if (clustersTreeNode){
                    if(movedParent->isDescendant(insertChild)){
                        clustersTreeNode->cluster->setAncestor( movedParent );
                    }
                }
            }
        }

        oldDistance = childNode->getParentDistance() +
                      movedParent->getParentDistance();
        movedParent->removeChild( childNode );
        parentNode->removeChild( movedParent );
        parentNode->addChild( childNode, oldDistance );
        //insert movedParent in insertBranch
        insertParent->removeChild(insertChild);
        insertParent->addChild(movedParent, frac*distance);
        movedParent->addChild(insertChild, (1.0-frac)*distance);
        //update...
        InferenceNode* ancestor = (InferenceNode*)
                getMostRecentAncestor( movedParent, parentNode );
        movedParent->updateNodesCountRecursively(ancestor);
        movedParent->invalidateRecursively( -1, ancestor);
        movedParent->updateIdenticalRecursively(ancestor);
        parentNode->updateNodesCountRecursively();
        parentNode->invalidateRecursively(-1);
        parentNode->updateIdenticalRecursively();
    }
    createIndex();
    return distance/oldDistance;
}

void HeterogeneousMCMCTree::loadDataAndModel
        ( SequenceTable * psequenceTable, Model * pmodel ){
    RootedTree::loadDataAndModel( psequenceTable, pmodel );
    pheterogeneous = dynamic_cast<Heterogeneous*>(pmodel);
    if ( !pheterogeneous ){
        cerr << "You must use a heterogeneous model with the"
             << " heterogeneous tree" << endl;
        exit(EXIT_FAILURE);
    }
    numberModels = pheterogeneous->getNumberModels();
    for ( unsigned int i = 0; i < numberModels; ++i ){
        modelMap[pheterogeneous->getModel(i)] = i;
    }

    //if the tree is already constructed
    if (nodeRefVector.size()){
        initialiseNodeModel();
    }
}


void HeterogeneousMCMCTree::update( UpdateMessage* message ){
    int cat=-1;
    int model=-1;
    assert(message->hasType(UpdateMessage::MODEL_TYPE));
    if (message->hasFlag(UpdateMessage::MODEL_PRIOR_FLAG)){
        //cout << "message prior from the model" << endl;
        return;
    }
    if (message->hasFlag(UpdateMessage::SYMBOL_CATEGORY)){
        cat = (int)message->modelMsg.symbolCategory;
    }
#ifdef DEBUG1
    cout << "Heterogeneous MolClk tree receive update from symbol " << cat << endl;
#endif
    if (message->hasFlag(UpdateMessage::MODEL_FLAG)){
        model = (int)message->modelMsg.model;
    }
#ifdef DEBUG1
    cout << "Heterogeneous tree receive update from model " << model << endl;
#endif
    //if one ancestral frequency was modified then return
    //(there is no mechanism yet to preserve computation at the root...)
    //so it will be computed again next time
    if ( model == (int)pheterogeneous->getNumberModels() ){
        return;
    }
    //if there is no model information invalidate the correct symbol category (might be all if -1)
    if ( model == -1 ){
        invalidateNodes( cat );
        return;
    }
    //else, if there is one specific model...
    vector<BasicNode*>::iterator iter = nodeRefVector.begin();
    ++iter;
    while ( iter != nodeRefVector.end() ){
        if ( modelMap[((InferenceNode*)(*iter))->getModel()] == model ){
            ((InferenceNode*)((*iter)->getParent()))->invalidateRecursively( cat );
        }
        ++iter;
    }
}

void HeterogeneousMCMCTree::retrieveNodesAux(){
    MCMCTreeBasic::retrieveNodesAux();
}

void HeterogeneousMCMCTree::saveNodesAux(){
    MCMCTreeBasic::saveNodesAux();
}


void HeterogeneousMCMCTree::getAllParameters( vector<double>& params ) const{
    params.clear();
    //get the branch lengths
    MCMCTreeBasic::getAllParameters( params );
    vector< double > modelId;
    vector<BasicNode*>::const_iterator iterNode = nodeRefVector.begin();
    ++iterNode;
    //get the model ID associated with the branch
    while ( iterNode != nodeRefVector.end() ){
        map< Model*, int>::const_iterator search =
                modelMap.find( ((InferenceNode*)(*iterNode))->getModel() );
        params.push_back( (double)((*search).second) );
        ++iterNode;
    }
    assert( params.size() == getNumberTreeParameters());
}

void HeterogeneousMCMCTree::setAllParameters( const vector<double>& params){
    assert( params.size() == getNumberTreeParameters());
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp;
    //first set basic parameters
    nbp = MCMCTreeBasic::getNumberTreeParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    MCMCTreeBasic::setAllParameters( p );
    //then set the model for each node
    vector<BasicNode*>::iterator iterNode = nodeRefVector.begin();
    ++iterNode;
    while ( iterNode != nodeRefVector.end() ){
        ((InferenceNode*)(*iterNode))->setModel(
                pheterogeneous->getModel( (int)(*iter) ) );
        ++iterNode;
        ++iter;
    }
}

unsigned int HeterogeneousMCMCTree::getNumberTreeParameters() const{
    unsigned int numberParameters = MCMCTreeBasic::getNumberTreeParameters();
    //add the number of branches to store the model associated with each of
    //them
    numberParameters += getNumberBranches();
    return numberParameters;
}

void HeterogeneousMCMCTree::getAllPerturbationParameters( vector<double>& params ) const{
    MCMCTreeBasic::getAllPerturbationParameters(params);
    params.push_back((double)numberSPRPerturbation);
    params.push_back((double)acceptedSPRPerturbation);
    params.push_back((double)numberNNIPerturbation);
    params.push_back((double)acceptedNNIPerturbation);
    params.push_back((double)numberLocalNNI);
    params.push_back((double)acceptedLocalNNI);
    vector<double> p;
    branchPerturbator->getAllPerturbationParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
    swapPerturbator->getAllPerturbationParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
    assert( params.size() == getNumberPerturbationParameters() );
}

void HeterogeneousMCMCTree::setAllPerturbationParameters( const vector<double>& params){
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
    numberSPRPerturbation = (unsigned int)(*iter);
    ++iter;
    acceptedSPRPerturbation = (unsigned int)(*iter);
    ++iter;
    numberNNIPerturbation = (unsigned int)(*iter);
    ++iter;
    acceptedNNIPerturbation = (unsigned int)(*iter);
    ++iter;
    numberLocalNNI = (unsigned int)(*iter);
    ++iter;
    acceptedLocalNNI = (unsigned int)(*iter);
    ++iter;
    nbp = branchPerturbator->getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    branchPerturbator->setAllPerturbationParameters ( p );
    nbp = swapPerturbator->getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    swapPerturbator->setAllPerturbationParameters ( p );
}

unsigned int HeterogeneousMCMCTree::getNumberPerturbationParameters() const{
    unsigned int numberPerturbationParameters = 6 + MCMCTreeBasic::getNumberPerturbationParameters();
    numberPerturbationParameters += branchPerturbator->getNumberPerturbationParameters();
    numberPerturbationParameters += swapPerturbator->getNumberPerturbationParameters();
    return numberPerturbationParameters;
}

void HeterogeneousMCMCTree::getAllPriorParameters( vector<double>& params ) const{
    MCMCTreeBasic::getAllPriorParameters(params);
    vector<double> p;
    branchPerturbator->getAllPriorParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
}

void HeterogeneousMCMCTree::setAllPriorParameters( const vector<double>& params){
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
    nbp = branchPerturbator->getNumberPriorParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    branchPerturbator->setAllPriorParameters( p );
    assert( iter == params.end() );
}

unsigned int HeterogeneousMCMCTree::getNumberPriorParameters() const{
    unsigned int numberPriorParameters = MCMCTreeBasic::getNumberPriorParameters();
    numberPriorParameters += branchPerturbator->getNumberPriorParameters();
    return numberPriorParameters;
}

double HeterogeneousMCMCTree::getLnPrior() const{
    double lnPrior = MCMCTreeBasic::getLnPrior();
    for ( unsigned int branchId = 0; branchId < getNumberBranches(); ++branchId ){
        lnPrior += branchPerturbator->getLnPriorValue(
                getBranchLength( branchId ) );
    }
    return lnPrior;
}

void HeterogeneousMCMCTree::constructFromString( string stringTree ){
    RootedTree::constructFromString( stringTree );
    //if the model has been loaded
    if(pheterogeneous){
        initialiseNodeModel();
    }
}

void HeterogeneousMCMCTree::constructRandomly( double maxBranchLength ){
    RootedTree::constructRandomly( maxBranchLength );
    //if the model has been loaded
    if(pheterogeneous){
        initialiseNodeModel();
    }
}

void HeterogeneousMCMCTree::initialiseNodeModel(){
    vector< BasicNode * >::iterator iter = nodeRefVector.begin();
    while ( iter != nodeRefVector.end() ){
        if(*iter!=root){
            ((InferenceNode*)(*iter))->setModel( pheterogeneous->getModel(
                    (int)(randBox.ran()*(double)numberModels) ) );
        }
        ++iter;
    }
}

