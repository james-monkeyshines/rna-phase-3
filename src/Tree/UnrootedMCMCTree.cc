#include "Tree/UnrootedMCMCTree.h"

#include "Util/ParametersSet.h"

#include <iostream>
#include <iomanip>
#include <assert.h>

#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

//register the tree to the tree factory with the name Unrooted MCMC tree
UnrootedMCMCTree UnrootedMCMCTree::prototype("Unrooted MCMC tree");

void UnrootedMCMCTree::constructFromString( string stringTree ){
    UnrootedTree::constructFromString( stringTree );
    checkBinary();
}


UnrootedMCMCTree::UnrootedMCMCTree( const string & registrationName ) :
InferenceTree(), MCMCTreeBasic(registrationName), UnrootedTree(){
}

UnrootedMCMCTree::UnrootedMCMCTree(ParametersSet& parameters) :
InferenceTree(), MCMCTreeBasic(parameters), UnrootedTree(parameters){
}

UnrootedMCMCTree::~UnrootedMCMCTree(){
    if (branchPerturbator){
        delete branchPerturbator;
    }
}

MCMCTree* UnrootedMCMCTree::clone(ParametersSet& parameters) const{
    return new UnrootedMCMCTree(parameters);
}

double UnrootedMCMCTree::changeBranchLength(){
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


bool UnrootedMCMCTree::validateBranchLength( bool validation ){
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




double UnrootedMCMCTree::globalTopologyChange(){
    assert( lastLocalChange == INVALID_LOCAL );
    assert( lastGlobalChange == INVALID_GLOBAL );

    double lnHastingRatio = 0.0;
    unsigned int branchId;
    unsigned int branchInsertId;
    InferenceNode* movedParent;
    InferenceNode* movedChild;
    InferenceNode* insertParent = NULL;
    InferenceNode* insertChild = NULL;
    bool found;
    InferenceNode* childNode;

    lastGlobalChange = (int)(randBox.ran() * 2);
    switch ( lastGlobalChange ){
      case SPR_PROP:                /* subtree pruning and regrafting */
        ++numberSPRPerturbation;
        //choose the branch to move and the insertion point
        found = false;
        while (!found){
            //the branch to move first
            branchId = (unsigned int)(getNumberBranches() * randBox.ran());
            movedParent = (InferenceNode*)getBranchOrig(branchId);
            movedChild = (InferenceNode*)getBranchDest(branchId);
            //if clusters are defined, collect a vector of possible insertion branches
            if (clustersTree){
                vector< pair<InferenceNode*, double> > branches;
                ClustersTreeNode * childCluster = findInsertion( movedChild );
                ClustersTreeNode * parentCluster = findInsertion( movedParent );
                //special case where the child is a cluster ancestor,
                //the branch can "move" in both child and parent clusters
                //(ie the rooting point of the child cluster can change OR
                //the child cluster can move in the parent cluster
                if (childCluster!=parentCluster){
                    //first add to branches all the possible rooting point
                    //for the childCluster
                    findClusterBranches( branches, childCluster );
                    /**                                            n1 --- > n2 ---> C
                     * example :    movedParent ---> movedChild  <
                     *                                             D
                     * we want to add the branches n1--->n2 and n2--->C to the possible
                     * moves (C defines another cluster and cannot be broken)
                     *                                                               */
                    //then add to branches all the possible branching point
                    //for the parentCluster
                    findClusterBranches( branches, parentCluster );
                    /**                                                             ...
                     * example :   C -->n1--->n2---> movedParent ---> movedChild  <
                     *                                                              ...
                     * we want to add the branches C--->n1 and n1--->n2 to the possible
                     * moves (C defines another cluster and cannot be broken)
                     *                                                               */
                }
                else{
                    findClusterBranches( branches, parentCluster );
                }
                vector< pair<InferenceNode*, double> >::iterator iter = branches.end();
                while (iter!=branches.begin()){
                    --iter;
                    if ( ((*iter).first == movedParent) ||
                         ((*iter).first == movedChild) ||
                         ((*iter).first->getParent() == movedParent) ||
                         ((*iter).first->getParent() == movedChild) ){
                        branches.erase(iter);
                    }
                }
                //if there is a place where we can move the branch then...
                found = (branches.size() != 0);
                if (found){
                    branchInsertId = ( int )( ( branches.size() * randBox.ran() ) );
                    insertChild = branches[branchInsertId].first;
                    insertParent = (InferenceNode*)insertChild->getParent();
                }
            }
            //usual case without clusters...
            else{
                //in the 4-taxon case there is a branch (the middle one) which cannot
                //move anywhere so becareful and assume we might have to return false
                //in this case too
                do{
                    branchInsertId = ( branchId + 1 +
                        ( int )( ( getNumberBranches() - 1 ) * randBox.ran() ) ) %
                        getNumberBranches();
                    insertParent = (InferenceNode*)getBranchOrig(branchInsertId);
                    insertChild = (InferenceNode*)getBranchDest(branchInsertId);
                    //is it an acceptable destination branch?
                    if ( (movedParent!=insertParent) &&
                         (movedParent!=insertChild) &&
                         (insertParent!=movedChild) ){
                        //there is no need to check whether moved branch == branch dest
                        //since insertChild cannot be equal to movedChild here
                        found = true;
                    }
                }while ( (!found) && (getNumberBranches() > 5) );
                //we must find a branch if this is not the 4 taxon case
            }
        }
        //save in pertSave the branch we are moving
        pertSave.clear();
        pertSave.push_back( movedParent );
        pertSave.push_back( movedChild );
        // save in pertSave the branch from where the moved branch
        // is extracted (it can be 1)the branch that contains movedChild,
        // 2)the branch that contains movedParent
        //special case 1) when the branch to move is between the root and the
        //insert branch
        InferenceNode* node;
        if ( insertParent->isDescendant(movedChild) ){
            pertSave.push_back( movedChild->getChild(0) );
            pertSave.push_back( movedChild->getChild(1) );
            //save the position of movedChild in the branch (child1----child2)
            oldFrac = pertSave[2]->getParentDistance() /
                  (pertSave[2]->getParentDistance()+pertSave[3]->getParentDistance());
        }
        else{
            //special case 2) when the branch to move is linked to the root
            if (!movedParent->getParent()){
                node = movedParent->getChild(0);
                if (node==movedChild){
                    node = movedParent->getChild(2);
                }
                pertSave.push_back( node );
                node = movedParent->getChild(1);
                if (node==movedChild){
                    node = movedParent->getChild(2);
                }
                pertSave.push_back( node );
                //save the position of movedParent (ie the root) in the branch defined by its two other children
                oldFrac = pertSave[2]->getParentDistance() /
                      (pertSave[2]->getParentDistance()+pertSave[3]->getParentDistance());

            }
            //normal case
            else{
                node = movedParent->getChild(0);
                if (node==movedChild){
                    node = movedParent->getChild(1);
                }
                pertSave.push_back( node );
                pertSave.push_back( (InferenceNode*)(movedParent->getParent()) );
                //save the position of movedParent in the branch (movedParent->getParent()---->childOrig)
                oldFrac = pertSave[2]->getParentDistance() /
                      (pertSave[2]->getParentDistance()+movedParent->getParentDistance());
            }
        }
        /*                           movedchild
         *   parent -> movedparent <
         *                           childOrig
         */
        lnHastingRatio = log( insertBranch( insertParent, insertChild,
                                            movedParent, movedChild, randBox.ran() ) );
      break;
      case NNI_PROP :                      /* Nearest Neighbor Interchange */
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
                if ( clustersTree &&
                    (clustersTree->findByAncestor(childNode) != NULL) ){
                    childNode = NULL;
                }
            }
        }
        pertSave.clear();
        pertSave.push_back( childNode );
        pertSave.push_back( (InferenceNode*)(getBranchOrig(branchId)) );
        chooseNNISwapNode( (InferenceNode*)getBranchDest(branchId) );
        NNI( pertSave[0], pertSave[1], pertSave[2], pertSave[3] );
        lnHastingRatio = 0.0;
      break;
      default :
        assert(0);
      break;
    }
    return lnHastingRatio;
}

void UnrootedMCMCTree::chooseNNISwapNode( InferenceNode* nodeChild ){

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
        //if nodeParent is root, another node will be chosen
        if (swapNode == nodeChild){
            swapNode = (InferenceNode*)(nodeParent->getParent());
        }
    }
    pertSave.push_back( swapNode ) ;
}

bool UnrootedMCMCTree::validateTopologyChange( bool validation ){
    assert( lastGlobalChange != INVALID_GLOBAL );
    switch ( lastGlobalChange ){
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
      default:
        assert(0);
      break;
    }
    lastGlobalChange = INVALID_GLOBAL;
    return true;
}


void UnrootedMCMCTree::stopBurn(){
    MCMCTreeBasic::stopBurn();
    numberSPRPerturbation = 0;
    acceptedSPRPerturbation = 0;
    numberNNIPerturbation = 0;
    acceptedNNIPerturbation = 0;
    numberLocalNNI = 0;
    acceptedLocalNNI = 0;
    branchPerturbator->stopBurn();
}

void UnrootedMCMCTree::initialiseMCMC( ParametersSet & parameters ){
    numberSPRPerturbation = 0;
    acceptedSPRPerturbation = 0;
    numberNNIPerturbation = 0;
    acceptedNNIPerturbation = 0;

    numberLocalNNI = 0;
    acceptedLocalNNI = 0;

    lastGlobalChange = INVALID_GLOBAL;
    lastLocalChange = INVALID_LOCAL;

    PriorField branchPriors(parameters.stringParameter("Branch lengths, prior"));
    branchPerturbator = new PerturbatorHelper( hyperPerturbator, *this, UpdateMessage::BRANCH,
                                               "Branch lengths", branchPriors,
                                               parameters, "Branch lengths",
                                               .05, 0.0, 0.5, "initial step" );
    MCMCTreeBasic::initialiseMCMC(parameters);
}


void UnrootedMCMCTree::printPerturbationParameters(ostream& outputStream){
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

    branchPerturbator->printPerturbationParameters( outputStream );
    MCMCTreeBasic::printPerturbationParameters( outputStream );
}


void UnrootedMCMCTree::NNI( InferenceNode* node1, InferenceNode* node2,
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
        //swapNodeParent might have been the outgroup
        if ( !installOutgroup() ){
            createIndex();
        }
        //the ancestor of the root cluster should have been updated
        //during installOutgroup if necessarry
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
        //update parameters recursively in the root direction from node2
        nodeParent->updateNodesCountRecursively();
        nodeParent->updateIdenticalRecursively();
        nodeParent->invalidateRecursively(-1);
        if ( !installOutgroup() ){
            createIndex();
        }
        else{
            assert(clustersTree);
        }
    }
}

double UnrootedMCMCTree::insertBranch( InferenceNode* insert1, InferenceNode* insert2,
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

    //special case when the branch to move is between the root and the
    //insert branch, in such a case, clusters are properly preserved
    //without intervention
    if ( insertParent->isDescendant(movedChild) ){
        /*                         subNodeY->..
         *root->....->movedChild <
         *                         subNodeX->...->insertParent->insertChild->....
         * ************************** CHANGED TO ****************************
         *                          insertParent->...->...->subNodeX->subNodeY->..
         *root->....->movedChild <
         *                          insertChild->....                         */
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

        //update parameters :
        //subNodeY is ok (it was detached and it remains a valid subtree
        //with count, identical and invalidate OK).
        //insertDest is OK
        //nodes movedChild (and ancestors) are not valid anymore
        //parent node (ie subNodeX) is not valid too
        parentNode->updateNodesCountRecursively();
        parentNode->invalidateRecursively( -1 );
        parentNode->updateIdenticalRecursively();
    }
    //otherwise, more traditional case
    else{
        //special treatment if the edge to move is linked to the root
        if (movedParent == root){
            /*                   movedChild->..
             * subNodeX<- root <
             *                   subNodeY->... ->insertParent->insertChild->..  */
            //detach the two other children of the root
            assert(root->getNumberChildren() == 3);
            InferenceNode* subNode1 = (InferenceNode*)(root->getChild(0));
            InferenceNode* subNode2 = (InferenceNode*)(root->getChild(1));
            if (movedChild==subNode1){
                subNode1 = (InferenceNode*)(root->getChild(2));
            }
            else{
                if (movedChild==subNode2){
                    subNode2 = (InferenceNode*)(root->getChild(2));
                }
            }
            oldDistance = subNode1->getParentDistance() +
                             subNode2->getParentDistance();
            root->removeChild( subNode1 );
            root->removeChild( subNode2 );
            //break the insert branch
            //root---->...---->insertParent ---/ /---->insertChild
            insertParent->removeChild( insertChild );
            //insertChild is a new child of the old root
            /*                                                         ->movedChild->..
             * ...<---subNodeX  // subNodeY->..->insertParent -> root <
             *                                                         ->insertChild->..  */
            root->addChild( insertChild, (1.0-frac)*distance );
            insertParent->addChild( root, frac*distance );
            //if we are moving the outgroup the old root must remain the root
            if (movedChild==outgroupNode){
                /*                                                        -> movedChild=outgroup
                 * ...<---subNodeX  // subNodeY<-..<-insertParent<-root <
                 *                                                        ->insertChild->..  */
                root->makeParent(NULL);
                InferenceNode* parentNode;
                InferenceNode* childNode;
                //subNodeY is the one which has been given a parent during root->makeParent(NULL);
                if ( subNode1->getParent() ){
                    parentNode = subNode1;
                    childNode = subNode2;
                }
                else{
                    parentNode = subNode2;
                    childNode = subNode1;
                }
                /*                                                       ->movedChild=outgroup
                 * ...<---subNodeX <- subNodeY<-..<-insertParent<-root <
                 *                                                       ->insertChild->..  */
                parentNode->addChild( childNode, oldDistance );
                parentNode->updateNodesCountRecursively();
                parentNode->invalidateRecursively(-1);
                parentNode->updateIdenticalRecursively();
            }
            //if the outgroup stays in place..
            else{
                //subNodeX is the outgroup
                /*                                                                    movedChild->..
                 * subNodeX<-newroot=subNodeY->..->insertParent->movedParent=oldRoot <
                 *                                                                    insertChild->..  */
                if (subNode1==outgroupNode){
                    subNode2->addChild( subNode1, oldDistance );
                    root = subNode2;
                }
                else{
                    assert(subNode2==outgroupNode);
                    subNode1->addChild( subNode2, oldDistance );
                    root = subNode1;
                }
                //the ancestor of the root cluster should be updated in installOutgroup
                //if (clustersTree){
                //    clustersTree->getRoot()->cluster->setAncestor(root);
                //}
                movedParent->updateNodesCountRecursively();
                movedParent->invalidateRecursively(-1);
                movedParent->updateIdenticalRecursively();
            }
        }
        //normal and simple case
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
            //                          movedChild
            //parentNode->movedParent <
            //                          childNode

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
            //unplug movedParent from the tree and join movedParent->parent and childNode
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
            parentNode->invalidateRecursively( -1 );
            parentNode->updateIdenticalRecursively();
        }
    }
    if (!installOutgroup()){
        createIndex();
    }
    return distance/oldDistance;
}

void UnrootedMCMCTree::getAllPerturbationParameters( vector<double>& params ) const{
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
    assert( params.size() == getNumberPerturbationParameters() );
}

void UnrootedMCMCTree::setAllPerturbationParameters( const vector<double>& params){
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
    assert(iter==params.end());
}

unsigned int UnrootedMCMCTree::getNumberPerturbationParameters() const{
    unsigned int numberPerturbationParameters = 6 + MCMCTreeBasic::getNumberPerturbationParameters();
    numberPerturbationParameters += branchPerturbator->getNumberPerturbationParameters();
    return numberPerturbationParameters;
}

void UnrootedMCMCTree::getAllPriorParameters( vector<double>& params ) const{
    MCMCTreeBasic::getAllPriorParameters(params);
    vector<double> p;
    branchPerturbator->getAllPriorParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
}

void UnrootedMCMCTree::setAllPriorParameters( const vector<double>& params ){
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

unsigned int UnrootedMCMCTree::getNumberPriorParameters() const{
    unsigned int numberPriorParameters = MCMCTreeBasic::getNumberPriorParameters();
    numberPriorParameters += branchPerturbator->getNumberPriorParameters();
    return numberPriorParameters;
}

double UnrootedMCMCTree::getLnPrior() const{
    double lnPrior = MCMCTreeBasic::getLnPrior();
    for ( unsigned int branchId = 0; branchId < getNumberBranches(); ++branchId ){
        lnPrior += branchPerturbator->getLnPriorValue(
                getBranchLength( branchId ) );
    }
    return lnPrior;
}


