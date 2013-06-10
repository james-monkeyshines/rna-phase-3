#include "Tree/SearchTreeInsertionBasic.h"

#include "Models/Model.h"

#include "Tree/OptimizerNode.h"

#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

#include <float.h>
#include <algorithm>

//define the DEBUG_MESSAGE flag to turn on information messages and follow
//the movements of clusters, tempoutgroup, rootclusters, ...

using namespace std;

template <class Pair>
struct myselect1st : public unary_function< Pair, typename Pair::first_type > {
    const typename Pair::first_type& operator()(const Pair& x) const {
        return x.first;
    }
};


SearchTreeInsertionBasic::SearchTreeInsertionBasic( const string & registrationName ) :
SearchTreeBasic(registrationName){
}


SearchTreeInsertionBasic::SearchTreeInsertionBasic(  ParametersSet& treeParameters ) :
InferenceTree(),OptimizerTree(),SearchTreeBasic(treeParameters){
    root = new OptimizerNode(this);
    createIndex();
    tempOutgroup = NULL;
    waitingClusters.clear();
    rootClusters.clear();
}

bool SearchTreeInsertionBasic::initialisation( SequenceTable * ptable, Model * pmodel ){
    return SearchTreeBasic::initialisation( ptable, pmodel );
}

void SearchTreeInsertionBasic::createInitialThreeLeavesTree( unsigned int spec1, unsigned int spec2, unsigned int spec3 ){

    Singleton< randombox > & randBox = Singleton< randombox >::instance();

    OptimizerNode * node = new OptimizerNode( ptable->species[spec1], this );
    //stack up the clusters for this node
    if (clustersTree){
        ClustersTreeNode* insertionPoint = clustersTree->findInsertion( spec1 );
        while ( insertionPoint ){
            waitingClusters[insertionPoint]=node;
            clustersQueue[node].push_back(insertionPoint);
            insertionPoint->cluster->setAncestor(node);
#ifdef DEBUG_MESSAGE
            cout << "NB: " << insertionPoint->cluster->getName() << " first insertion, waiting at " << node << ':' << node->toString() << endl;
#endif
            insertionPoint=insertionPoint->parent;
        }
    }
    cout << "Adding: " << ptable->species[spec1] << endl;
    assert(root);
    add((OptimizerNode*)node,(OptimizerNode*)root,randBox.ran()*.6);

    //creating the tree with the two last species
    for (unsigned int spec = 1; spec < 3; ++spec){
        //specId = spec2 if spec==1; specId = spec3 if spec==2
        unsigned int specId;
        specId = (spec == 1) ? spec2 : spec3;
        //create the node for the species
        cout << "Adding: " << ptable->species[specId] << endl;
        OptimizerNode * node = new OptimizerNode( ptable->species[specId], this );
        if (clustersTree){
            //look for the smallest cluster the species might be in
            ClustersTreeNode* insertionPoint = clustersTree->findInsertion( specId );
            assert(insertionPoint);
            //stack up the new clusters until a cluster already declared is found (root must have been declared)
            while ( waitingClusters.find(insertionPoint)==waitingClusters.end() ){
                assert(insertionPoint);
                waitingClusters[insertionPoint]=node;
                //temporarily set the ancestor node
                insertionPoint->cluster->setAncestor(node);
                clustersQueue[node].push_back(insertionPoint);
#ifdef DEBUG_MESSAGE
                cout << "NB: " << insertionPoint->cluster->getName() << " first insertion, waiting at "
                     << node << ':' << node->toString() << endl;
#endif
                insertionPoint=insertionPoint->parent;
            }
            //find the node in the cluster insertionPoint
            OptimizerNode* oldNode = waitingClusters[insertionPoint];
#ifdef DEBUG_MESSAGE
            cout << "NB: " << insertionPoint->cluster->getName() << " detected at " << oldNode << ':' << oldNode->toString() << endl;
#endif
            assert(insertionPoint->cluster->getAncestor() == oldNode);
            if (oldNode==root){
                assert(spec==2);
                //insertionPoint in the root, prepare rootClusters
                deque<ClustersTreeNode*>::iterator iter = clustersQueue[(OptimizerNode*)root].begin();
                assert(rootClusters.empty());
                while ( *iter != insertionPoint ){
#ifdef DEBUG_MESSAGE
                    cout << "NB: " << insertionPoint->cluster->getName() << " added to rootClusters " << endl;
#endif
                    rootClusters.push_back(*iter);
                    ++iter;
                }
                //this node is selected as a temporary outgroup if the two others are in separate clusters
                if (!rootClusters.empty()){
                    tempOutgroup = node;
#ifdef DEBUG_MESSAGE
                    cout << "NB: " << tempOutgroup << ':' << tempOutgroup->toString() << " declared tempOutgroup" << endl;
#endif
                }
            }
            else{
                deque<ClustersTreeNode*>& oldStack = clustersQueue[oldNode];
                deque<ClustersTreeNode*>& newStack = clustersQueue[(OptimizerNode*)root];
                do{
#ifdef DEBUG_MESSAGE
                    cout << "NB: " << oldStack.back()->cluster->getName() << " waiting at " << oldNode << ':'
                          << oldNode->toString() << " sent to root" << endl;
#endif
                    newStack.push_front(oldStack.back());
                    waitingClusters[oldStack.back()]=(OptimizerNode*)root;
                    //if we are adding the third node and he belongs to a cluster with one
                    //and only one of the two first nodes then we can set up a temporary outgroup too
                    if(spec==2){
                        rootClusters.push_front(oldStack.back());
#ifdef DEBUG_MESSAGE
                      cout << "NB: " << ".... also sent to rootClusters" << endl;
#endif
                    }
                    oldStack.back()->cluster->setAncestor((OptimizerNode*)root);
                    oldStack.pop_back();
                } while( newStack.front() != insertionPoint );
                //initialise the "temporary outgroup" (ie one of the two first nodes)
                if(spec==2){
                    assert(root->getNumberChildren()==2);
                    tempOutgroup = root->getChild(0);
                    if (tempOutgroup==oldNode){
                        tempOutgroup = root->getChild(1);
                    }
#ifdef DEBUG_MESSAGE
                    cout << "NB: " << tempOutgroup << ':' << tempOutgroup->toString() << " declared outgroup" << endl;
#endif
                }
            }
        }
        add((OptimizerNode*)node,(OptimizerNode*)root,randBox.ran()*.6);
    }
}


void SearchTreeInsertionBasic::retrievePossibleBranches( vector< InferenceNode* >& possibleBranches,
                                   unsigned int speciesId ){
    possibleBranches.clear();
    vector< pair< InferenceNode*, double > > temp;
    temp.clear();

    if (clustersTree){
        ClustersTreeNode* insertionPoint = clustersTree->findInsertion( speciesId );
        while ( waitingClusters.find(insertionPoint)==waitingClusters.end() ){
            insertionPoint=insertionPoint->parent;
        }
        //the root of the cluster tree must have been found, insertionPoint should not be NULL
        assert(insertionPoint);
        OptimizerNode * oldNode = waitingClusters[insertionPoint];
        deque<ClustersTreeNode*>& oldStack = clustersQueue[oldNode];
#ifdef DEBUG_MESSAGE
        cout << "NB: " << insertionPoint->cluster->getName() << " detected at " << oldNode << ':' << oldNode->toString() << endl;
#endif
        if (oldNode!=root){
            //general case: we have to insert ourselves under a cluster
            assert(oldNode->getParent());
            //the cluster parent branch is allowed, we will have to "split"
            //the oldNode in two nodes parent and child and add it as the
            //second children of the new parent if this branch is chosen
            possibleBranches.push_back( (InferenceNode*)oldNode );
            // we just have to find the allowed branches under this one
            // if oldNode does not define some clusters which are "smaller"
            if(oldStack.front()==insertionPoint){
#ifdef DEBUG_MESSAGE
                cout << "NB: insertion under cluster " << insertionPoint->cluster->getName() << endl;
#endif
                assert(insertionPoint->cluster->getAncestor()==oldNode);
                findClusterBranches( temp, insertionPoint );
            }
        }
        //oldNode == root
        else{
#ifdef DEBUG_MESSAGE
            cout << "NB: insertion at cluster " << insertionPoint->cluster->getName() << " waiting at the root" << endl;
#endif
            /* if insertion point is at the root we have two special cases since
             * one of the branch may or may not be allowed in the cluster */
            //rootClusters is used to record the clusters for the 2 selected children
            if ( find(rootClusters.begin(),rootClusters.end(),insertionPoint) != rootClusters.end() ){
                assert(tempOutgroup);
#ifdef DEBUG_MESSAGE
                cout << "NB: ....and found in rootClusters ROOT ->branch to tempoutgroup possible "
                     << tempOutgroup << ':' << tempOutgroup->toString() << endl;
#endif
                possibleBranches.push_back( (InferenceNode*)tempOutgroup );
                //we can insert ourselves under rootClusters
                if(rootClusters.front()==insertionPoint){
                    assert(rootClusters.front()==oldStack.front());
#ifdef DEBUG_MESSAGE
                    cout << "NB: ... and found 1st in rootClusters" << endl;
#endif
                    //one branch is not allowed, this cluster is only possible for
                    //two children of the root (the third child must be tempOutgroup)
                    bool out = false;
                    for ( list< BasicNode* >::iterator iter = root->getChildrenList().begin();
                        iter != root->getChildrenList().end(); ++iter ){
                        if ( *iter != tempOutgroup ){
                            ((OptimizerNode*)(*iter))->findIncludedBranches( temp, insertionPoint );
                        }
                        else{
                            assert(!out);
                            out=true;
                        }
                    }
                    assert(out);
                }
            }
            else{
                //if the insertion point is not in root clusters it might call for a rerooting later
                if(!rootClusters.empty()){
                    assert(tempOutgroup);
#ifdef DEBUG_MESSAGE
                    cout << "NB: ... and NOT found in rootClusters" << endl;
#endif
                    //we check whether insertionPoint is the smallest cluster not in rootCluster waiting at the root
                    if ( rootClusters.back()->parent == insertionPoint ){
#ifdef DEBUG_MESSAGE
                        cout << "NB: insertion above/under " << tempOutgroup << ':' << tempOutgroup->toString() << endl;
#endif
                        ((OptimizerNode*)(tempOutgroup))->findIncludedBranches( temp, insertionPoint );
                    }
                    else{
#ifdef DEBUG_MESSAGE
                        cout << "NB: insertion above " << tempOutgroup << ':' << tempOutgroup->toString()
                             << " waiting for rerooting..." << endl;
#endif
                        possibleBranches.push_back( (InferenceNode*)tempOutgroup );
                    }
                }
                else{
#ifdef DEBUG_MESSAGE
                    cout << "NB: no outgroup/rootClusters" << endl;
#endif
                    assert(!tempOutgroup);
                    //we check whether insertionPoint is the smallest cluster waiting at the root
                    if ( clustersQueue[(OptimizerNode*)root].front() == insertionPoint ){
                        ((OptimizerNode*)(root))->findClusterBranches( temp, insertionPoint );
                    }
                    else{
#ifdef DEBUG_MESSAGE
                        cout << "NB: insertion at root children waiting for rerooting" << endl;
#endif
                        for ( list< BasicNode* >::iterator iter = root->getChildrenList().begin();
                            iter != root->getChildrenList().end(); ++iter ){
                                possibleBranches.push_back( (InferenceNode*)(*iter) );
                        }
                    }
                }
            }
        }
    }
    else{
        ((OptimizerNode*)(root))->findClusterBranches( temp, NULL );
    }

    //push the branches in temp to possibleBranches
    unsigned int output = possibleBranches.size();
    possibleBranches.resize( output + temp.size() );
    vector< InferenceNode* >::iterator iter = possibleBranches.begin();
    advance( iter, output);
    transform( temp.begin(), temp.end(), iter, myselect1st< pair< InferenceNode*, double > >() );
}


void SearchTreeInsertionBasic::insertion( OptimizerNode * node, OptimizerNode * newParent, OptimizerNode * insertNode ){
#ifdef DEBUG_MESSAGE
    cout << "NB: insertion" << endl;
#endif
    double dist = insertNode->getParentDistance();
    //for "security", if dist == 0.00 (it can happen if insertNode is the intermediate
    //state between its parent and its relative) then likelihood might start with -inf.
    dist = MAX(dist,.006);
    //insert new parent between allowedBranches[i].first and its parent
    OptimizerNode* oldParent = (OptimizerNode*)(insertNode->getParent());
    remove(insertNode);
    add(insertNode, newParent, dist/2.0 );
    //we add the new species as the second child of newParent using the same branch length
    add(node, newParent, dist/2.0 );
    //plug back the subtree to the tree
    add(newParent, oldParent, dist/2.0 );
}

void SearchTreeInsertionBasic::removal( OptimizerNode * node, OptimizerNode * newParent, OptimizerNode * insertNode ){
#ifdef DEBUG_MESSAGE
    cout << "NB: removal" << endl;
#endif
    double dist = 2.0*insertNode->getParentDistance();
    //insert new parent between allowedBranches[i].first and its parent
    OptimizerNode* oldParent = (OptimizerNode*)(newParent->getParent());
    assert(node->getParent()== newParent);
    assert(insertNode->getParent()== newParent);
    //tempOutgroup should have been put back in removalValidate if it was necessary
    assert(tempOutgroup!=newParent);
    remove(newParent);
    remove(insertNode);
    remove(node);
    add(insertNode, oldParent, dist );
}

void SearchTreeInsertionBasic::insertionValidate( OptimizerNode * node, OptimizerNode * newParent, OptimizerNode * insertNode ){

    assert(node->isLeaf());
    assert(node->cnumber!=-1);
    assert(rootSave.size()==(getNumberBranches()-5)/2);

    rootSave.push( pair< OptimizerNode*, OptimizerNode* >(node, (OptimizerNode*)root) );
#ifdef DEBUG_MESSAGE
    cout << "NB: PUSHING ROOT " << root << endl;
#endif
    //update the state of the tree (clusters)
    if (clustersTree){
        assert(node->isLeaf());
        assert(node->cnumber!=-1);
        ClustersTreeNode* insertionPoint = clustersTree->findInsertion( node->cnumber - 1 );
        assert(insertionPoint);
#ifdef DEBUG_MESSAGE
        cout << "NB: insertionValidate for " << node->getLabel() << ", id=" << node->cnumber
             << " from " << insertionPoint->cluster->getName() << endl;
#endif
        //assign the new clusters to this leaf
        while ( waitingClusters.find(insertionPoint)==waitingClusters.end() ){
            assert(insertionPoint);
            waitingClusters[insertionPoint]=node;
            //temporarily set the ancestor node
            insertionPoint->cluster->setAncestor(node);
            clustersQueue[node].push_back(insertionPoint);
#ifdef DEBUG_MESSAGE
            cout << "NB: " << insertionPoint->cluster->getName() << " first insertion, waiting at " << node << ':' << node->toString() << endl;
#endif
            insertionPoint=insertionPoint->parent;
        }
        assert(insertionPoint);
        OptimizerNode * oldNode = waitingClusters[insertionPoint];
        if (oldNode!=root){
            assert(oldNode->getParent());
            //general case: we inserted the node under a cluster (can be root cluster)
            //and a split occured (when insert at the parent branch of the cluster
            //ancestor
            if(insertNode==oldNode){
                // oldNode was "splitted in two to accomodate the new species or was
                // replaced as a cluster ancestor.
                // we must transmit the clusters bigger or equal than insertionPoint
                // to newParent
                deque<ClustersTreeNode*>& oldStack = clustersQueue[oldNode];
                deque<ClustersTreeNode*>& newStack = clustersQueue[newParent];
                while( newStack.front() != insertionPoint ){
                    newStack.push_front(oldStack.back());
                    oldStack.back()->cluster->setAncestor(newParent);
                    waitingClusters[oldStack.back()]=newParent;
                    oldStack.pop_back();
#ifdef DEBUG_MESSAGE
                    cout << "NB: " << newStack.front()->cluster->getName() << " waiting at " << insertNode << ':' << oldNode->toString()
                                                                           << " sent to new parent " << newParent << ':'
                                                                           << newParent->toString() << endl;
#endif
                }
                if (insertNode==tempOutgroup){
#ifdef DEBUG_MESSAGE
                    cout << "NB: ...tempOutgroup follow" << endl;
#endif
                    tempOutgroup = newParent;
                }
            }
        }
        // oldNode == root, insertionPoint at the root...
        else{
            //insertionPoint in rootClusters...
            if ( find(rootClusters.begin(),rootClusters.end(),insertionPoint) != rootClusters.end() ){
                //rerooting?
                if (insertNode==tempOutgroup){
                    //we install the new root at that point
#ifdef DEBUG_MESSAGE
                    cout << "NB: rerooting on the tempOutgroup chosen" << endl;
#endif
                    OptimizerNode* oldRoot = (OptimizerNode*)root;
                    makeRoot( (InferenceNode*)newParent );
                    createIndex();
                    assert(root==newParent);
                    reroot( oldRoot, newParent, insertionPoint, node );
                }
            }
            else{
                //insertionPoint NOT in rootClusters...
                if(!rootClusters.empty()){
                    assert(tempOutgroup);
                    if ( rootClusters.back()->parent == insertionPoint ){
                        //change tempOutgroup if necessary
                        if (insertNode==tempOutgroup){
                            tempOutgroup = newParent;
#ifdef DEBUG_MESSAGE
                            cout << "NB: insertion in the outgroup branch, tempOutgroup moved to " << newParent
                                 << ':' << newParent->toString() << endl;
#endif
                         }
                    }
                    else{
                        OptimizerNode* oldRoot = (OptimizerNode*)root;
                        makeRoot( (InferenceNode*)newParent );
                        createIndex();
                        assert(root==newParent);
                        reroot( oldRoot, newParent, insertionPoint, node );
                    }
                }
                else{
                    assert(!tempOutgroup);
                    if ( clustersQueue[(OptimizerNode*)root].front() == insertionPoint ){
                        //nothing special, this is a traditional insertion
                    }
                    else{
                        OptimizerNode* oldRoot = (OptimizerNode*)root;
                        makeRoot( (InferenceNode*)newParent );
                        createIndex();
                        assert(root==newParent);
                        reroot( oldRoot, newParent, insertionPoint, node );
                    }
                }
            }
        }
    }
    //if no clusters there is nothing special to do


    //update the state of the tree (outgroup and rerooting)
    OptimizerNode* outSpecies = (OptimizerNode*)findOutgroup();
    //if the outgroup is present in the tree
    if (outSpecies != NULL){
        OptimizerNode* oldRoot = (OptimizerNode*)root;
        //(in the case of a cluster) if the outgroup
        //is at the root special root case with clusters
        if (outSpecies == root){
            assert(clustersTree);
            //if tempOutgroup exists and can receive the root
            if ( tempOutgroup && (!tempOutgroup->isLeaf()) ){
                assert(rootClusters.back()->cluster==outgroupCluster);
#ifdef DEBUG_MESSAGE
                cout << "NB: outgroup species differentiated... rerooting" << endl;
#endif
                //put an iterator in clustersQueue[oldRoot] to point at the first clusters not in rootClusters
                //(ie the smallest clusters which contains ALL the species currently in the tree)
                deque<ClustersTreeNode*>::iterator iter = find( clustersQueue[(OptimizerNode*)root].begin(),
                                                                clustersQueue[(OptimizerNode*)root].end(),
                                                                rootClusters.back());
                ++iter;
                //at least the all-species cluster must be there
                assert(iter!=clustersQueue[(OptimizerNode*)root].end());
                //change rootClusters, the future root clusters will be all the clusters currently waiting
                //at tempOutgroup (which is, in some sense, the opposite of the outgroup we want)
                rootClusters=clustersQueue[(OptimizerNode*)tempOutgroup];
#ifdef DEBUG_MESSAGE
                cout << "NB: root clusters destroyed and recontructed with node waiting at false outgroup" << endl;
#endif
                //the clusters waiting at oldRoot and containing everything are sent to the
                //future new root (ie tempOutgroup)
                deque<ClustersTreeNode*>::iterator iter2 = iter;
                while( iter2 != clustersQueue[(OptimizerNode*)root].end() ){
                    waitingClusters[*iter2]=(OptimizerNode*)tempOutgroup;
                    (*iter2)->cluster->setAncestor(tempOutgroup);
                    clustersQueue[(OptimizerNode*)tempOutgroup].push_back( *iter2 );
                    ++iter2;
                }
#ifdef DEBUG_MESSAGE
                cout << "NB: clusters comprising all the species currently in the tree were sent to the future root" << endl;
#endif
                //remove clusters bigger than outgroup clusters in the old root (future tempOutgroup and outgroup)
                clustersQueue[(OptimizerNode*)root].erase(iter, clustersQueue[(OptimizerNode*)root].end());
                //swap root
                makeRoot( (InferenceNode*)tempOutgroup );
                createIndex();
                //new tempoutgroup
                if (rootClusters.empty()){
#ifdef DEBUG_MESSAGE
                    cout << "NB: temp outgroup destroyed" << endl;
#endif
                    tempOutgroup = NULL;
                }
                else{
#ifdef DEBUG_MESSAGE
                    cout << "NB: temp outgroup moved" << endl;
#endif
                    tempOutgroup = outSpecies;
                }
            }
        }
        else{
            if (installOutgroup()){
#ifdef DEBUG_MESSAGE
                cout << "NB: outgroup species added... rerooting: " << toString() << endl;
#endif
                if(clustersTree){
                    //standard case, all the global clusters are moved to the new root
                    deque<ClustersTreeNode*>::iterator iter = clustersQueue[oldRoot].begin();
                    if(!rootClusters.empty()){
                        while ( *iter != rootClusters.back() ){
                            ++iter;
                        }
                        ++iter;
                    }
                    reroot( oldRoot, (OptimizerNode*)root, *iter, outSpecies );
                }
            }
        }
    }

    //if there is an outgroup (!NULL and !root) then tempOutgroup = this outgroup or = null (if rootClusters is empty)
    assert( (outSpecies == NULL) || (outSpecies == root) || (outSpecies == tempOutgroup) || (tempOutgroup==NULL) );

}

void SearchTreeInsertionBasic::removalValidate( OptimizerNode * node, OptimizerNode * newParent, OptimizerNode * insertNode ){
    pair< OptimizerNode*, OptimizerNode* >  old = rootSave.top();
    rootSave.pop();
    assert(rootSave.size()==(getNumberBranches()-5)/2);
    OptimizerNode* oldRoot = old.second;
    assert( node == old.first );
    assert(node->isLeaf());
    assert(node->cnumber!=-1);


    //the tree was updated...
    //either we added a outgroup species (or included in the outgroup cluster)
    //either the outgroup cluster was added before but could not be used as tempoutgroup
    //because it was in rootClusters till the node was added
    //either the species we added splitted rootClusters in two
#ifdef DEBUG_MESSAGE
    cout << "NB: POPING ROOT " << oldRoot << endl;
#endif
    if(root!=oldRoot){
#ifdef DEBUG_MESSAGE
        cout << "NB: REROOTING " << toString() << endl;;
#endif
        //the root changed when this node was inserted...
        //it could be because the node which was inserted was the outgroup species
        //or a node which splitted the clusters waiting at the root
        //in any cases, newParent must have been elected to be the new root.
        assert(newParent==root);
        //since there was a rerooting, the previous tempOutgroup IF IT EXISTED
        //could only be the actual parent of the oldRoot
        OptimizerNode* oldOutgroup = (OptimizerNode*)(oldRoot->getParent());
        //we memorize whether the new node has been added as a tempOutgroup
        //(it imply the new node did not split the previous rootClusters but forced us to create a new one,
        //since it was added in-between the global clusters waiting at the root not in rootClusters)
        //Remember rootClusters are the clusters that apply to two node out of three at the root,
        //tempOutgroup being the root child excluded
        bool wastempOutgroup = (tempOutgroup==node);
        //reroot the tree
        makeRoot(oldRoot);
        createIndex();
        //if clusters are used, take care of restoring their state
        if (clustersTree){
            OptimizerNode* cancelledNewRoot = (OptimizerNode*)newParent;
            deque<ClustersTreeNode*>::iterator iter = clustersQueue[cancelledNewRoot].begin();
            //oldInsertionPoint can be the first element in rootClusters if the old rootClusters
            //was splitted. If rootClusters exist and was not splitted but overwritten
            //it means the oldInsertionPoint if the first element of the cancelledQueue not
            //in rootClusters. We can know this is the case if the node we added was declared tempOutgroup
            if( !rootClusters.empty() && wastempOutgroup ){
                while ( *iter != rootClusters.back() ){
                    ++iter;
                }
                ++iter;
            }
            cancelReroot( cancelledNewRoot, oldRoot, *iter, oldOutgroup );
        }
#ifdef DEBUG_MESSAGE
        cout << "NB: REROOTED TO " << toString() << endl;
        if (tempOutgroup) cout << "NB: tempOutgroup is " << tempOutgroup->toString() << endl;
#endif
    }

    if(clustersTree){
#ifdef DEBUG_MESSAGE
            cout << "NB: removalValidate for " << node->getLabel() << ", id=" << node->cnumber << endl;
#endif
        //node will disappear... removal is easy
        deque<ClustersTreeNode*>& nodeStack = clustersQueue[node];
        while( !nodeStack.empty() ){
            waitingClusters.erase(nodeStack.back());
            nodeStack.back()->cluster->setAncestor(NULL);
            nodeStack.pop_back();
        }
        clustersQueue.erase(node);

        //parent will disappear too, but it is harder to remove it...
        //...the simple case, there has been no rerooting during the process
        assert( insertNode->getParent() == newParent );
        assert( newParent != root );
        //removal "easy"
        deque<ClustersTreeNode*>& parentStack = clustersQueue[newParent];
        deque<ClustersTreeNode*>& insertStack = clustersQueue[insertNode];
        //sent back the cluster in newParent to the insertNode
        while( !parentStack.empty() ){
            waitingClusters[ parentStack.front() ] = insertNode;
            parentStack.front()->cluster->setAncestor( insertNode );
            insertStack.push_back( parentStack.front() );
#ifdef DEBUG_MESSAGE
            cout << "NB: send back " << parentStack.front()->cluster->getName() << " to " << insertNode->toString() << endl;
#endif
            parentStack.pop_front();
        }
        clustersQueue.erase(newParent);
        if (tempOutgroup==newParent){
            tempOutgroup = insertNode;
#ifdef DEBUG_MESSAGE
            cout << "NB: tempOutgroup sent back to " << insertNode->toString() << endl;
#endif
        }
    }
}



void SearchTreeInsertionBasic::reroot( OptimizerNode* oldRoot, OptimizerNode* newRoot, ClustersTreeNode* insertionPoint, OptimizerNode* outgroup ){
    //we had rootClusters < clustersQueue[root] and insertionPoint in clustersQueue[root] and maybe in rootClusters;
    deque<ClustersTreeNode*>& oldStack = clustersQueue[oldRoot];
    deque<ClustersTreeNode*>& newStack = clustersQueue[newRoot];
    
    //iter2 is an iterator on oldStack: the smallest cluster in oldStack to be sent to newStack
    deque<ClustersTreeNode*>::iterator iterBeginNew;
    
    //first reconstruct rootClusters for the new root:
    //if insertionPoint in rootClusters, (it means the root is moving to one of its child)
    //remove the clusters smaller than insertionPoint from rootClusters (they stay with oldRoot)
    deque<ClustersTreeNode*>::iterator iter = find(rootClusters.begin(),rootClusters.end(),insertionPoint);
    if ( iter != rootClusters.end() ){
        assert(oldRoot->getParent()==newRoot);
        //remove element from rootClusters, tempOutgroup remains the same
#ifdef DEBUG_MESSAGE
        cout << "NB: removing " << (*rootClusters.begin())->cluster->getName()
             << " to " << (*iter)->cluster->getName() << "(excluded) from rootClusters." << endl;
#endif
        rootClusters.erase( rootClusters.begin(), iter );
        //clusters from insertionPoint are sent to the new root
        iterBeginNew = find( oldStack.begin(), oldStack.end(), insertionPoint );
    }
    //if insertionPoint NOT in rootCluster it can be a long-distance change
    //insertionPoint is in the last part of oldStack, the range of clusters waiting at the root
    //that are not in rootClusters. The old Root will conserve all the clusters in rootClusters,
    //other clusters are sent to the new root
    else{
        iterBeginNew = oldStack.begin();
        if (!rootClusters.empty()){
            assert(tempOutgroup);
            while(*iterBeginNew!=rootClusters.back()){
                ++iterBeginNew;
                assert(iterBeginNew!=oldStack.end());
            }
            ++iterBeginNew;
            assert(iterBeginNew!=oldStack.end());
        }
        //rootClusters is reconstructed for the new root.
        //locate insertionPoint in oldStack and create the new rootClusters with clusters
        //smaller than insertionPoint that will not stay with the old root
        deque<ClustersTreeNode*>::iterator iterIns = find(iterBeginNew,oldStack.end(),insertionPoint);
#ifdef DEBUG_MESSAGE
        cout << "NB: overwrite rootClusters, with " << (*iterBeginNew)->cluster->getName()
             << " to " << (*iterIns)->cluster->getName() << "(excluded)." << endl;
#endif
        rootClusters.clear();
        rootClusters.insert( rootClusters.begin(), iterBeginNew, iterIns );
        if (!rootClusters.empty()){
#ifdef DEBUG_MESSAGE
            cout << "NB: the added species is new tempOutgroup" << endl;
#endif
            //new tempOutgroup
            tempOutgroup = outgroup;
        }
        else{
#ifdef DEBUG_MESSAGE
            cout << "NB: rootClusters empty, destroy tempOutgroup" << endl;
#endif
            tempOutgroup = NULL;
        }
    }
    
    //the clusters bigger than insertionPoint (case insertionPoint<rootClusters.back())
    //OR the clusters that were not in rootClusters (case insertionPoint not in rootClusters) are sent to the new root
    //at least the root cluster should be there
    for( deque<ClustersTreeNode*>::iterator iter = iterBeginNew; iter != oldStack.end(); ++iter ){
        waitingClusters[*iter]=(OptimizerNode*)newRoot;
        (*iter)->cluster->setAncestor(newRoot);
#ifdef DEBUG_MESSAGE
        cout << "NB: sending " << (*iter)->cluster->getName() << " from the old root to the new root" << endl;
#endif
    }
    newStack.insert( newStack.end(), iterBeginNew, oldStack.end() );
    oldStack.erase( iterBeginNew, oldStack.end() );
}

void SearchTreeInsertionBasic::cancelReroot( OptimizerNode* cancelledNewRoot, OptimizerNode* restoredRoot, ClustersTreeNode* oldInsertionPoint, OptimizerNode* outgroup ){
    //we have rootClusters < clustersQueue[cancelledNewRoot] and insertionPoint in clustersQueue[root] and maybe in rootClusters;
    deque<ClustersTreeNode*>& cancelledStack = clustersQueue[cancelledNewRoot];
    deque<ClustersTreeNode*>& restoredStack = clustersQueue[restoredRoot];
    
    assert( find(cancelledStack.begin(),cancelledStack.end(),oldInsertionPoint) != cancelledStack.end());
    
    //first reconstruct rootClusters for the restored root:
    deque<ClustersTreeNode*>::iterator iter = find(rootClusters.begin(),rootClusters.end(),oldInsertionPoint);
    //if oldInsertionPoint is still in rootClusters, (it means the root moved to one of its child)
    //send back the clusters of the restoredRoot in rootClusters
    if ( iter != rootClusters.end() ){
        assert(iter==rootClusters.begin());
        //the rerooting was already done, check that the restored root is already parent of the cancelled one
        assert(cancelledNewRoot->getParent()==restoredRoot);
        if (!restoredStack.empty()){
            rootClusters.insert( rootClusters.begin(), restoredStack.begin(), restoredStack.end() );
#ifdef DEBUG_MESSAGE
             cout << "NB: restoring rootClusters " << rootClusters.front()->cluster->getName()
                  << " to " << rootClusters.back()->cluster->getName() << endl;
#endif
        }
        else{
            //nothing this is the case were all the elements form the old rootclusters have been
            //used for the new root clusters (ie the node added inserted itself in the branch leading
            //to temp outgroup but still was an element of the first cluster in rootClusters
        }
    }
    //if insertionPoint was NOT in rootCluster it could have been a long-distance change
    //oldInsertionPoint is the first cluster waiting at the cancelled root
    //that is not in not in rootClusters
    //the cancelled Root will conserve all the clusters in rootClusters
    else{
        //rootClusters is reconstructed for the new root.
        if (!restoredStack.empty()){
#ifdef DEBUG_MESSAGE
            cout << "NB: overwrite rootClusters, with " << restoredStack.front()->cluster->getName()
                 << " to " << restoredStack.back()->cluster->getName() << "(included); tempoutgroup was " << outgroup->toString() << endl;
#endif
            rootClusters = restoredStack;
            tempOutgroup = outgroup;
        }
        else{
#ifdef DEBUG_MESSAGE
            cout << "NB: overwrite rootClusters (old rootClusters did not exists, clear empty and destroy tempOutgroup)" << endl;
#endif
            rootClusters.clear();
            tempOutgroup = NULL;
        }
    }
    
    //send the clusters back to the restoredStack
    for( deque<ClustersTreeNode*>::iterator iter = cancelledStack.begin(); iter != cancelledStack.end(); ++iter ){
        waitingClusters[*iter]=(OptimizerNode*)restoredRoot;
        (*iter)->cluster->setAncestor(restoredRoot);
#ifdef DEBUG_MESSAGE
        cout << "NB: sending " << (*iter)->cluster->getName() << " from the cancelled root to the restored root" << endl;
#endif
    }
    restoredStack.insert( restoredStack.end(), cancelledStack.begin(), cancelledStack.end() );
    cancelledStack.clear();
}

