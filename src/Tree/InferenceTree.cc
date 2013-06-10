#include "Tree/InferenceTree.h"

#include <math.h>
#include <assert.h>
#include <time.h>

#include <string>
#include <iostream>



#include "Util/statlib.h"

#include "Util/array2D.h"
#include "Util/array_bag.h"

#include "PatternDesign/Singleton.h"
#include "Util/randombox.h"

#include "Tree/InferenceNode.h"
#include "Tree/TreeMap.h"
#include "Tree/ClustersTree.h"
#include "Tree/ClustersTreeNode.h"

#include "Models/Model.h"


class SequenceTable;


using namespace std;



InferenceTree::InferenceTree() : Tree() {
    pmodel = NULL;
    ptable = NULL;
    clustersStream = NULL;
    clustersTree = NULL;
    root = NULL;
}


// Creates a tree from its string representation
InferenceTree::InferenceTree( string stree ) : Tree() {
    pmodel = NULL;
    ptable = NULL;
    clustersStream = NULL;
    clustersTree = NULL;
    // Create a depth map of the tree
    TreeMap treeMap( stree );
    // Create the nodes
    root = new InferenceNode( treeMap, this );
    createIndex();
}

// Constructs a tree by copying another tree
InferenceTree::InferenceTree( InferenceTree & src ) : Tree() {

    pmodel = NULL;
    ptable = NULL;
    clustersStream = NULL;
    clustersTree = NULL;

    //outgroup = src.outgroup;
    // since nodes have a reference to the tree the copy will be different
    root = new InferenceNode( this );
    if ( src.getTable() ) {
        assert( src.getModel() );
        // same model and sequence table (can be NULL)
        loadDataAndModel( src.getTable(), src.getModel() );
        *clustersTree = *(src.clustersTree);
    }
    ((InferenceNode*)root)->recursiveCopy( (InferenceNode*)(src.root), NULL );
    createIndex();
}

void InferenceTree::lookupResize(){
    assert(root);
    assert ( this->ptable != NULL );
    assert ( this->pmodel != NULL );
    unsigned int multiTree = ((InferenceNode*)root)->lookupResize();
    //resize the tree lookup array to make it suitable for all
    //internal nodes
    lookup.resize(multiTree);
    for ( unsigned int childId = 0; childId < multiTree; ++childId ) {
        lookup[childId] = new array3D<double> [ptable->getNumberCategories()];
        for( unsigned int cat = 0; cat < ptable->getNumberCategories(); ++cat ){
            lookup[childId][cat].resize(
                pmodel->getNumberRatesCategories(cat),
                pmodel->getNumberStates(cat),
                pmodel->getNumberStates(cat) );
        }
    }
}


InferenceTree & InferenceTree::quickCopy( const InferenceTree & src, bool withPartialLikelihood ){
        assert( ( src.pmodel == pmodel ) && ( src.ptable == ptable ) &&
         ( pmodel ) && ( ptable ) );
        assert( nodeRefVector.size() == src.nodeRefVector.size() );
        assert( *clustersTree == *(src.clustersTree) );
        //nodes are copied in reverse order (leaves -> root)
        vector<BasicNode*>::const_reverse_iterator iterSrc =
                src.nodeRefVector.rbegin();
        vector<BasicNode*>::reverse_iterator iterDest =
                nodeRefVector.rbegin();
        //to perform a perfect copy of the tree we have to restore the links
        //between node manually, we store the children visited with
        //the copy we made of them in the following map
        map< const BasicNode*, BasicNode* > mapNode;

        //loop invariant : nodes from rbegin to iterDest are made "equal"
        //to the nodes from src.rbegin() to iterSrc
        while ( iterSrc != src.nodeRefVector.rend() ){
            // Efficient rearrangement of the leaves, a simple swap is enough
            // since the partialLikelihood vector is invariant for them
            if ( (*iterSrc)->isLeaf() ){
                //look for the same leaf in the nodeRefVector to swap it
                vector<BasicNode*>::reverse_iterator iter = iterDest;
                while ( ((InferenceNode*)(*iter))->cnumber !=
                               ((InferenceNode*)(*iterSrc))->cnumber ) {
                    ++iter;
                    assert( iter != nodeRefVector.rend() );
                } // end while
                // swap if necessary
                if ( iter != iterDest ) {
                    BasicNode* temp = *iter;
                    *iter = *iterDest;
                    *iterDest = temp;
                }
#ifdef DEBUG2
                assert ( ((InferenceNode*)(*iterDest))->partialLikelihood ==
                         ((InferenceNode*)(*iterSrc))->partialLikelihood );
                assert( ((InferenceNode*)(*iterDest))->identical ==
                        ((InferenceNode*)(*iterSrc))->identical );
#endif
                //store the node in mapNode for future link of the new copy
                //to its father copied later in the process
                mapNode[*iterSrc] = *iterDest;
            }
            else{
                //look for the next internal node in the nodeRefVector to
                //swap it
                vector<BasicNode*>::reverse_iterator iter = iterDest;
                while ( ((InferenceNode*)(*iter))->isLeaf() ){
                    ++iter;
                    assert( iter != nodeRefVector.rend() );
                } // end while
                // swap if necessary
                if ( iter != iterDest ) {
                    BasicNode* temp = *iter;
                    *iter = *iterDest;
                    *iterDest = temp;
                }
                //copy calculated, likelihood, pnodeModel and identical
                //nbInternalNodes, nbInternalLeaves, actual children list
                //are cleaned
                ((InferenceNode*)(*iterDest))->quickCopy
                    (*((InferenceNode*)(*iterSrc)), withPartialLikelihood);
                //add the previously duplicated children now
                for ( list< BasicNode * >::const_iterator iter =
                           ((*iterSrc)->getChildrenList()).begin();
                      iter != ((*iterSrc)->getChildrenList()).end(); ++iter ){
                    //find the src child in the map
                    BasicNode* child = mapNode[*iter];
                    //link the child copy to the node
                    assert(child != NULL);
                    (*iterDest)->addChild(child,(*iter)->getParentDistance());
                }
                //store the node in mapNode for future link of the new copy
                //to its father copied later in the process
                mapNode[*iterSrc] = *iterDest;
            }
            ++iterSrc;
            ++iterDest;
        }
        root = nodeRefVector[0];
        return ( * this );
}


InferenceTree & InferenceTree::operator = ( const InferenceTree & src ) {

    // No self assignment
    if ( root == src.getRoot() )
        return ( * this );

    //the = operator is designed for trees dealing with the same set of species
    //with the same number of leaves and internal nodes
    if ( ( src.pmodel == pmodel ) && ( src.ptable == ptable ) &&
         ( pmodel ) && ( ptable ) ) {
        return quickCopy(src, true);
    }
    //general case, a full copy is required
    assert( root == NULL );
    root = new InferenceNode( this );
    if ( src.getTable() ) {
        assert( src.getModel() );
        // same model and sequence table (can be NULL)
        loadDataAndModel( src.getTable(), src.getModel() );
        *clustersTree = *(src.clustersTree);
    }
    ((InferenceNode*)root)->recursiveCopy( (InferenceNode*)(src.root), NULL );
    createIndex();
    return ( * this );
}

// Destructor
InferenceTree::~InferenceTree() {
    if (!lookup.empty()){
        for ( unsigned int childId = 0; childId < lookup.size(); ++childId ) {
            if (lookup[childId] != NULL){
                delete [] lookup[childId];
            }
        }
    }
}


void InferenceTree::loadDataAndModel( SequenceTable * psequenceTable,
Model * pmodel ) {
    assert ( this->ptable == NULL );
    assert ( this->pmodel == NULL );
    this->ptable = psequenceTable;
    this->pmodel = pmodel;

    if (ptable->getNumberCategories() != pmodel->getNumberSymbolCategory()){
        cerr << "ERROR: " << ptable->getNumberCategories()
             << " categories of symbols have been found in your sequence table"
             << " whereas the substitution model defined can deal with "
             << pmodel->getNumberSymbolCategory() << " categories" << endl;
        exit(EXIT_FAILURE);
    }
    invalidateNodes(-1);
    //create lookup table
    assert( lookup.empty() );
    lookup.resize(NUMBER_LOOKUP_TABLE);
    for ( unsigned int childId = 0; childId < NUMBER_LOOKUP_TABLE; ++childId ) {
        lookup[childId] = new array3D<double> [ptable->getNumberCategories()];
        for( unsigned int cat = 0;
             cat < ptable->getNumberCategories(); ++cat ){
            lookup[childId][cat].resize(
                pmodel->getNumberRatesCategories(cat),
                pmodel->getNumberStates(cat),
                pmodel->getNumberStates(cat) );
        }
    }

    //WARNING, the order to compute the identical is IMPORTANT
    //WARNING load data and model assume that nodes without children are leaves
    //this might cause some problem when this function is called on a partially
    //constructed tree (tree with just the root node accepted)
    for ( vector< BasicNode*>::reverse_iterator iter = nodeRefVector.rbegin();
          iter != nodeRefVector.rend(); ++iter ){
        if ( ((*iter)==root) || (*iter)->getNumberChildren() ){
            ((InferenceNode*)(*iter))->initialisation( false );
            //initialisation will compute the identical by itself if possible
            //((InferenceNode*)(*iter))->computeIdentical();
        }
        else{
            ((InferenceNode*)(*iter))->initialisation( true );
        }
    }
    if (clustersStream){
        processClusters();
    }
}

bool InferenceTree::makeRoot( InferenceNode* newRoot ){
    InferenceNode* oldRoot = (InferenceNode*)root;
    if ( Tree::makeRoot( newRoot ) ){
        if (pmodel){
            oldRoot->invalidateRecursively(-1);
            oldRoot->updateIdenticalRecursively();
            return true;
        }
    }
    return false;
}


string InferenceTree::toStringNumbered( bool topologyOnly ){
    assert( pmodel && ptable );
    return ( ((InferenceNode*)root)->toStringNumbered(topologyOnly) );
}



void InferenceTree::move( InferenceNode* src, InferenceNode* dest ){
    InferenceNode* parent = (InferenceNode*)(src->getParent());
    if (parent != dest){
        Tree::move( src, dest);
        if(pmodel){
            InferenceNode* ancestor =
                   (InferenceNode*)getMostRecentAncestor(parent, dest);
            parent->updateIdenticalRecursively(ancestor);
            parent->invalidateRecursively(-1,ancestor);
            dest->updateIdenticalRecursively();
            dest->invalidateRecursively(-1);
        }
    }
}

void InferenceTree::swap( InferenceNode* node1, InferenceNode* node2 ){
    InferenceNode* parent1 = (InferenceNode*)(node1->getParent());
    InferenceNode* parent2 = (InferenceNode*)(node2->getParent());
    //cannot swap the root
    assert( parent1 );
    assert( parent2 );
    assert( parent1 != parent2 );
    InferenceNode* ancestor = (InferenceNode*)getMostRecentAncestor(node1, node2);
    assert( (ancestor!=node1) && (ancestor!=node2) );
    if (parent1 != parent2){
        double distance1 = node1->getParentDistance();
        double distance2 = node2->getParentDistance();
        parent1->removeChild( node1 );
        parent2->removeChild( node2 );
        parent1->addChild( node2, distance2 );
        parent2->addChild( node1, distance1 );
    }
    if(pmodel){
        //update parameters recursively in the root direction from parent1
        //stop just before ancestor
        parent1->updateNodesCountRecursively(ancestor);
        parent1->updateIdenticalRecursively(ancestor);
        parent1->invalidateRecursively(-1,ancestor);
        //update parameters recursively until the root is reached
        //(cross ancestor)
        parent2->updateNodesCountRecursively();
        parent2->updateIdenticalRecursively();
        parent2->invalidateRecursively(-1);
    }
    createIndex();
}


void InferenceTree::remove( InferenceNode* src ){
    InferenceNode* parent = (InferenceNode*)src->getParent();
    Tree::remove( src );
    if(pmodel){
        parent->updateIdenticalRecursively();
        parent->invalidateRecursively(-1);
    }
}


void InferenceTree::add( InferenceNode* src, InferenceNode* dest, double distance ){
    Tree::add( src, dest, distance );
    if(pmodel){
        //check if the new node has been initialised
        assert( src->partialLikelihood.size() );
        dest->updateIdenticalRecursively();
        dest->invalidateRecursively(-1);
    }
}




double InferenceTree::loglikelihood() {
    assert ( ptable && pmodel );

    //BEWARE, the order of processing is important
    //compute the children likelihood before
    for ( vector< BasicNode*>::reverse_iterator iter = nodeRefVector.rbegin();
          iter != nodeRefVector.rend(); ++iter ){
        if (!(*iter)->isLeaf()){
            ((InferenceNode*)(*iter))->likelihood();
        }
    }


    //ofstream ost("temp.l");


    double likelihood = 0.0;
    for( unsigned int cat = 0; cat < ptable->getNumberCategories(); ++cat ){
        int sequenceLength = ptable->getSequencesLength(cat);
        int numberRateCategories = pmodel->getNumberRatesCategories(cat);
        double rateProb[numberRateCategories];
        for ( int rateCat = 0; rateCat < numberRateCategories; ++rateCat ) {
            rateProb[rateCat] =
            pmodel->getRateCategoryProbability( rateCat, cat );
        }
        int numberStates = pmodel->getNumberStates(cat);
        double frequency[numberStates][numberRateCategories];
        for ( int state = 0; state < numberStates; ++state ) {
            for ( int rateCat = 0; rateCat < numberRateCategories; ++rateCat ) {
                frequency[state] [rateCat] = pmodel->getFrequency( state,
                                                         rateCat, cat );
            }
        }

        for ( int site = 0; site < sequenceLength; ++site ) {
            double tempr = 0.0;
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                double tempn = 0.0;
                for ( int nucleotide = 0; nucleotide < numberStates;
                      ++nucleotide ) {
                    tempn += frequency[nucleotide][rate]*
                        (*(((InferenceNode*)root)->
                           partialLikelihood[cat]))(site, nucleotide, rate);
                }
                tempr += ( rateProb[rate] * tempn );
            } // end rate
            likelihood += log( tempr );
           // ost << log( tempr ) << ' ';
        }//end site


        //invariant
        const vector< pair< string, unsigned int > > &inv =
            ptable->getInvariantBases(cat);
        for ( unsigned int invId = 0; invId < inv.size(); ++invId ) {
            double tempr = 0.0;
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                double tempn = 0.0;
                for ( int nucleotide = 0; nucleotide < numberStates;
                      ++nucleotide ) {
                    tempn += frequency[nucleotide] [rate]  *
                    (*(((InferenceNode*)root)->
                           partialLikelihood[cat]))
                                  (sequenceLength+invId, nucleotide, rate);
                }
                tempr += ( rateProb[rate] * tempn );
            }
            likelihood += ( log( tempr ) * inv[invId].second );

         //   for (unsigned int tot = 0; tot < inv[invId].second; ++tot)
         //       ost << log( tempr ) << ' ';
        } // end invId

       // ost << endl << "---------------" << endl;
    }
    return ( likelihood );
}

void InferenceTree::invalidateNodes(int cat) {
    for ( vector<BasicNode*>::iterator iter = nodeRefVector.begin();
          iter != nodeRefVector.end(); ++iter ){
        if( !(*iter)->isLeaf() ){
            //invalidate all
            if (cat == -1){
                ((InferenceNode*)(*iter))->partialLikelihood =
                         ((InferenceNode*)(*iter))->partialLikelihoodWork;
            }
            //invalidate just one category
            else{
                ((InferenceNode*)(*iter))->partialLikelihood[cat] =
                         ((InferenceNode*)(*iter))->partialLikelihoodWork[cat];
            }
        }
    }
}

void InferenceTree::saveNodes(int cat) {
    for ( vector<BasicNode*>::iterator iter = nodeRefVector.begin();
          iter != nodeRefVector.end(); ++iter ){
        if( !(*iter)->isLeaf() ){
            //save all
            if (cat == -1){
                for ( unsigned int cat = 0;
                      cat < getTable()->getNumberCategories(); ++cat ){
                    //if work is used, swap the array with save
                    if (((InferenceNode*)(*iter))->partialLikelihoodSave[cat] != ((InferenceNode*)(*iter))->partialLikelihood[cat]){
                        ((InferenceNode*)(*iter))->partialLikelihoodWork[cat] = ((InferenceNode*)(*iter))->partialLikelihoodSave[cat];
                        ((InferenceNode*)(*iter))->partialLikelihoodSave[cat] = ((InferenceNode*)(*iter))->partialLikelihood[cat];
                    }
                }
            }
            //retrieve just one category
            else{
                if (((InferenceNode*)(*iter))->partialLikelihoodSave[cat] != ((InferenceNode*)(*iter))->partialLikelihood[cat]){
                    ((InferenceNode*)(*iter))->partialLikelihoodWork[cat] = ((InferenceNode*)(*iter))->partialLikelihoodSave[cat];
                    ((InferenceNode*)(*iter))->partialLikelihoodSave[cat] = ((InferenceNode*)(*iter))->partialLikelihood[cat];
                }
            }
        }
    }
}

void InferenceTree::retrieveNodes(int cat) {
    for ( vector<BasicNode*>::iterator iter = nodeRefVector.begin();
          iter != nodeRefVector.end(); ++iter ){
        if( !(*iter)->isLeaf() ){
            //retrieve all
            if (cat == -1){
                ((InferenceNode*)(*iter))->partialLikelihood = ((InferenceNode*)(*iter))->partialLikelihoodSave;
            }
            //retrieve just one category
            else{
                ((InferenceNode*)(*iter))->partialLikelihood[cat] = ((InferenceNode*)(*iter))->partialLikelihoodSave[cat];
            }
        }
    }
}


void InferenceTree::setBranchLength( unsigned int branchId, double newDist ) {
    Tree::setBranchLength( branchId, newDist);
    if(pmodel){
        ((InferenceNode*)(nodeRefVector[branchId+1]->getParent()))->invalidateRecursively(-1);
    }
}

void InferenceTree::constructRandomly( double maxBranchLength ) {

    assert(ptable);

    unsigned int numberSpecies = ptable->getNumberSpecies();
    vector < double > branchLength;
    vector < unsigned int > array( numberSpecies );
    string * leaf = new string[numberSpecies];

    Singleton< randombox > & randBox = Singleton< randombox >::instance();
    // Create a tree by addition of species to randomly selected branches
    for ( unsigned int i = 0; i < numberSpecies; ++i ) {
        array[i] = i;
    }
    array = statlib::random_permutation( array );

    // Read in all species names
    for ( unsigned int i = 0; i < numberSpecies; ++i ) {
        leaf[i] = ptable->species[i];
    }

    root = new InferenceNode(this);
    if (clustersTree){
        const vector< ClustersTreeNode* > clusters = clustersTree->getList();
        //create clusters ancestral node
        vector< ClustersTreeNode* >::const_reverse_iterator iter = clustersTree->getList().rbegin();
        //the biggest cluster contains all the species and should be associated the root node
        assert( (*iter)->cluster->getNumberSpecies() == ptable->getNumberSpecies() );
        (*iter)->cluster->setAncestor(root);
        ++iter;
        //but we can create the other nodes
        while (iter != clustersTree->getList().rend()){
            InferenceNode* node = new InferenceNode(this);
            (*iter)->cluster->setAncestor(node);
            ++iter;
        }
        //and link them (in a "random" binary tree)
        //reset iter... obviously, we do not link the root (hence the iter=iter+1)
        iter = clustersTree->getList().rbegin();
        ++iter;
        while (iter != clustersTree->getList().rend()){
            //if the parent cluster does not have two children yet
            if ( (*iter)->parent->cluster->getAncestor()->getNumberChildren() < 2 ){
                add( (InferenceNode*)(*iter)->cluster->getAncestor(),
                     (InferenceNode*)(*iter)->parent->cluster->getAncestor(),
                     maxBranchLength * randBox.ran() );
            }
            //otherwise create a new node to be new parent of two clusters
            else{
                InferenceNode* oldParent = (InferenceNode*)(*iter)->parent->cluster->getAncestor();
                assert( oldParent->getNumberChildren() == 2 );
                InferenceNode* movedNode = oldParent->getChild((int)(2*randBox.ran()));
                double oldDistance = movedNode->getParentDistance();
                InferenceNode* newNode = new InferenceNode(this);
                remove( movedNode );
                add( newNode, oldParent, maxBranchLength * randBox.ran() );
                add( (InferenceNode*)(*iter)->cluster->getAncestor(), newNode, maxBranchLength * randBox.ran() );
                add( movedNode, newNode, oldDistance );
            }
            ++iter;
        }
    }
    //now add the leaves to the tree, respect the clusters if they are defined
    for ( unsigned int i = 0; i < numberSpecies; ++i ) {
        //for each species we look for possible insertion points
        vector< pair<InferenceNode*,double> > branchVector;
        double sum = 0.0;
        ClustersTreeNode* clustersTreeNode = NULL;
        //if clusters are defined
        if(clustersTree){
            //look for the smallest clusters which contains the species i
            clustersTreeNode = clustersTree->findInsertion( array[i] );
            //retrieve all the possible branches (ie branches in clustersTreeNode but not in a smallest
            //included cluster)
            assert((InferenceNode*)clustersTreeNode);
            assert((InferenceNode*)clustersTreeNode->cluster->getAncestor());
            sum = ((InferenceNode*)clustersTreeNode->cluster->getAncestor())->
                          findClusterBranches(branchVector, clustersTreeNode);
        }
        else{
            sum = ((InferenceNode*)root)->findClusterBranches(branchVector, NULL);
        }
        //if one branch or less was found then just plug the leaf directly where it should be
        if (branchVector.size() <= 1){
            if(clustersTreeNode){
                add( new InferenceNode(leaf[array[i]], this ),
                     (InferenceNode*)clustersTreeNode->cluster->getAncestor(),
                     maxBranchLength * randBox.ran() );
            }
            else{
                add( new InferenceNode(leaf[array[i]], this ),
                     (InferenceNode*)root,
                     maxBranchLength * randBox.ran() );
            }
        }
        else{
            //branchVector contains possible insertion branches with length information
            //choose the insertion point with probability proportionnal to the "free space"
            double ins = sum * randBox.ran();
            vector< pair<InferenceNode*,double> >::iterator iter = branchVector.begin();
            while(ins>0.0){
                assert(iter!=branchVector.end());
                ins -= (*iter).second;
                ++iter;
            }
            //the selected branch is defined by the node in --iter
            --iter;
            InferenceNode* movedNode = (InferenceNode*)((*iter).first);
            InferenceNode* oldParent = (InferenceNode*)movedNode->getParent();
            double oldDistance = (*iter).second;
            assert(oldDistance==movedNode->getParentDistance());
            InferenceNode* newNode = new InferenceNode(this);
            remove( movedNode );
            add( newNode, oldParent, ins + oldDistance );
            add( new InferenceNode(leaf[array[i]],this), newNode, -ins );
            add( movedNode, newNode, maxBranchLength * randBox.ran() );
        }
    }
    //createIndex(); useless because it is done in Tree::add()
    delete [] leaf;
}

void InferenceTree::constructFromString( string stringTree ){
    assert(ptable);
    // Create a depth map of the tree
    TreeMap treeMap( stringTree );
    // Create the nodes
    root = new InferenceNode( treeMap, this );
    createIndex();
    if (clustersTree){
        assignClusters();
    }
}

#ifdef DEBUG1
void InferenceTree::checkIdentical( ){
    cout << "WARNING: checking identical" << endl;
    for ( vector<BasicNode*>::iterator iter = nodeRefVector.begin();
          iter != nodeRefVector.end(); ++iter ){
        ((InferenceNode*)(*iter))->checkIdentical();
    }
}
#endif

void InferenceTree::loadClusters( const string& fileName ){
    clustersStream = new ifstream( fileName.c_str() );
    if (!clustersStream){
        cerr << "The clusters file: " << fileName << " seems invalid" << endl;
        exit(EXIT_FAILURE);
    }
    if (ptable){
        processClusters();
    }
}

void InferenceTree::processClusters(){
    assert(clustersStream != NULL);
    assert(ptable);
    clustersTree = new ClustersTree(*clustersStream, ptable);
}

void InferenceTree::assignClusters(){
    Cluster * cluster = ((InferenceNode*)root)->assignClustersRec( clustersTree );
    assert(cluster->getNumberSpecies()==ptable->getNumberSpecies());
    delete cluster;
    //check whether the assignment is complete
    for ( ClustersTree::const_iterator iter = clustersTree->begin();
          iter != clustersTree->end(); ++iter ){
        if ( (*iter).second.first->getAncestor() == NULL){
            cerr << "Sorry, the cluster " << (*iter).second.first->getName()
                 << " could not be assigned to a node in the tree:" << endl;
            cerr << toString(true) << ';' << endl;
            exit(EXIT_FAILURE);
        }
    }
}

#ifdef DEBUG1
void InferenceTree::checkClusters(){
    cout << "WARNING: checking cluster" << endl;
    if (clustersTree){
        for ( vector< ClustersTreeNode* >::const_iterator iter = clustersTree->getList().begin();
              iter != clustersTree->getList().end(); ++iter ){
            Cluster* cluster = ((InferenceNode*)((*iter)->cluster->getAncestor()))->getCluster();
            if( !(*iter)->cluster->matches(*cluster) ){
                cerr << "Cluster mismatch: " << (*iter)->cluster->getName() << endl;
                cerr << ((InferenceNode*)((*iter)->cluster->getAncestor()))->toString(true) << endl;
                vector<unsigned int> mask = (*iter)->cluster->getMask();
                cerr << "expected mask=";
                for ( unsigned int i = 0; i < mask.size(); ++i ){
                    cerr << mask[i] << '-';
                }
                cerr << endl;
                mask = cluster->getMask();
                cerr << "observed mask=";
                for ( unsigned int i = 0; i < mask.size(); ++i ){
                    cerr << mask[i] << '-';
                }
                cerr << endl;
                exit(EXIT_FAILURE);
            }
            delete cluster;
        }
    }
}
#endif

InferenceNode* InferenceTree::findClusterAncestor(InferenceNode* node) const{
    ClustersTreeNode* clustersTreeNode = findInsertion( node );
    return ((InferenceNode*)clustersTreeNode->cluster->getAncestor());
}

double InferenceTree::findClusterBranches( vector< pair<InferenceNode*,double> > & branchVector,
                                           ClustersTreeNode* clustersTreeNode ) const{
    assert(clustersTree);
    assert((InferenceNode*)clustersTreeNode->cluster->getAncestor());
    return ((InferenceNode*)clustersTreeNode->cluster->getAncestor())->findClusterBranches(branchVector, clustersTreeNode);
}

ClustersTreeNode* InferenceTree::findInsertion( InferenceNode* node ) const{
    assert(clustersTree);
    return findInsertionRec( node, clustersTree->getRoot() );
}

ClustersTreeNode* InferenceTree::findInsertionRec( InferenceNode* node, ClustersTreeNode* ancestorCluster ) const{
    assert( ancestorCluster->cluster->getAncestor() );
    for ( list<ClustersTreeNode*>::const_iterator iter = ancestorCluster->begin();
          iter != ancestorCluster->end(); ++ iter ){
        if ( node == (*iter)->cluster->getAncestor() ){
            return (*iter);
        }
        if ((*iter)->cluster->getAncestor()){
            if ( node->isDescendant( (*iter)->cluster->getAncestor() ) ){
                return findInsertionRec( node, *iter );
            }
        }
    }
    //node is a descendant of ancestorCluster but is not the descendant of any cluster included
    //in ancestorCluster... therefore, ancestorCluster is the insertion point
    return ancestorCluster;
}

double InferenceTree::getClusterDistance( InferenceNode* node, ClustersTreeNode*& orig ) const{
    //find stoppingPoint the ancestor of node we are interested in
    //and return the distance stoppingPoint -----> node
    BasicNode* stoppingPoint;
    //if clusters were defined we want the distance to node from
    //the most recent cluster ancestor which contains it.
    if ( clustersTree ){
        orig = findInsertion( node );
        stoppingPoint = orig->cluster->getAncestor();
        //if node is already a cluster ancestor we are interested in
        //the distance to the parent cluster
        if (stoppingPoint==node){
            if ( !orig->getParent() ){
                assert(node==root);
                //stoppingPoint = root;
            }
            else{
                orig = orig->getParent();
                stoppingPoint = orig->cluster->getAncestor();
            }
        }
    }
    //if no cluster has been defined we are interested in the distance to the root
    else{
        orig = NULL;
        stoppingPoint = root;
    }

    double dist = 0.0;
    while (node != stoppingPoint){
        dist += node->getParentDistance();
        node =(InferenceNode*)(node->getParent());
    }
    return dist;
}

double InferenceTree::getRootDistance( BasicNode* node ) const{
    double dist = 0.0;
    while (node != root){
        assert(node);
        dist += node->getParentDistance();
        node = node->getParent();
    }
    return dist;
}
