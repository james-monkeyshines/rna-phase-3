#include "Tree/ClustersTree.h"

#include "Tree/Cluster.h"
#include "Tree/ClustersTreeNode.h"
#include "Tree/Tree.h"

#include "Sequence/SequenceTable.h"
#include "Util/FileParser.h"

#include <assert.h>

using namespace std;

ClustersTree::ClustersTree(SequenceTable* table):ClustersSet( table->getNumberSpecies() ){
     //save a pointer to the sequence table
    this->table = table;

    //create the root cluster which contains all species
    Cluster* clusterRoot = new Cluster(totalNumberSpecies, "root");
    clusterRoot->invert();
    ClustersSet::insert(clusterRoot);
    root = new ClustersTreeNode(clusterRoot);
    makeList();
}

ClustersTree::ClustersTree( ifstream& fileStream, SequenceTable* table ):
ClustersSet( table->getNumberSpecies() ){

    //save a pointer to the sequence table
    this->table = table;

    //create the root cluster which contains all species
    Cluster* clusterRoot = new Cluster(totalNumberSpecies, "root");
    clusterRoot->invert();
    ClustersSet::insert(clusterRoot);
    root = new ClustersTreeNode(clusterRoot);

    // Read in clusters
    string clusterName;
    char speciesName[100];
    fileStream >> ws;
    if (fileStream.eof()){
        cerr << "Sorry, invalid clusters files" << endl;
        exit(EXIT_FAILURE);
    }
    cout << endl;
    while ( !fileStream.eof() ) {
        fileStream >> clusterName;
        fileStream >> ws;
        if (fileStream.peek()!='('){
            cout << "Cluster " << clusterName << ": ";
            if (table->findSequence(clusterName)!=-1){
                cerr << "Sorry, clusters cannot be given the name of existing species" << endl;
                exit(EXIT_FAILURE);
            }
            char c;
            Cluster* newCluster = new Cluster(totalNumberSpecies, clusterName);
            do{
                char* pName = speciesName;
                fileStream.get( c );
                while ( ( c != ';' ) && ( !isspace(c) ) && !fileStream.eof() ){
                    *pName = c;
                    ++pName;
                    fileStream.get( c );
                }
                *pName = 0;
                int speciesId = table->findSequence(speciesName);
                if (speciesId==-1){
                    cerr << endl << "Unknown species: " << speciesName
                         << " found in " << clusterName << " (missing ';' ?)...aborting..." << endl;
                    exit(EXIT_FAILURE);
                }
                newCluster->addSpecies( (unsigned int)speciesId );
                fileStream >> ws;
                //c=fileStream.peek();
            }while ( ( c != ';' ) && !fileStream.eof() );
            if ( c != ';' ){
                cerr << endl << "Error while reading cluster file, (missing ';' after cluster "
                     << clusterName << "?)" << endl;
                exit(EXIT_FAILURE);
            }
            cout << newCluster->getNumberSpecies() << " species";
            if ( newCluster->getNumberSpecies() == 1 ){
                cerr << "Cluster must have more than one species" << endl;
                exit(EXIT_FAILURE);
            }
            insert( newCluster );
            fileStream >> ws;
            cout << "; created..." << endl;
        }
        //if first character is a '(', create the cluster for the topology
        else{
            cout << "Topology " << clusterName << ": ";
            string stringTree;
            FileParser::readTree( fileStream, stringTree );
            unsigned int numberClusters = 0;
            addTreeClusters(stringTree, clusterName, numberClusters, table);
            fileStream >> ws;
            cout << numberClusters << " clusters created..." << endl;
        }
    }
    makeList();
}

Cluster* ClustersTree::addTreeClusters( const string& stringTree, const string& clusterName,
                                              unsigned int& clusterId, SequenceTable* table ){

    char name[100];

    Tree tree(stringTree);
    list< BasicNode* > children = tree.getRoot()->getChildrenList();
    if (clusterId!=0){
        sprintf( name, "%s_%d", clusterName.c_str(), clusterId );
    }
    else{
        sprintf( name, "%s", clusterName.c_str() );
    }
    ++clusterId;
    Cluster* newCluster = new Cluster(totalNumberSpecies, name);
    for ( list< BasicNode* >::iterator iter = children.begin();
          iter != children.end(); ++iter ){
        if (!(*iter)->isLeaf()){
            Cluster* childCluster = addTreeClusters( (*iter)->toString(true), clusterName, clusterId, table );
            newCluster->join( *childCluster );
        }
        else{
            int speciesId = table->findSequence((*iter)->getLabel());
            if (speciesId==-1){
                cerr << endl << "Unknown species: " << (*iter)->getLabel()
                     << "...aborting..." << endl;
                exit(EXIT_FAILURE);
            }
            newCluster->addSpecies( (unsigned int)speciesId );
        }
    }
    insert( newCluster );
    return ( newCluster );
}

bool ClustersTree::unrootCheck(){
    //there is a risk of complementary cluster if there are two children at the root
    if (root->size() == 2){
        unsigned int numberSpecies = root->front()->cluster->getNumberSpecies() +
                                     root->back()->cluster->getNumberSpecies();
        if ( table->getNumberSpecies() == numberSpecies ){
            return false;
        }
    }
    return true;
}

ClustersTree::~ClustersTree(){
    delete root;
}

bool ClustersTree::insert(Cluster* cluster){
    if (!ClustersSet::insert(cluster)){
        cerr << "Error, duplicated cluster definition, cannot add " << cluster->getName() << endl;
        exit(EXIT_FAILURE);
    }
    //find the insertion point for the new cluster
    ClustersTreeNode* insertionPoint = root->findInsertion( cluster );
    assert(insertionPoint);
    if (!insertionPoint->cluster->compatible( *cluster ) ){
        cerr << "Cluster " << cluster->getName() << " conflicts with " << (insertionPoint->cluster)->getName() << endl;
        exit(EXIT_FAILURE);
        return false;
    }
    assert(cluster);
    ClustersTreeNode* clustersNode = new ClustersTreeNode(cluster);
    list<ClustersTreeNode*>::iterator iter = insertionPoint->begin();
    while ( iter != insertionPoint->end() ){
        assert((*iter)->cluster);
        if (cluster->contains(*((*iter)->cluster))){
            clustersNode->add( *iter );
            insertionPoint->erase(iter++);
        }
        else{
            ++iter;
        }
    }
    insertionPoint->add(clustersNode);
    return true;
}


void ClustersTree::makeList(){
    listNode.clear();
    root->makeListRec(listNode);
}



ClustersTreeNode* ClustersTree::findInsertion( unsigned int speciesId ) const{
    return root->findInsertion( speciesId );
}

ClustersTreeNode* ClustersTree::findInsertion( Cluster* cluster ) const{
    return root->findInsertion( cluster );
}

ClustersTreeNode* ClustersTree::findByAncestor( InferenceNode* node ) const{
    for ( vector< ClustersTreeNode* >::const_iterator iter = listNode.begin();
          iter != listNode.end(); ++iter ){
        if ( (*iter)->cluster->getAncestor() == node ){
            return *iter;
        }
    }
    return NULL;
}
