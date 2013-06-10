#include "Tree/ConsensusTree.h"

#include "Util/ParametersSet.h"
#include "Util/FileParser.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Sequence/SequenceTable.h"
#include "Tree/BasicNode.h"
#include "Tree/Cluster.h"
#include "Tree/LengthCluster.h"
#include "Tree/ClustersSet.h"

#include <algorithm>


//Convention used in consensus : the name must be the name of the
//corresponding MCMC tree with consensus instead of MCMC
ConsensusTree ConsensusTree::prototype("Unrooted consensus tree");


ConsensusTree::ConsensusTree( const string & registrationName ):
consensusFileParameters("consensus") {
    Singleton< Factory<ConsensusTree> > & consensusFactory =
            Singleton< Factory<ConsensusTree> >::instance();
    consensusFactory.subscribe( this, registrationName );
}


ConsensusTree::ConsensusTree( ParametersSet& parameters ):
consensusFileParameters("consensus"){
    string outputFileBaseName( parameters.stringParameter( "Output file" ) );
    string samplesFileName = outputFileBaseName + ".smp";
    string branchLengthsFileName = outputFileBaseName + ".bl";
    samplesFile.open(samplesFileName.c_str(), ifstream::in );
    branchLengthsFile.open(branchLengthsFileName.c_str(), ifstream::in );
    if (!samplesFile.is_open()){
        cerr << "Error while opening the samples file "
             << samplesFileName << "...abort..." << endl;
        exit(EXIT_FAILURE);
    }
    if (!branchLengthsFile.is_open()){
        cerr << "Error while opening the branch lengths file "
             << branchLengthsFileName << "...abort..." << endl;
        exit(EXIT_FAILURE);
    }
    expectedNumberSamples = parameters.intParameter("Sampling iterations");
    expectedNumberSamples = expectedNumberSamples /
                    parameters.intParameter("Sampling period");
    FileParser::confirmOpenFile( consensusFile, outputFileBaseName + ".cns", true );
    FileParser::confirmOpenFile( consensusTreeFile, outputFileBaseName + ".cnt", true );
    //the filename to use in the .cns file for the .cnt file (remove the path)
    //look for the last '/' or '\' (can be both for windows systems)
    string::reverse_iterator iter = find(outputFileBaseName.rbegin(), outputFileBaseName.rend(), '/');
    iter = find(outputFileBaseName.rbegin(), iter, '\\');
    string inFileBaseName( iter.base(), outputFileBaseName.end() );
    consensusFileParameters["Trees file"] = inFileBaseName + ".cnt";
    seq = new SequenceTable(parameters( "DATAFILE" ));
}


ConsensusTree::~ConsensusTree(){
    if(seq){
        delete seq;
    }
    if(clades){
        delete clades;
    }
}

ConsensusTree* ConsensusTree::clone( ParametersSet& parameters ) const{
    return new ConsensusTree( parameters );
}

void ConsensusTree::readTree( istream& inputFile, string & stringTree ){
    FileParser::readTree( inputFile, stringTree );
    inputFile >> ws;
    if (stringTree.length() == 0 ){
        cerr << "error while reading the \".samples\" file, sample no:"
              << foundNumberSamples + 1 << endl;
        exit(EXIT_FAILURE);
    }
    inputFile >> ws;
}

void ConsensusTree::readBranchLengths( istream& inputFile,
                                        vector<double> & branchLengths,
                                        unsigned int expectedNumber ){
    static char line[BLFILE_LINELENGTH];
    char* linep;
    int numberChar;
    double length;

    inputFile >> ws;
    inputFile.getline(line, BLFILE_LINELENGTH);
    inputFile >> ws;
    linep = line;
    branchLengths.clear();
    while ( sscanf( linep, "%lf%n", &length, &numberChar ) == 1 ){
        linep += numberChar;
        branchLengths.push_back(length);
    }
    if (branchLengths.size() != expectedNumber){
        cerr << "Error while processing the sample " << foundNumberSamples
             << '.' << endl
             << "The number of branch lengths (" << branchLengths.size()
             << ") does not match the corresponding tree ("
             << expectedNumber << " branches expected)"
             << endl;
        exit(EXIT_FAILURE);
    }
}


unsigned int ConsensusTree::process(){
    string sampledTreeString;
    foundNumberSamples = 0;

    //initialise the clusters set used to contain all the future clades found
    clades = new ClustersSet(seq->getNumberSpecies());

    vector< Cluster* > treeClades;

    vector< double > branchLengths;
    cout << "Expected number of samples: " << expectedNumberSamples << endl;
    while ( !samplesFile.eof() ) {
        ++foundNumberSamples;
        if (foundNumberSamples%1000==0){
            cout << foundNumberSamples << endl;
        }
        //read the next tree
        readTree(samplesFile, sampledTreeString);
        Tree sampledTree(sampledTreeString);
        //read the next set of branch lengths
        readBranchLengths(branchLengthsFile, branchLengths,
                                    sampledTree.getNumberBranches() );
        //transform the specie's IDs into the specie's names (because the sorting mechanism will
        //might mess up the branches otherwise) and store the ID in the appropriate field
        for( vector < BasicNode * >::iterator iter = sampledTree.nodeRefVector.begin();
             iter != sampledTree.nodeRefVector.end(); ++iter ){
            if ( (*iter)->isLeaf() ){
                (*iter)->cnumber = atoi( ((*iter)->getLabel()).c_str() );
                (*iter)->changeLabel( seq->species[(*iter)->cnumber-1] );
            }
        }
        //resort the tree accordingly and include the branches
        sampledTree.createIndex();
        sampledTree.setBranchVector(branchLengths);
        //prepare the tree and process its clade
        processClades( sampledTree, treeClades );
        for ( vector< Cluster* >::iterator iter = treeClades.begin();
              iter != treeClades.end(); ++iter ){
            clades->add((*iter));
            delete (*iter);
        }
        treeClades.clear();
    }
    if (foundNumberSamples != expectedNumberSamples){
        cerr << "Warning: the sample size and the expected size are different"
             << endl;
    }
    if ( !branchLengthsFile.eof() ){
        cerr << "Error: There are more samples in the branch lengths file than "
             << "in the samples file. Abort..." << endl;
        exit(EXIT_FAILURE);
    }
    return foundNumberSamples;
}


void ConsensusTree::produceConsensus(){

    multimap< unsigned int, Cluster* > cladesList = clades->listing();
    vector< pair<unsigned int, Cluster*> > included;
    multimap< unsigned int, Cluster* > treeClades;
    vector< pair<unsigned int, Cluster*> > excluded;

    multimap< unsigned int, Cluster* >::reverse_iterator iterMap;
    vector< pair<unsigned int, Cluster*> >::iterator iterVector;

    //find the clades to insert in the final consensus tree
    //and order them (smaller ones first)
    iterMap = cladesList.rbegin();
    while( ( (*iterMap).first > (foundNumberSamples + 1) / 2 ) &&
           ( iterMap != cladesList.rend() ) ){
        treeClades.insert( pair<unsigned int, Cluster*>
              ( ((*iterMap).second)->getNumberSpecies(), (*iterMap).second ) );
        ((*iterMap).second)->setSupport(100.0 *(double)((*iterMap).first)/(double)foundNumberSamples);
        unsigned int cladeSize = ((*iterMap).second)->getNumberSpecies();
        if ( ( cladeSize != 1 ) && ( cladeSize != seq->getNumberSpecies() )
            //&& ( cladeSize != seq->getNumberSpecies() - 1 ) //ok for unrooted
                ){
            included.push_back( (*iterMap) );
        }
        ++iterMap;
    }
    // below 50%, add clusters only if they are not incompatible with
    // previous ones
    while( iterMap != cladesList.rend() ){
        bool compatible = true;
        iterVector = included.begin();
        while ( compatible && iterVector != included.end() ){
            //if intersection...
            if (((*iterMap).second)->intersect(*((*iterVector).second))){
                //...without inclusion => incompatible
                if ( !(((*iterVector).second)->contains(*((*iterMap).second))) &&
                     !(((*iterMap).second)->contains(*((*iterVector).second))) ){
                    compatible = false;
                }
             }
            ++iterVector;
        }
        if (compatible){
            treeClades.insert( pair<unsigned int, Cluster*>
              ( ((*iterMap).second)->getNumberSpecies(), (*iterMap).second ) );
            ((*iterMap).second)->setSupport(100.0 *(double)((*iterMap).first)/(double)foundNumberSamples);
            included.push_back( (*iterMap) );
        }
        else{
            excluded.push_back( (*iterMap) );
        }
        ++iterMap;
    }
    vector<unsigned int> spec;
    cout << "Sets included in the consensus" << endl;
    for ( iterVector = included.begin(); iterVector != included.end();
          ++iterVector ){
        ((*iterVector).second)->getSpecies( spec );
        vector<unsigned int>::iterator iterSpec = spec.begin();
        unsigned int i = 0;
        while(iterSpec!=spec.end()){
            while(i!=*iterSpec){
                cout << '.';
                ++i;
            }
            cout << '*';
            ++i;
            ++iterSpec;
        }
        while(i < seq->getNumberSpecies() ){
            cout << '.';
            ++i;
        }
        cout << "   "
             << 100.0 * (double)(*iterVector).first/(double)foundNumberSamples
             << '%' << endl;
    }
    cout << endl << "Sets NOT included in the consensus" << endl;
    for ( iterVector = excluded.begin(); iterVector != excluded.end();
          ++iterVector ){
        ((*iterVector).second)->getSpecies( spec );
        vector<unsigned int>::iterator iterSpec = spec.begin();
        unsigned int i = 0;
        while(iterSpec!=spec.end()){
            while(i!=*iterSpec){
                cout << '.';
                ++i;
            }
            cout << '*';
            ++i;
            ++iterSpec;
        }
        while(i < seq->getNumberSpecies() ){
            cout << '.';
            ++i;
        }
        cout << "   "
             << 100.0 * (double)(*iterVector).first/(double)foundNumberSamples
             << '%' << endl;
    }
    produceTree( treeClades );
}

void ConsensusTree::produceTree(multimap< unsigned int, Cluster* > & treeClades){
    vector<unsigned int> species;
    multimap< unsigned int, Cluster* >::iterator iterMap1;
    multimap< unsigned int, Cluster* >::iterator iterMap2;
    multimap< unsigned int, Cluster* >::iterator iterErase;
    BasicNode* father;
    BasicNode* child;

    multimap< unsigned int, Cluster* > treeCladesCopy = treeClades;
    for ( iterMap1 = treeClades.begin(); iterMap1 != treeClades.end();
            ++iterMap1 ){
        //create leaves
        if ((*iterMap1).first==1){
            ((*iterMap1).second)->getSpecies( species );
            unsigned int speciesId = species[0];
            ((*iterMap1).second)->setAncestor( new BasicNode(
                                             seq->species[speciesId]) );
        }
        else{
            father = new BasicNode();
            ((*iterMap1).second)->setAncestor( father );
            iterMap2 = treeClades.begin();
            while (iterMap2 != iterMap1){
                if ( ((*iterMap1).second)->contains( *((*iterMap2).second)) ){
                    child = ((*iterMap2).second)->getAncestor();
                    father->addChild( child, getDist(((*iterMap1).second), ((*iterMap2).second)) );
                    iterErase = iterMap2;
                    ++iterMap2;
                    treeClades.erase( iterErase );
                }
                else{
                    ++iterMap2;
                }
            }
        }
    }
    setRoot( ((*(--treeClades.end())).second)->getAncestor() );


    for ( iterMap1 = treeCladesCopy.begin(); iterMap1 != treeCladesCopy.end();
            ++iterMap1 ){
        //create leaves
        if ((*iterMap1).first==1){
            ((*iterMap1).second)->getSpecies( species );
            unsigned int speciesId = species[0];
            ((*iterMap1).second)->setAncestor( new BasicNode(
                                             seq->species[speciesId]) );
        }
        else{
            father = new BasicNode();
            ((*iterMap1).second)->setAncestor( father );
            iterMap2 = treeCladesCopy.begin();
            while (iterMap2 != iterMap1){
                if ( ((*iterMap1).second)->contains( *((*iterMap2).second)) ){
                    child = ((*iterMap2).second)->getAncestor();
                    father->addChild( child, ((*iterMap2).second)->getSupport() );
                    iterErase = iterMap2;
                    ++iterMap2;
                    treeCladesCopy.erase( iterErase );
                }
                else{
                    ++iterMap2;
                }
            }
        }
    }
    supportTree.setRoot( ((*(--treeCladesCopy.end())).second)->getAncestor() );
}


Cluster* ConsensusTree::processClades( const Tree& sampledTree,
                                       vector< Cluster* > & treeClades,
                                       BasicNode* startPoint ){

    if (startPoint == NULL){
        startPoint = sampledTree.getRoot();
    }

    Cluster* cluster = createCluster( startPoint );
    treeClades.push_back(cluster);
    if (startPoint->isLeaf()){
         unsigned int id = startPoint->cnumber - 1;
         cluster->addSpecies( id );
    }
    else{
        list< BasicNode* >::const_iterator iter;
        for ( iter = startPoint->getChildrenList().begin();
              iter != startPoint->getChildrenList().end(); ++iter ){
            Cluster* retCluster = processClades( sampledTree, treeClades, *iter );
            cluster->join( *retCluster );
        }
    }
    return cluster;
}

Cluster* ConsensusTree::createCluster(BasicNode* node){
    return new LengthCluster( seq->getNumberSpecies(), node->getParentDistance() );
}


void ConsensusTree::prepareWriting(){
     consensusFileParameters["Number of trees"] = "2";
     consensusTreeFile << toString() << ';' << endl;
     consensusTreeFile << supportTree.toString() << ';' << endl;
     consensusFileParameters("TREE1")["Parameter name"] = "Branch lengths";
     consensusFileParameters("TREE1")["Branch lengths"] = "yes";
     consensusFileParameters("TREE1")["Branch labels"] = "yes";
     consensusFileParameters("TREE2")["Parameter name"] = "Clade supports";
     consensusFileParameters("TREE2")["Clades indices"] = "yes";
}

void ConsensusTree::writeConsensusFile(){
     prepareWriting();
     consensusFileParameters.saveToFile(consensusFile);
     consensusFile.close();
     consensusTreeFile.close();
}

double ConsensusTree::getDist( Cluster* father, Cluster* child ){
    return ((LengthCluster*)child)->getLength()/
           (double)((LengthCluster*)child)->getNumber();
    father = father;
}
