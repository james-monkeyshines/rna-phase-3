#include "Tree/HeterogeneousConsensusTree.h"

#include "Util/ParametersSet.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Sequence/SequenceTable.h"
#include "Tree/BasicNode.h"
#include "Tree/Cluster.h"
#include "Tree/LengthCluster.h"
#include "Tree/HeterogeneousCluster.h"
#include "Tree/ClustersSet.h"

#include "Models/Heterogeneous.h"

//Convention used in consensus : the name must be the name of the 
//corresponding MCMC tree with consensus instead of MCMC
HeterogeneousConsensusTree HeterogeneousConsensusTree::prototype("Heterogeneous consensus tree");


HeterogeneousConsensusTree::HeterogeneousConsensusTree( const string& registrationName ):
ConsensusTree(registrationName){}


HeterogeneousConsensusTree::~HeterogeneousConsensusTree(){
}


HeterogeneousConsensusTree::HeterogeneousConsensusTree(ParametersSet& parameters):
ConsensusTree( parameters ) {
    string fileName;
    string outputFileBaseName( parameters.stringParameter( "Output file" ) );
    fileName = outputFileBaseName + ".bm";
    branchModelsFile.open(fileName.c_str(), ifstream::in );    
    if (branchModelsFile.eof()){
        cerr << "Error while opening the branch models file " 
             << fileName << "...abort..." << endl;
        exit(EXIT_FAILURE);
    }
    fileName = outputFileBaseName + ".mp";
    modelParametersFile.open(fileName.c_str(), ifstream::in );    
    if (modelParametersFile.eof()){
        cerr << "Error while opening the branch models file " 
             << fileName << "...abort..." << endl;
        exit(EXIT_FAILURE);
    }
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model> >::instance();
    model = dynamic_cast<Heterogeneous*>(modelFactory.create(
        parameters( "MODEL" ).stringParameter("Model"), parameters( "MODEL" )));
    if (model==NULL){
        cerr << "Unexpected error." << endl
             << "The model in the control file cannot be "
             << "interpreted as a heterogeneous model. " << endl
             << "However you are using an heterogeneous tree." << endl
             << "Did you change the control file after the MCMC run?" << endl;
    }
    model->initialisation( NULL );
}


ConsensusTree* HeterogeneousConsensusTree::clone( ParametersSet& parameters ) const{
    return new HeterogeneousConsensusTree( parameters );
}


void HeterogeneousConsensusTree::readBranchModels( istream& inputFile,
                                        vector<unsigned int> & branchModels,
                                        unsigned int expectedNumber ){
    static char line[BLFILE_LINELENGTH];
    char* linep;
    int numberChar;
    unsigned int model;
    
    inputFile >> ws;
    inputFile.getline(line, BLFILE_LINELENGTH);
    inputFile >> ws;
    linep = line;
    branchModels.clear();
    while ( sscanf( linep, "%d%n", &model, &numberChar ) == 1 ){
        linep += numberChar;
        branchModels.push_back(model);
    }
    if (branchModels.size() != expectedNumber){
        cerr << "Error while processing the sample " << foundNumberSamples
             << '.' << endl
             << "The number of branches models (" << branchModels.size()
             << ") does not match the corresponding tree ("
             << expectedNumber << " branches expected)"
             << endl;
        exit(EXIT_FAILURE);
    }
}

void HeterogeneousConsensusTree::readModelParameters( istream& inputFile,
                                         vector<double> & averageSubstitutionRate,
                                         vector< vector<double> > & modelParameters){
    static char line[MPFILE_LINELENGTH];
    char* linep;
    int numberChar;
    double value;
    
    inputFile >> ws;
    inputFile.getline(line, MPFILE_LINELENGTH);
    inputFile >> ws;
    linep = line;
    
    averageSubstitutionRate.clear();
    modelParameters.clear();
    modelParameters.resize(model->getNumberModels());
    
    unsigned int numberParameters = model->getModel(0)->getNumberLineParameters();
    
    for( unsigned int i = 0; i < model->getNumberSymbolCategory(); ++i ){
        for( unsigned int j = 0; j < model->ancestralFrequencies[i]->getNumberLineParameters(); ++j ){
            if (sscanf( linep, "%lf%n", &value, &numberChar ) != 1){
                cerr << "Error while processing the sample " << foundNumberSamples
                     << '.' << endl
                     << "Cannot read an ancestral frequencies for the category "
                     << i+1 << " in the .mp file" << endl;
                exit(EXIT_FAILURE);
            }
            linep += numberChar;
        }
    }
            
    for( unsigned int i = 0; i < model->getNumberModels(); ++i ){
        if (i==0){
            //by default the average substitution rate of the first model is 1
            averageSubstitutionRate.push_back(1.000);
        }
        else{
            if (sscanf( linep, "%lf%n", &value, &numberChar ) != 1){
                cerr << "Error while processing the sample " << foundNumberSamples
                     << '.' << endl
                     << "Cannot read the average substitution rate of model "
                     << i+1 << " in the .mp file" << endl;
                exit(EXIT_FAILURE);
            }
            linep += numberChar;
            averageSubstitutionRate.push_back(value);
        }
        for( unsigned int j = 0; j < numberParameters; ++j ){
            if (sscanf( linep, "%lf%n", &value, &numberChar ) != 1){
                cerr << "Error while processing the sample " << foundNumberSamples
                     << '.' << endl
                     << "Cannot read a parameter of model "
                     << i+1 << " in the .mp file" << endl;
                exit(EXIT_FAILURE);
            }
            linep += numberChar;
            modelParameters[i].push_back(value);
        }
    }
    if (sscanf( linep, "%lf%n", &value, &numberChar ) == 1){
        cerr << "Error while processing the sample " << foundNumberSamples
             << '.' << endl
             << "Too many parameters on the line" << endl;
        exit(EXIT_FAILURE);
    }
}


unsigned int HeterogeneousConsensusTree::process(){
    string sampledTreeString;
    foundNumberSamples = 0;
    
    //initialise the clusters set used to contain all the future clades found
    clades = new ClustersSet(seq->getNumberSpecies());
    
    vector< Cluster* > treeClades;
            
    vector< double > branchLengths;
    vector< unsigned int > branchModels;
    
    //to store separatly the parameters of each model of the
    //heterogeneous model
    vector < double > averageSubstitutionRate;
    vector < vector<double> > modelParameters;
    
    cout << "Expected number of samples: " << expectedNumberSamples << endl;
    while ( !samplesFile.eof() ) {
        ++foundNumberSamples;
        if ( foundNumberSamples%1000 == 0 ){
            cout << foundNumberSamples << endl;
        }
        //read the next tree
        readTree(samplesFile, sampledTreeString);
        Tree sampledTree(sampledTreeString);
        //transform the specie's IDs into the specie's names (because the sorting mechanism will
        //might mess up the branches otherwise) and store the ID in the appropriate field
        for( vector < BasicNode * >::iterator iter = sampledTree.nodeRefVector.begin();
             iter != sampledTree.nodeRefVector.end(); ++iter ){
            if ( (*iter)->isLeaf() ){
                (*iter)->cnumber = atoi( ((*iter)->getLabel()).c_str() );
                (*iter)->changeLabel( seq->species[(*iter)->cnumber-1] );
            }
        }
        //resort the tree        
        sampledTree.createIndex();
        
        //read the next set of branch lengths
        readBranchLengths(branchLengthsFile, branchLengths,
                                    sampledTree.getNumberBranches() );
        //read the associated models
        readBranchModels(branchModelsFile, branchModels,
                                    sampledTree.getNumberBranches() );
        //read the associated model parameters
        readModelParameters( modelParametersFile, averageSubstitutionRate,
                           modelParameters );
        
        //prepare the tree and find its clades
        sampledTree.setBranchVector(branchLengths);
        processClades( sampledTree, treeClades );
        
        
        vector< Cluster* >::iterator cl = treeClades.begin();
        vector< unsigned int >::iterator bm = branchModels.begin();
                           
        //add model information to each clades
        processModelClades( cl, averageSubstitutionRate,
                            modelParameters, bm, sampledTree );
        
        //clades should be
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
    //transform the specie's IDs into the specie's names
    //    cf ConsensusTree for memory, deleted here
    return foundNumberSamples;
}





void HeterogeneousConsensusTree::processModelClades( vector< Cluster* >::iterator& cl,
    const vector < double > & averageSubstitutionRate,
    const vector < vector<double> > & modelParameters,
    vector< unsigned int >::iterator& bm, Tree& sampledTree,
    BasicNode* startPoint ){
    
    // treeClades and bm are not supposed to be in the same order (makes
    // things harder)
    
    if (startPoint == NULL){
        startPoint = sampledTree.getRoot();
        //no parameters for the root cluster
       ((HeterogeneousCluster*)(*cl))->addModelParameters(
               0.0, vector<double>(), 0 );
        ++cl;
    }
    if (!startPoint->isLeaf()){
        //save an iterator at the current position, oldBM and oldBM + 1 are the
        //children of the node
        vector< unsigned int >::iterator oldBM = bm;
        for (unsigned int i = 0; i < startPoint->getNumberChildren(); ++i ){
            ++bm;
        }
        for ( list< BasicNode* >::const_iterator iter = startPoint->getChildrenList().begin();
              iter != startPoint->getChildrenList().end(); ++iter ){
            ((HeterogeneousCluster*)(*cl))->addModelParameters( averageSubstitutionRate[*oldBM],
                                      modelParameters[*oldBM], *oldBM + 1 );
            ++oldBM;
            ++cl;
            processModelClades( cl, averageSubstitutionRate, modelParameters,
                                bm,  sampledTree, *iter );
        }        
    }
}



void HeterogeneousConsensusTree::produceTree(multimap< unsigned int, Cluster* > & treeClades){
    multimap< unsigned int, Cluster* > treeCladesCopy = treeClades;
    vector<unsigned int> species;
    multimap< unsigned int, Cluster* >::iterator iterMap1;
    multimap< unsigned int, Cluster* >::iterator iterMap2;
    multimap< unsigned int, Cluster* >::iterator iterErase;
    BasicNode* father;
    BasicNode* child;
    
    ConsensusTree::produceTree( treeCladesCopy );
    //on top of the two trees built by the base class, HeteogeneousConsensus
    //build a tree for each parameter of the model 
    parametersTree.resize(model->getModel(0)->getNumberLineParameters()+1+2);
    for( unsigned int i = 0; i < parametersTree.size(); ++i ){
        treeCladesCopy = treeClades;
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
                        father->addChild( child, getValue((*iterMap2).second, i ));
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
        parametersTree[i].setRoot( ((*(--treeClades.end())).second)->getAncestor() );
    }
}


Cluster* HeterogeneousConsensusTree::createCluster(BasicNode* node){
    return new HeterogeneousCluster( seq->getNumberSpecies(),
                              node->getParentDistance() );
}

void HeterogeneousConsensusTree::prepareWriting(){
            
    char numberString[10];
    ConsensusTree::prepareWriting();
    //there are 2 + numberTrees trees.
    //the two first trees are created by the base class ConsensusTree.
    //They are for the branch lenghts and the clades supports
    unsigned int numberTrees = parametersTree.size();
    sprintf( numberString, "%d", numberTrees + 2);
    consensusFileParameters["Number of trees"] = numberString;
    consensusTreeFile << parametersTree[0].toString() << ';' << endl;
    consensusTreeFile << parametersTree[1].toString() << ';' << endl;
    consensusFileParameters("TREE3")["Parameter name"] = "MAP model";
    consensusFileParameters("TREE3")["Branch labels"] =  "yes";
    consensusFileParameters("TREE4")["Parameter name"] = "MAP model BPP";
    consensusFileParameters("TREE4")["Horizontal branch colors"] = "yes";
    consensusFileParameters("TREE4")["Branch labels"] =  "yes";
    for (unsigned int j = 2; j < numberTrees; ++j){
        consensusTreeFile << parametersTree[j].toString() << ';' << endl;
        sprintf( numberString, "TREE%d", j + 3);
        consensusFileParameters(numberString)["Parameter name"] = numberString;
        consensusFileParameters(numberString)["Horizontal branch colors"] = "yes";
        consensusFileParameters(numberString)["Branch labels"] =  "yes";
    }
    
}

double HeterogeneousConsensusTree::getValue( Cluster* child, unsigned int i ){
    map<unsigned int, unsigned int>::const_iterator iter = ((HeterogeneousCluster*)child)->getModelCount().begin();
    unsigned int max = 0;
    unsigned int idmax = 0;
    unsigned int tot = 0;
    switch(i){
        case 2: return ((HeterogeneousCluster*)child)->getAverageSubstitutionRate()/
               (double)((LengthCluster*)child)->getNumber();
        case 0: case 1:
            while (iter != ((HeterogeneousCluster*)child)->getModelCount().end()){
                tot += (*iter).second;
                if (max<(*iter).second){
                    idmax = (*iter).first;
                    max = (*iter).second;
                }
                ++iter;
            }
            if (i == 0){
                return idmax;
            }
            else{
                return (double)max/(double)tot;
            }
        default: return (((HeterogeneousCluster*)child)->getParameters())[i-3]/
                (double)((LengthCluster*)child)->getNumber();        
    }
}
