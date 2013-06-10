#include "Simulate.h"

// Molecular sequence table
#include "Sequence/SequenceTable.h"

// Phylogenetic Tree
#include "Tree/SimulationTree.h"
#include "Tree/InferenceTree.h"

// Statistics library
#include "Util/statlib.h"

// Util
#include "Util/array2D.h"
#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

// model prototype
#include "Models/Model.h"
#include "PatternDesign/Factory.h"

#include <math.h>
#include <time.h>
#include <assert.h>

//STL streams and algo
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <algorithm>
#include <numeric>

using namespace std;


int Simulate::run( int argc, char * argv[] ) {

    Phase::run( argc, argv );
    
    if (argc != 2){
        cerr << "usage: " <<  argv[0] << " simulate_control_file" << endl;
        exit(EXIT_FAILURE);
    }
    
    // Retrieve the control-file parameters
    FileParser * fileParser = new FileParser( argv[1] );
    ParametersSet* ret = fileParser->retrieveParametersSet();
    if (!ret){
        cerr << "The file " << argv[1] << " is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    ParametersSet parameters = *ret;
    delete fileParser;
    delete ret;
    
    // substitution model creation
    Singleton < Factory< Model > > & modelFactory = Singleton < Factory< Model > >::instance();
    Model* model = modelFactory.create(
                      parameters( "MODEL" ).stringParameter( "Model" ),
                      parameters( "MODEL" ) );
    model->initialisation();

    //in both use of simulate, a model file is either parsed or written
    string modelFileName = parameters.stringParameter("Model parameters file");

    /** *************************************************************************
     * simulate first use : retrieve the name of the model parameters with the
     * given construction parameter
     ************************************************************************** */
    if ( parameters.boolParameter("Retrieve the name of the model's parameters") ){
        ofstream modelParametersInstance;
        FileParser::confirmOpenFile( modelParametersInstance, modelFileName, true );
        (model->getModelParameters()).saveToFile(modelParametersInstance);
        modelParametersInstance.close();
        cout << "An instance of a model parameters file has been written,"
             << " file name : " << modelFileName << endl;
        exit(EXIT_SUCCESS);
    }

    
    /** *************************************************************************
     * simulate second use : simulate sequences
     ************************************************************************** */
    // random seed
    time_t seed;
    if ( parameters.findParameter( "Random seed" ) ){
        seed = parameters.intParameter("Random seed");
    }
    else{
        seed = time( 0 );
        cout << "Random seed not found, CPU time used: seed=" << seed << endl;
    }
    //initialise the randombox
    Singleton < randombox > & rand = Singleton < randombox >::instance();
    rand.setSeed( seed );

    // initialise the model according to the model parameters file
    fileParser = new FileParser( modelFileName );
    ret = fileParser->retrieveParametersSet();
    if (!ret){
        cerr << "The file " << modelFileName << " is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    ParametersSet modelParameters = *ret;
    delete fileParser;
    delete ret;

    model->setModelParameters( modelParameters );
    model->validChange();
    
    // Initialise the tree
    SimulationTree* tree = new SimulationTree( parameters( "TREE" ) );
    
    //prepare the stream in the output file
    string sequencesFileName = parameters.stringParameter("Output file");
    ofstream sequencesOutput;
    FileParser::confirmOpenFile( sequencesOutput, sequencesFileName, true );

    sequencesOutput << "#Model parameters file: " << modelFileName << endl;
    sequencesOutput << "#Random seed: " << seed << endl;
    sequencesOutput << "#Tree: " << tree->toString( false ) << endl;
    sequencesOutput << endl << endl;
    
    
    vector< unsigned int > length;
    vector< string > structure;
    char tmpString1[50];
    char tmpString2[50];
    sprintf(tmpString1,"Number of symbols from class ");
    sprintf(tmpString2,"Structure for the elements of class ");
    char * pchar1 = tmpString1;
    char * pchar2 = tmpString2;
    while (*(++pchar1));
    while (*(++pchar2));
    //get the number of symbols for each model
    unsigned int numberModels = 0;
    bool contin = true;
    while( contin ){
        sprintf(pchar1,"%d", numberModels+1);
        if ( parameters.findParameter(tmpString1) ){
            ++numberModels;
            sprintf(pchar2,"%d", numberModels);
            length.push_back((unsigned int)(parameters.intParameter(tmpString1)));
            if ( parameters.findParameter(tmpString2) ){
                structure.push_back(parameters.stringParameter(tmpString2));
            }
            else{
                structure.push_back("");
            }
        }
        else{
            contin = false;
        }
    }

    if ( numberModels!=model->getNumberSymbolCategory() ){
        cerr << "Error, your model contains " << model->getNumberSymbolCategory()
             << " class(es) of symbol." << endl;
        cerr << "Your control file should consequently declare the fields "
             << "\"Number of symbols from class x\" for all x between 1 and "
             << model->getNumberSymbolCategory() << '.' << endl;
        exit(EXIT_FAILURE);
    }

    // introduction to the file to be rewritten
    sequencesOutput << tree->getNumberTips() << " ";
    bool complete = true;
    vector<unsigned int> rawLength;
    sprintf(tmpString1, "Number of nucleotides from class ");
    pchar1 = tmpString1;
    while (*(++pchar1));
    for (unsigned int i = 0; complete && (i<numberModels); ++i){
        sprintf(pchar1,"%d", i+1);
        if ( parameters.findParameter(tmpString1) ){
            rawLength.push_back(parameters.intParameter(tmpString1));
        }
        else{
            complete = false;
            rawLength.clear();
        }
    }
    //if all the lengths (in nucleotides) are there then output total length
    if (complete){
        //sum the elements of rawLength to have the total length
        sequencesOutput << accumulate( rawLength.begin(), rawLength.end(), 0.0 );
    }
    else{
        cout << "One or more parameters \"Number of nucleotides from class x\" are missing, "
             << "total length of alignment left blank in the sequence file" << endl;
        sequencesOutput << "??fill_number_nucleotides??";
    }
    sequencesOutput << "  ";
    string dataFileType;
    if ( parameters.findParameter("Data file type") ){
        dataFileType = parameters.stringParameter("Data file type");
        sequencesOutput << dataFileType;
    }
    else{
        complete = false;
        cout << "parameters \"Data file type\" not found, left blank in the sequence file" << endl;
        sequencesOutput << "??fill_sequence_type??";
    }
    //a sequence with STRUCT is complete only if the structure was properly provided
    complete = complete && ( (dataFileType!="STRUCT") ||
           (find(structure.begin(),structure.end(),"")==structure.end()) );
    
    
    if (!parameters.checkAllUsed()){
        cerr << "Press 'C' and return to confirm..." << endl;
        char c;
        do{
            cin >> c;
        } while (c!='C');
    }
        
    tree->loadModel( model );
    tree->simulate( length );
    

    sequencesOutput << endl << endl << "#structure" << endl;
    for( unsigned int i = 0; i < numberModels; ++i ){
        for( unsigned int j = 0; j < length[i]; ++j ){
            sequencesOutput << structure[i];
        }
        sequencesOutput << endl;
    }
    sequencesOutput << endl;
    
    tree->printSequenceTable( sequencesOutput );
    sequencesOutput << endl;
    
    //to output the model line we need rawLength to be filled
    if (!rawLength.empty()){
        for( unsigned int i = 0; i < numberModels; ++i ){
            for( unsigned int j = 0; j < rawLength[i]; ++j ){
                 sequencesOutput << i+1 << " ";
            }
            sequencesOutput << endl;
        }
    }
    cout << endl << "Sequences written in " << sequencesFileName << endl;
    
    sequencesOutput.close();

    cout.setf(ios::fixed);
    cout << setprecision(4);
    if (complete){
        cout << "Optional parameters present, parse the created sequence file and compute likelihood" << endl;
        InferenceTree inferenceTree( tree->toString(false) );
        ParametersSet sequencesFileParameter("DATAFILE");
        sequencesFileParameter["Data file"] = sequencesFileName;
        sequencesFileParameter["Interleaved data file"] = "no";
        sequencesFileParameter["Heterogeneous data models"] = "yes";
        SequenceTable* seqTable = new SequenceTable( sequencesFileParameter );
        inferenceTree.loadDataAndModel( seqTable, model );
        assert( inferenceTree.toString( false ) == tree->toString( false ) );
        inferenceTree.invalidateNodes( -1 );
        cout << "LnLikelihood = " << inferenceTree.loglikelihood() << endl;
    }
    return(EXIT_SUCCESS);
}

int main( int argc, char * argv[] ) {
    Simulate sim;
    return sim.run( argc, argv );
}




