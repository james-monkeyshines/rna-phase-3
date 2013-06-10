#include "DistPhase.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/randombox.h"
#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include "Models/Model.h"
#include "Tree/DistTree.h"

#include "Sequence/SequenceTable.h"


using namespace std;


int DistPhase::run( int argc, char * * argv ) {

    Phase::run( argc, argv );

    if ( argc != 2 ) {
        cerr << "usage : " << argv[0] << " control_file" << endl;
        exit(EXIT_FAILURE);
    }

    // Create a FileParser from the control-file
    FileParser * fileParser = new FileParser( argv[1] );

    // Retrieve the control-file parameters
    ParametersSet* ret = fileParser->retrieveParametersSet();
    delete fileParser;

    //check whether the control file was properly opened
    if (!ret){
        cerr << "The file " << argv[1] << " is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    ParametersSet parameters = *ret;
    delete ret;


    /** ************************************************************************
     * Output File Initialization
     ************************************************************************ */
    ofstream outputFile;
    string outputFileName = parameters.stringParameter( "Output file" );
    FileParser::confirmOpenFile( outputFile, outputFileName, true );

    /** initialise the random number generator */
    Singleton<randombox>& randBox = Singleton<randombox>::instance();
    int seed;
    if (parameters.findParameter("Random seed") ){
        seed = parameters.intParameter("Random seed");
    }
    else{
        seed = time(0);
        cout << "Random seed not found, CPU time used: seed=" << (int)seed << endl;
        outputFile << "Random seed not found, CPU time used: seed=" << (int)seed << endl;
    }
    randBox.setSeed(seed);

    /** ************************************************************************
     * Data file Initialization
     ************************************************************************ */
    SequenceTable * sequenceTable =  new SequenceTable( parameters( "DATAFILE" ) );


    /** ************************************************************************
     * Model Initialization
     ************************************************************************ */
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model>  >::instance();
    Model * model = modelFactory.create( parameters( "MODEL" ).stringParameter( "Model" ),
         parameters( "MODEL" ) );

    //initialisation of the model
    model->initialisation( NULL );
    // Create a FileParser from the model parameters-file
    string modelFileName = parameters("MODEL").stringParameter("Model parameters file");
    FileParser * modelFileParser = new FileParser( modelFileName );
    // Retrieve the model parameters
    ParametersSet * modelParameters = modelFileParser->retrieveParametersSet();
    delete modelFileParser;
    if (!modelParameters){
        cerr << "The file " << modelFileName << " is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    //initialise the model with them
    model->setModelParameters(*modelParameters);
    model->validChange();
    delete modelParameters;

    /** ************************************************************************
     * Tree Initialization
     ************************************************************************ */
    ParametersSet treeParameters("TREE");
    DistTree* distTree = new DistTree( treeParameters );


    /** ************************************************************************
     * output format
     ************************************************************************ */
    DistTree::MatrixFormat frmt = DistTree::SQUARE;
    if (parameters.findParameter("Output format")){
        string format = parameters.stringParameter("Output format");
        if ( format == "lower-triangular" ){
            frmt=DistTree::LOWER;
        }
        else{
            if ( format == "upper-triangular" ){
                frmt=DistTree::UPPER;
            }
        }
        //frmt still equal to SQUARE?
        if (frmt==DistTree::SQUARE){
            if ( format != "square" ){
                cerr << "Unrecognized output format: \"" << format << "\". Use \"lower-triangular\", \"upper-triangular\" or \"square\" (case-sensitive)." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }

    if (!parameters.checkAllUsed()){
        cerr << "Press 'C' and return to confirm..." << endl;
        char c;
        do{
            cin >> c;
        } while (c!='C');
    }

    cout << "Number of taxa  = " << sequenceTable->getNumberSpecies() << endl ;
    cout << "Sequence length = " << sequenceTable->getInitialLength() << endl ;
    cout << "Model Name      = " << model->getName() << endl;
    
    //init
    distTree->initialisation( sequenceTable, model );
    
    //go!
    distTree->computePairwiseDist();
    
    
    
    cout.setf(ios::fixed);
    cout << setprecision(4);
    cout << "Results:" << endl;
    distTree->printResults( cout, frmt );
    
    outputFile.setf(ios::fixed);
    outputFile << setprecision(4);
    distTree->printResults( outputFile, frmt  );
    outputFile.close();
    
    return (EXIT_SUCCESS);
}


int main( int argc, char * argv[] ) {
    DistPhase distPhase;
    int res = distPhase.run( argc, argv );
    return res;
}

