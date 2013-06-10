#include "MLphase.h"

#include <iostream>
#include <iomanip>
#include <fstream>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/randombox.h"
#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include "Models/Model.h"
#include "Tree/SearchTree.h"

#include "Sequence/SequenceTable.h"


using namespace std;


int MLphase::run( int argc, char * * argv ) {

    Phase::run( argc, argv );

    if ( argc != 2 ) {
        cerr << "usage : " << argv[0] << " control_file" << endl;
        exit(EXIT_FAILURE);
    }

    // time structure to record the beginning and the end of the simulation
    time_t start;
    time_t end;

    //declare output files
    ofstream outputTreeFile;
    ofstream outputModelFile;
    ofstream outputFile;

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
    string outputFileName = parameters.stringParameter( "Output file" );
    FileParser::confirmOpenFile( outputTreeFile,
                                 outputFileName + string(".tre"), true );
    FileParser::confirmOpenFile( outputModelFile,
                                 outputFileName + string(".mod"), true );
    FileParser::confirmOpenFile( outputFile,
                                 outputFileName + string(".out"), true );

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
    bool optimizeModel;
    if ( parameters("MODEL").findParameter("Starting model parameters file") ){
        model->initialisation( NULL );
        // Create a FileParser from the model parameters-file
        string startingModelFileName = parameters("MODEL").stringParameter("Starting model parameters file");
        FileParser * modelFileParser = new FileParser( startingModelFileName );
        // Retrieve the model parameters
        ParametersSet * modelParameters = modelFileParser->retrieveParametersSet();
        delete modelFileParser;
        if (!modelParameters){
            cerr << "The file " << startingModelFileName << " is not valid" << endl;
            exit(EXIT_FAILURE);
        }
        //initialise the model with them
        model->setModelParameters(*modelParameters);
        model->validChange();
        delete modelParameters;
        if (parameters("MODEL").findParameter("Optimize model parameters")){
            optimizeModel = parameters("MODEL").boolParameter("Optimize model parameters");
        }
        else{
            optimizeModel = true;
        }
    }
    else{
        model->initialisation( sequenceTable );
        //if initialisation from the sequence, the model should be optimized
        if ( parameters("MODEL").findParameter("Optimize model parameters") &&
             !parameters("MODEL").boolParameter("Optimize model parameters") ){
            //leave it just a warning, it can be case for AA model or JC69
            cerr << "WARNING: if the model is not initialized with proper model "
                 << "parameters (using the field \"Starting model parameters "
                 << "file\"), you probably should optimize the model parameters. "
                 << "You can safely ignore this warning if using an amino acid "
                 << "model, or the JC69 nucleotide model."
                 << endl;
            //cerr << "Press 'C' and return to confirm..." << endl;
            //char c;
            //do{
            //    cin >> c;
            //} while (c!='C');
            optimizeModel = false;
        }
        else{
            optimizeModel = true;
        }
    }

	bool empiricalFreqs;
    if ( parameters("MODEL").findParameter("Empirical frequencies") ){
		empiricalFreqs = parameters("MODEL").boolParameter("Empirical frequencies");
	}
	else{
		empiricalFreqs = false;
	}

    cout.setf(ios::fixed);
    cout << setprecision(4);
    outputFile.setf(ios::fixed);
    outputFile << setprecision(4);

    /** ************************************************************************
     * Tree Initialization
     ************************************************************************ */
    Singleton < Factory<SearchTree> > & treeFactory = Singleton < Factory<SearchTree> >::instance();
    SearchTree* tree = treeFactory.create( parameters( "TREE" ).stringParameter( "Search algorithm" ),
         parameters( "TREE" ) );

    /** ************************************************************************
     * Initial model state initialisation, best result initialisation
     ************************************************************************ */
    vector<double> initialModelParameters;
    if (optimizeModel){
        model->getAllParameters(initialModelParameters);
        tree->optimizeFlag(true);
    }
    else{
        tree->optimizeFlag(false);
    }
    tree->empiricalFreqsFlag(empiricalFreqs);
    if (parameters.findCategory("PENALTY_MODEL")){
        model->initialiseML(parameters("PENALTY_MODEL"));
    }

    start = time( 0 );

    cout << "Number of taxa  = " << sequenceTable->getNumberSpecies() << endl ;
    cout << "Sequence length = " << sequenceTable->getInitialLength() << endl ;
    cout << "Model Name      = " << model->getName() << endl ;
    cout << "ML search started" << endl;
    bool first = tree->initialisation( sequenceTable, model );
    if (first){
        double loglik = tree->loglikelihood();
        cout << "initial stage" << endl;
        cout << tree->toString() << ';' << endl;
        cout << "ML=" << loglik << endl;
        cout << "********************************************************************************" << endl;
        outputFile << "initial stage" << endl;
        outputFile << tree->toString() << ';' << endl;
        outputFile << "ML=" << loglik << endl;
        outputFile << "********************************************************************************" << endl;
    }
    unsigned int stage = 0;
    while(tree->processNext()){
        ++stage;
        double loglik = tree->loglikelihood();
        cout << "*************************" << endl;
        cout << "stage " << stage << endl;
        cout << tree->toString() << ';' << endl;
        cout << "ML=" << loglik << endl;
        cout << "********************************************************************************" << endl;
        outputFile << "stage " << stage << endl;
        outputFile << tree->toString() << ';' << endl;
        outputFile << "ML=" << loglik << endl;
        outputFile << "********************************************************************************" << endl;
    }
    end = time( 0 );
    cout << "Results:" << endl;
    tree->printResults( cout );
    cout << "Time Taken = " << ( end - start ) << endl;
    outputFile << "Results:" << endl;
    tree->printResults( outputFile );
    outputFile << "Time Taken = " << ( end - start ) << endl;
    outputFile.close();

    outputTreeFile.setf(ios::fixed);
	tree->printTree( outputTreeFile );
    //outputTreeFile << bestTree << ";" << endl;
    outputTreeFile.close();

    outputModelFile.setf(ios::fixed);
	ParametersSet bestModelParameters = model->getModelParameters();
    bestModelParameters.saveToFile( outputModelFile );
    outputModelFile.close();

    if (!parameters.checkAllUsed()){
        cerr << "(Double-check that these parameters were intended to be ignored.)" << endl;
    }

    return (EXIT_SUCCESS);
}


int main( int argc, char * argv[] ) {
    MLphase mlphase;
    int res = mlphase.run( argc, argv );
    return res;
}

