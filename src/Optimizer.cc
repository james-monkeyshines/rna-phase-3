#include "Optimizer.h"

#include "Util/Optimise.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <ctype.h>
#include <functional>
#include <float.h>


#include "Tree/OptimizerTree.h"
#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"
#include "Util/randombox.h"
#include "Util/FileParser.h"
#include "Util/ParametersSet.h"
#include "Models/Model.h"

class SequenceTable;

using namespace std;

int Optimizer::run( int argc, char * * argv ) {

    Phase::run( argc, argv );

    if ( argc != 2 ) {
        cerr << "usage : " << argv[0] << " control_file" << endl;
        exit(EXIT_FAILURE);
    }
    
    //declare output files
    ofstream outputTreeFile;
    ofstream outputModelFile;
    ofstream outputFile;

    // Create a FileParser from the control-file
    FileParser * fileParser = new FileParser( argv[1] );

    // Retrieve the control-file parameters
    ParametersSet * ret = fileParser->retrieveParametersSet();
    if (!ret){
        cerr << "The file " << argv[1] << " is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    ParametersSet parameters = *ret;
    delete fileParser;
    delete ret;


    time_t seed;

    if ( parameters.findParameter( "Random seed" ) ){
        seed = parameters.intParameter("Random seed");
    }
    else{
        seed = time( 0 );
        cout << "Random seed not found, CPU time used: seed=" << seed << endl;
    }
    //initialise the randombox, it will be used unless initial model parameters
    //and initial branch lengths are provided
    Singleton < randombox > & randBox = Singleton < randombox >::instance();
    randBox.setSeed(seed);

    /** ************************************************************************
     * Data File Initialization
     ************************************************************************ */
    // Read in molecular data
    SequenceTable * table = new SequenceTable( parameters( "DATAFILE" ) );

    /** ************************************************************************
     * Model Initialization
     ************************************************************************ */
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model> >::instance();
    Model* model = modelFactory.create( parameters( "MODEL" ).stringParameter( "Model" ),
         parameters( "MODEL" ) );
    //model optimisation...
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
        model->initialisation( table );
        //if initialisation from the sequence, the model should be Optimized
        if ( parameters("MODEL").findParameter("Optimize model parameters") &&
             !parameters("MODEL").boolParameter("Optimize model parameters") ){
            //leave it just a warning, it can be case for AA model or JC69
            cerr << "WARNING: if the model is not initialized with proper model "
                 << "parameters (using the field \"Starting model parameters "
                 << "file\"), you probably should Optimize the model parameters. "
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

    cout.setf(ios::fixed);
    cout << setprecision(4);
    outputFile.setf(ios::fixed);
    outputFile << setprecision(4);

    /** ************************************************************************
     * Tree Initialization
     ************************************************************************ */
    Singleton < Factory<OptimizerTree> > & treeFactory = Singleton < Factory<OptimizerTree> >::instance();
    string userTreeFileName = parameters("TREE").stringParameter( "Tree file" );
    ifstream userTreeFile( userTreeFileName.c_str() );


    /** ************************************************************************
     * Initial model state initialisation, best result initialisation
     ************************************************************************ */
    vector<double> initialModelParameters;
    bool optimizeModelPenalty = false;
    vector<double> initialModelPenaltyParameters;
    if (optimizeModel){
        model->getAllParameters(initialModelParameters);
    }
    if (parameters.findCategory("PENALTY_MODEL")){
        model->initialiseML(parameters("PENALTY_MODEL"));
    }
    

    model->getAllPenaltyParameters(initialModelPenaltyParameters);
    optimizeModelPenalty = (model->getLnPenalty() != 0.0) || (initialModelPenaltyParameters.size()!=0);
    double bestLikelihood = -DBL_MAX/2.5;
    double bestModelPenalty = -DBL_MAX/2.5;
    ParametersSet bestModelParameters = model->getModelParameters();
    vector<double> bestModelPenaltyParameters;
    string bestTree;
	int freeParams;
	int branches;
	double lnLAdjEq;
	double lnLAdjEmp;
	int modelSize;
    unsigned int index = 0;
    while(!userTreeFile.eof()){
        ++index;
        cout << "Tree " << index << ':' << endl;
        outputFile << "Tree " << index << ':' << endl;
        string stringTree;
        FileParser::readTree(userTreeFile, stringTree);
        if ( stringTree.length() == 0 ){
            cerr << "Error while reading the tree in the file : "
                 << userTreeFileName << endl;
            exit(EXIT_FAILURE);
        }
        parameters("TREE")["Tree"] = stringTree;
        OptimizerTree* tree = treeFactory.create( "ML optimizer tree", parameters( "TREE" ) );

        if (optimizeModel){
            model->setAllParameters(initialModelParameters);
        }
        if (optimizeModelPenalty){
            model->setAllPenaltyParameters(initialModelPenaltyParameters);
        }
        if ( optimizeModel || optimizeModelPenalty ){
            model->validChange();
        }

        // Load the model and data
        tree->initialisation(table, model);

        cout << tree->toString() << endl;
        cout << "initial likelihood = " << tree->loglikelihood();
        if (optimizeModelPenalty){
             cout << "      initial penalty = " <<  model->getLnPenalty();
        }
        cout << endl;

        Optimise::optimiseQuasiNewton( *tree, *model, optimizeModel, empiricalFreqs,
                optimizeModelPenalty ? 1e-9 : 1e-8 , true );

        double maxLikelihood = tree->loglikelihood();
        double maxModelPenalty = 0.0;
        if (optimizeModelPenalty){
            maxModelPenalty = model->getLnPenalty();
        }
        if ( maxLikelihood + maxModelPenalty > bestLikelihood + bestModelPenalty ){
            bestLikelihood = maxLikelihood;
            bestModelPenalty = maxModelPenalty;
            bestTree = tree->toString(false);
			branches = tree->getNumberBranches();
			freeParams = model->getNumberFreeParameters();
			lnLAdjEq = model->getlnLAdjustmentEq();
			lnLAdjEmp = model->getlnLAdjustmentEmp();
			modelSize = model->getNumberSymbolCategory();
			if (optimizeModel){
                bestModelParameters = model->getModelParameters();
            }
            if (optimizeModelPenalty){
                model->getAllPenaltyParameters(bestModelPenaltyParameters);
            }
        }

        if (optimizeModel){
            model->printParameters( cout );
            model->printParameters( outputFile );
        }
        cout << tree->toString(false) << ";" << endl;
        outputFile << tree->toString(false) << ";" << endl;
        
        cout.setf(ios::fixed);
        cout << setprecision(4);
        outputFile.setf(ios::fixed);
        outputFile << setprecision(4);
        if(optimizeModelPenalty){
            cout << "optimized value=" << maxLikelihood + maxModelPenalty << "     (";
            outputFile << "optimized value=" << maxLikelihood + maxModelPenalty << "     (";
        }
        cout << "maxLikelihood = " << maxLikelihood;
        outputFile << "maxLikelihood = " << maxLikelihood;
        if(optimizeModelPenalty){
            cout << ", maxPenalty=" << maxModelPenalty << ')';
            outputFile << ", maxPenalty=" << maxModelPenalty << ')';
        }
        cout << endl << "-----------------------" << endl;
        outputFile << endl << "-----------------------" << endl;
        delete tree;
        userTreeFile >> ws;
    }
    cout << endl << endl;
    outputFile << endl << endl;
    if(optimizeModelPenalty){
        cout << "best result=" << bestLikelihood + bestModelPenalty << "     (";
        outputFile << "best result=" << bestLikelihood + bestModelPenalty << "     (";
    }
	cout << "free parameters = model + branches = " << freeParams << " + " << branches << " = " << freeParams+branches << endl;
    outputFile << "free parameters = model + branches = " << freeParams << " + " << branches << " = " << freeParams+branches << endl;
	if (lnLAdjEq <= 0 & modelSize > 1) {
		cout << "likelihood adjustment (equal frequencies): " << lnLAdjEq << endl;
		cout << "likelihood adjustment (empirical frequencies => +9 free parameters): " << lnLAdjEmp << endl;
		outputFile << "likelihood adjustment (equal frequencies): " << lnLAdjEq << endl;
		outputFile << "likelihood adjustment (empirical frequencies => +9 free parameters): " << lnLAdjEmp << endl;
	}
    cout << "best likelihood = " << bestLikelihood;
    outputFile << "best likelihood = " << bestLikelihood;
    if(optimizeModelPenalty){
        cout << ", best penalty=" << bestModelPenalty << ')';
        outputFile << ", best penalty=" << bestModelPenalty << ')';
    }
    cout << endl << bestTree << ";" << endl;
    outputFile << endl << bestTree << ";" << endl;
    outputFile.close();
    //save the best tree into its file...
    outputTreeFile.setf(ios::fixed);
    outputTreeFile << bestTree << ";" << endl;
    outputTreeFile.close();
    //... and the best substitution model
    outputModelFile.setf(ios::fixed);
    bestModelParameters.saveToFile( outputModelFile );
    outputModelFile.close();

    //check whether all parameters and categories have been used to issue some warnings.
    if (!parameters.checkAllUsed()){
        cerr << "(Double-check that these parameters were intended to be ignored.)" << endl;
    }

    delete table;
    delete model;

    return (EXIT_SUCCESS);
}

int main( int argc, char* argv[]){
    Optimizer optimizer;
    return optimizer.run( argc, argv );
}
