#include "MCMCphase.h"

#include <algorithm>
#include <assert.h>
/*
#include <dlfcn.h>
*/

#include <iostream>
#include <iomanip>
#include <iterator>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/randombox.h"
#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include "Models/Model.h"

#include "Tree/MCMCTree.h"

using namespace std;

void MCMCphase::runMCMC( MCMCTree * tree, Model * model,
unsigned int perturbationTree, unsigned int perturbationModel,
int burninCycles, int samplingCycles, int samplingPeriod,
int currentIteration, string terminator ) {

    Singleton < randombox > & randBox = Singleton < randombox >::instance();

    double oldLnLikelihood = tree->loglikelihood();
    double newLnLikelihood = oldLnLikelihood;

     //IMPORTANT: we will assume that modification in the model do not affect the prior of the tree
     //Consequently, the tree is responsible to compute the possible priors related to the variation
     //of the model across the tree
     double oldLnTreePrior = tree->getLnPrior();
     double newLnTreePrior = oldLnTreePrior;
     double oldLnModelPrior = model->getLnPrior();
     double newLnModelPrior = oldLnModelPrior;

    //if this is beginning of the run or after a restore, we have to write
     //the current tree or the content of the bestParams vector as the best tree
    bool newBest = true;

    //storage vector to fill the rescue file
    vector< double > modelParams;
    vector< double > treeParams;
    vector< double > modelPertParams;
    vector< double > treePertParams;
    vector< double > modelPriorParams;
    vector< double > treePriorParams;
    vector< double > randomBoxParams;


    /**
     * save the first state
     * we save a state when it is accepted and we come
     * back to it before each perturbation
     */
    //tree->saveNodesAux();

    model->attach( *tree );

    tree->saveNodesAux();

    double lnHastingRatio = 0.0;
    double lnProb;

   // int perturbCat;
   // int perturbModel;
    time_t t;

    //precision 4 for the output of the likelihood during the run
    cout.setf(ios::fixed);
    cout << setprecision(4);

    double probPerturbationTree = (double)perturbationTree/
                              (double)(perturbationTree+perturbationModel);
    //launch the chain for burninCycles + samplingCycles
    for ( int cycle = currentIteration; cycle < ( burninCycles + samplingCycles ); ++cycle ){

        assert (oldLnLikelihood == newLnLikelihood);
        assert (oldLnTreePrior == newLnTreePrior);
        assert (oldLnModelPrior == newLnModelPrior);

        //if the burnin period is over, warn the model (to stop the
        //dynamical modification of the perturbation steps)
        if ( cycle == burninCycles){
            model->stopBurn();
            tree->stopBurn();
            acceptedModelPerturbation = 0;
            numberModelPerturbation = 0;
            acceptedTreePerturbation = 0;
            numberTreePerturbation = 0;
        }

#ifdef DEBUG3
        float testPrior = (float)tree->getLnPrior();
        if ((float)oldLnTreePrior != testPrior){
            cerr << "tree prior saved value = " << oldLnTreePrior << endl;
            cerr << "tree prior returned value = " << testPrior << endl;
            exit(EXIT_FAILURE);
        }
        testPrior = (float)model->getLnPrior();
        if ((float)oldLnModelPrior != testPrior){
            cerr << "model prior saved value = " << oldLnModelPrior << endl;
            cerr << "model prior returned value = " << testPrior << endl;
            exit(EXIT_FAILURE);
        }
        tree->checkIdentical();
        tree->checkNodesCount();
        tree->checkIndex();
        tree->checkClusters();
        tree->invalidateNodes(-1);
        if (fabs(newLnLikelihood - tree->loglikelihood()) > 1e-6 ){
                        exit(EXIT_FAILURE);
        }
#endif
        /** tree perturbation *************************************/
        if ( randBox.ran() < probPerturbationTree ) {
            lnHastingRatio = tree->perturb();
            ++numberTreePerturbation;
            newLnLikelihood = tree->loglikelihood();
            //get new tree prior
            newLnTreePrior = tree->getLnPrior();
#ifdef DEBUG3
            tree->checkClusters();
            tree->checkIdentical();
            tree->checkIndex();
            tree->invalidateNodes(-1);
            assert(newLnLikelihood == tree->loglikelihood());
            //model prior should not be affected by changes in the tree
            assert(oldLnModelPrior == model->getLnPrior());
#endif
            //prob acceptance = MIN ( 1.0,  PP(newState) / PP(oldState) )
            lnProb = newLnLikelihood - oldLnLikelihood + lnHastingRatio + newLnTreePrior - oldLnTreePrior;
            // Accept new state
            // WARNING: experiment cannot be repeated if lnProb is (too often) close to 0.
            // Because the call too randBox.ran() will become random (tiny numerical imprecision)
            if ( (lnProb>=0.0) || (log(randBox.ran()) < lnProb) ){
                tree->validatePerturbation( true );
                ++acceptedTreePerturbation;
                oldLnLikelihood = newLnLikelihood;
                oldLnTreePrior = newLnTreePrior;
                if ( newLnLikelihood > bestLikelihood ) {
                    newBest = true;
                    bestStringTree = tree->toString( false );
                    model->getAllParameters( bestModelParams );
                    bestLikelihood = newLnLikelihood;
                }
            }
            // Reject new state
            else {
                //?usually we do not revert to the old state with this call?
                tree->validatePerturbation( false );
                newLnLikelihood = oldLnLikelihood;
                newLnTreePrior = oldLnTreePrior;
            }
        }
        /** model perturbation *************************************/
        else{
            lnHastingRatio = model->perturb();
            //part of the saved computation was invalidated
            //tree->invalidateNodesAux(  );
            //should be automatic now
            ++numberModelPerturbation;
            newLnLikelihood = tree->loglikelihood();
            //get new priors (the prior on the tree might have changed if there is a prior on model variation across the tree)
            newLnModelPrior = model->getLnPrior();
            newLnTreePrior = tree->getLnPrior();
            //prob acceptance = MIN ( 1.0,  PP(newState) / PP(oldState) )
            lnProb = newLnLikelihood - oldLnLikelihood + lnHastingRatio + newLnModelPrior - oldLnModelPrior + newLnTreePrior - oldLnTreePrior;
            // Accept new state
            // WARNING: experiment cannot be repeated if lnProb is (too often) close to 0.
            // Because the call too randBox.ran() will become random (tiny numerical imprecision)
            if ( (lnProb>=0.0) || (log(randBox.ran()) < lnProb) ){
                model->validatePerturbation( true );
                ++acceptedModelPerturbation;
                oldLnLikelihood = newLnLikelihood;
                oldLnTreePrior = newLnTreePrior;
                oldLnModelPrior = newLnModelPrior;
                //save the new state (computation)
                tree->saveNodesAux();
                if ( newLnLikelihood > bestLikelihood ) {
                    newBest = true;
                    bestStringTree = tree->toString( false );
                    model->getAllParameters( bestModelParams );
                    bestLikelihood = newLnLikelihood;
                }
            }
            // Reject new state
            else {
                model->validatePerturbation( false );
                tree->retrieveNodesAux();
                newLnLikelihood = oldLnLikelihood;
                newLnTreePrior = oldLnTreePrior;
                newLnModelPrior = oldLnModelPrior;
            }
        }

        // output the plot and save the current state
        if ( !((cycle+1)%samplingPeriod) ) {
           /* try to perform a maximum of tasks before saving
            * to minimize the writing time (and avoid corrupt
            * save)                                           */
            model->getAllParameters( modelParams );
            tree->getAllParameters( treeParams );
            model->getAllPerturbationParameters( modelPertParams );
            tree->getAllPerturbationParameters( treePertParams );
            model->getAllPriorParameters( modelPriorParams );
            tree->getAllPriorParameters( treePriorParams );
            t = time(0) - start;
            if (newBest){
                newBest = false;
                /* start writing */
                rescueBest.seekp( 0 );
                for (unsigned int i =0; i < bestModelParams.size(); ++i ){
                    rescueBest.write( (char*)&bestModelParams[i], sizeof(double));
                }
                bestTreeFile.seekp( 0 );
                bestTreeFile << bestStringTree << terminator << endl;
            }
            rescue.seekp( 0 );
            randBox.getState( randomBoxParams );
            for (unsigned int i =0; i < randomBoxParams.size(); ++i ){
                rescue.write( (char*)&randomBoxParams[i], sizeof(double));
            }
            for (unsigned int i =0; i < modelParams.size(); ++i ){
                rescue.write( (char*)&modelParams[i], sizeof(double));
            }
            rescue << tree->toString( true );
            rescue.put(';');
            for (unsigned int i =0; i < treeParams.size(); ++i ){
                rescue.write( (char*)&treeParams[i], sizeof(double) );
            }
            for (unsigned int i =0; i < modelPertParams.size(); ++i ){
                rescue.write( (char*)&modelPertParams[i], sizeof(double));
            }
            for (unsigned int i =0; i < treePertParams.size(); ++i ){
                rescue.write( (char*)&treePertParams[i], sizeof(double));
            }
            for (unsigned int i =0; i < modelPriorParams.size(); ++i ){
                rescue.write( (char*)&modelPriorParams[i], sizeof(double));
            }
            for (unsigned int i =0; i < treePriorParams.size(); ++i ){
                rescue.write( (char*)&treePriorParams[i], sizeof(double));
            }
            rescue.write( (char*)&cycle, sizeof(int));
            rescue.write( (char*)&bestLikelihood, sizeof(double));
            rescue.write( (char*)&acceptedModelPerturbation, sizeof(unsigned int) );
            rescue.write( (char*)&numberModelPerturbation, sizeof(unsigned int) );
            rescue.write( (char*)&acceptedTreePerturbation, sizeof(unsigned int) );
            rescue.write( (char*)&numberTreePerturbation, sizeof(unsigned int) );
            rescue.write( (char*)&t, sizeof(time_t));
            rescue << flush;
            plotfile << newLnLikelihood << ' '
                     << newLnTreePrior << ' ' << newLnModelPrior << endl;
            //sampling
            if (cycle>burninCycles){
                tree->sample();
                model->printLine( mpFile );
                mpFile << endl;
                model->printPriorLine( priorFile );
                priorFile << endl;

            }
        }

#ifdef DEBUG1
        cout << (cycle+1) << ".LnLikelihood = " << newLnLikelihood
                          << "    ;Prior on tree = " << newLnTreePrior
                          << "    ;Prior on subs. model = " << newLnModelPrior << endl;
#endif
#ifdef DEBUG2
        cout << (cycle+1) << ".Tree = " << tree->toString(false) << endl;
        tree->invalidateNodes(-1);
        if (newLnLikelihood!=tree->loglikelihood()){
            cerr << "saved likelihood: " << newLnLikelihood << endl;
            cerr << "computed likelihood: " << tree->loglikelihood() << endl;
        }
        assert (fabs((float)newLnLikelihood-(float)tree->loglikelihood())<1e-6);
#endif

        if ( !((cycle+1)%100) ) {
            cout << (cycle+1) << ".LnLikelihood = " << newLnLikelihood
                     << "    ;Prior on tree = " << newLnTreePrior
                     << "    ;Prior on subs. model = " << newLnModelPrior << endl;
            cout << tree->toString( false ) << endl;
        }
    }

    // Frequency Proposal Phase (End)
    cout.setf(ios::fixed);
    output.setf(ios::fixed);
    cout << setprecision(2);
    output << setprecision(2);
    cout << "Tree acceptance Rate (during sampling period)   =  "
         << (double)acceptedTreePerturbation*100.0 / (double)numberTreePerturbation
         << "%" << endl;
    cout << "Model acceptance Rate (during sampling period)  =  "
         << (double)acceptedModelPerturbation*100.0 / (double)numberModelPerturbation
         << "%" << endl;
    output << "Tree acceptance Rate (during sampling period)   =  "
         << (double)acceptedTreePerturbation*100.0 / (double)numberTreePerturbation
         << "%" << endl;
    output << "Model acceptance Rate (during sampling period)  =  "
         << (double)acceptedModelPerturbation*100.0 / (double)numberModelPerturbation
         << "%" << endl;

    tree->printPerturbationParameters(cout);
    tree->printPerturbationParameters(output);

    cout << setprecision(4);
    output << setprecision(4);
    cout << "Best Loglikelihood = " << bestLikelihood << endl;
    cout << "Best Tree" << endl;
    cout << bestStringTree << endl;
    output << "Best Tree" << endl;
    output << bestStringTree << endl;
    output << "Best Loglikelihood = " << bestLikelihood << endl;

    model->setAllParameters( bestModelParams );
    model->validChange();
    model->printParameters( cout );
    model->printParameters( output );
    model->printPerturbationParameters( cout );
    model->printPerturbationParameters( output );

    model->getModelParameters().saveToFile( bestModelParametersFile );
    bestTreeFile << bestStringTree << terminator << endl;

#ifdef DEBUG2
    InferenceTree checktree( bestStringTree );
    checktree.loadDataAndModel( tree->getTable(), model );
    assert( checktree.toString( true ) == bestStringTree );
    //branch lengths are different (rounding errors)
    //so we expect a different likelihood
    if( fabs(checktree.loglikelihood() - bestLikelihood) > 10 ){
        cerr << "Warning, big likelihood difference found for the saved best state" << endl;
    }

#endif
}



int MCMCphase::run( int argc, char * * argv ) {

    Phase::run( argc, argv );
    
    if ( argc != 2 ) {
        cerr << "usage : " << argv[0] << " control_file" << endl;
        exit(EXIT_FAILURE);
    }

    // random seed
    int seed = 0;

    // sequenceTable object : to read the data file and store the molecular
    // data
    SequenceTable * table = NULL;

    // Create a FileParser from the control-file
    FileParser * fileParser = new FileParser( argv[1] );
    // Retrieve the control-file parameters
    ParametersSet * ret = fileParser->retrieveParametersSet();
    delete fileParser;
    if (!ret){
        cerr << "The file " << argv[1] << " is not valid" << endl;
        exit(EXIT_FAILURE);
    }
    ParametersSet parameters = *ret;
    delete ret;

    string terminatorChar;

    int currentIteration;


   /** ************************************************************************
    * Rescue mode?
    ************************************************************************ */
    string filename;
    string outputFileBaseName( parameters.stringParameter( "Output file" ) );
    filename = outputFileBaseName + string( ".rsb" );
    ifstream rescueBestRead( filename.c_str(), ifstream::in | ifstream::binary);
    filename = outputFileBaseName + string( ".rsc" );
    ifstream rescueRead( filename.c_str(), ifstream::in | ifstream::binary);
    bool overwrite = true;
    string answer;
    if ( rescueRead.is_open() ) {
        cout << "Rescue file found, it seems a previous experiment did "
             << "not finish properly." << endl;
        cout << "Do you want to try to recover it? ";
        cin >> answer  ;
        if ( ( answer[0] == 'y' ) || ( answer[0] == 'Y' ) ) {
            overwrite = false;
        }
        else{
            rescueRead.close();
        }
    }
    
    // Find out if the ".output" file already exists and if it does
    // confirm overwrite
    filename = outputFileBaseName + string( ".out" );
    FileParser::confirmOpenFile( output, filename, overwrite );

   /** ************************************************************************
    * Initialise/restore the random generator
    ************************************************************************ */
    Singleton < randombox > & randBox = Singleton < randombox >::instance();
    //this vector is used to restore the randombox state in case of recovery
    vector<double> randomParams;
    if (rescueRead.is_open()){
        unsigned int nbParams = randBox.getNumberParams();
        randomParams.resize( nbParams );
        for (unsigned int i = 0; i < nbParams; ++i){
            rescueRead.read( (char*)&randomParams[i], sizeof(double) );
        }
        if ( !rescueRead.good() || rescueRead.fail() ){
            cerr << "Sorry, cannot restore the saved state... rescue file corrupted" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        if ( parameters.findParameter( "Random seed" ) ){
            seed = parameters.intParameter("Random seed");
        }
        else{
            seed = (int)time( 0 );
            cout << "Random seed not found, CPU time used: seed=" << seed << endl;
        }
        randBox.setSeed(seed);
    }

   /** ************************************************************************
    * Data File Initialization
    ************************************************************************ */
    cout << "Data file initialisation..." << flush;
    // Read in molecular data
    string dataFile = parameters( "DATAFILE" ).stringParameter( "Data file" );
    table = new SequenceTable( parameters( "DATAFILE" ) );
    // print out molecular data to output file and standard output
    //table->print( output );
    //table->print( cout );
    //output << endl;
    //cout << endl;
    cout << "done" << endl;
   /** ************************************************************************
    * End Data File Initialization
    ************************************************************************ */

   /** ************************************************************************
    * Model Initialization
    ************************************************************************ */
    cout << "Model initialisation..." << flush;
    bool randStartingModel = false;
    string startingModelFileName;
    Model * startingModel = NULL;
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model> >::instance();
    startingModel = modelFactory.create( parameters( "MODEL" ).stringParameter( "Model" ),
         parameters( "MODEL" ) );
    if (rescueRead.is_open()){
        startingModel->initialisation( NULL );
        //before restoring the saved state we want to recover the vector
        //bestModelVector
        unsigned int numberModelParameters = startingModel->getNumberFreeParameters();
        cout << "recovering the best model visited..." << flush;
        bool readOK = true;
        bestModelParams.resize(numberModelParameters);
        for (unsigned int j = 0; j < numberModelParameters && readOK; ++j){
            if ( !rescueBestRead.good() ) readOK = false;
            rescueBestRead.read( (char*)&bestModelParams[j], sizeof(double) );
            if ( rescueBestRead.fail() ) readOK = false;
        }
        if (!readOK){
            bestModelParams.resize(0);
            cerr << "failure, best state lost" << endl;
        }
        //now we restore the current state
        cout << "restoring the saved state..." << flush;
        vector<double> modelParametersRescue;
        modelParametersRescue.resize(numberModelParameters);
        readOK = true;
        for (unsigned int j = 0; j < numberModelParameters && readOK; ++j){
            if ( !rescueRead.good() ) readOK = false;
            rescueRead.read( (char*)&modelParametersRescue[j], sizeof(double) );
            if ( rescueRead.fail() ) readOK = false;
        }
        if (!readOK){
            cerr << "failure, cannot restore the state...sorry." << endl;
            exit(EXIT_FAILURE);
        }
        startingModel->setAllParameters( modelParametersRescue );
        startingModel->validChange();
        parameters( "MODEL" ).touch("Starting model parameters file");
    }
    else{
        randStartingModel = !parameters( "MODEL" ).findParameter("Starting model parameters file");
        if ( randStartingModel ){
            cout << "initialise with empirical sequences..." << flush;
            startingModel->initialisation( table );
        }
        else{
            cout << "initialise with a model file..." << flush;
            startingModel->initialisation( NULL );
            startingModelFileName =
                parameters( "MODEL" ).stringParameter("Starting model parameters file");
            // Create a FileParser from the model parameters-file
            FileParser * modelFileParser = new FileParser( startingModelFileName );
            // Retrieve the model parameters
            ret = modelFileParser->retrieveParametersSet();
            delete modelFileParser;
            if (!ret){
                cerr << "The file " << startingModelFileName << " is not a valid model file" << endl;
                exit(EXIT_FAILURE);
            }
            ParametersSet modelParameters = *ret;
            delete ret;
            //initialise the model with them
            startingModel->setModelParameters(modelParameters);
            startingModel->validChange();
        }
    }
    cout << "done" << endl;


   /** ************************************************************************
    * Tree Initialization
    ************************************************************************ */
    cout << "Tree initialisation..." << flush;
    bool randStartingTree = false;
    string startingTreeFileName;
    Singleton< Factory<MCMCTree> > & treeFactory =
            Singleton< Factory<MCMCTree> >::instance();
    MCMCTree * startingTree =
        treeFactory.create( parameters("TREE").stringParameter( "Tree" ),
                            parameters( "TREE" ));
    if (rescueRead.is_open()){
        cout << "recovering the best tree visited..." << flush;
        ifstream bestTreeF( (outputFileBaseName + string( "-best.tre" )).c_str() );
        if (!bestTreeF.is_open()){
            cerr << "Warning, cannot retrieve the best tree visited" << endl;
            //first the best state (with length)
        }
        else{
            FileParser::readTree(bestTreeF, bestStringTree);
        }
        cout << "restoring the saved state..." << flush;
        string stringTree;
        char c = 0;
        while ( c != ';' ){
            rescueRead.get(c);
            stringTree += c;
        }
        // Load the model and data
        startingTree->loadDataAndModel( table, startingModel );
        startingTree->constructFromString( stringTree );
        //now we restore the current tree (branch lengths are after in the rescue file)
        unsigned int numberTreeParameters = startingTree->getNumberTreeParameters();
        vector<double> treeParametersRescue;
        treeParametersRescue.resize(numberTreeParameters);
        for (unsigned int j = 0; j < numberTreeParameters; ++j){
            rescueRead.read( (char*)&treeParametersRescue[j], sizeof(double) );
        }
        startingTree->setAllParameters( treeParametersRescue );
        parameters("TREE").touch("Starting tree file");
    }
    else{
        // Load the model and data
        startingTree->loadDataAndModel( table, startingModel );
        randStartingTree = !parameters("TREE").findParameter("Starting tree file");
        if ( randStartingTree ){
            cout << "random construction..." << flush;
            startingTree->constructRandomly( .6 );
        }
        else{
            cout << "initialise from file..." << flush;
            startingTreeFileName = parameters("TREE").stringParameter("Starting tree file");
            string stringTree;
            ifstream startingTreeFile(startingTreeFileName.c_str());
            FileParser::readTree(startingTreeFile, stringTree);
            if ( stringTree.length() == 0 ){
                cerr << endl << "Error while reading the tree in the file : "
                     << startingTreeFileName << endl;
                exit(EXIT_FAILURE);
            }
            startingTree->constructFromString( stringTree );
        }
    }
    //Invalidate partial likelihood calculations(to force fresh recalculations)
    //startingTree->invalidateNodesAux();
    cout << "done..." << endl;


   /** ************************************************************************
    * recover burnin/sampling time
    ************************************************************************ */
    int numberBurninIterations = parameters.intParameter("Burnin iterations");
    int numberSamplingIterations = parameters.intParameter("Sampling iterations");
    int samplingPeriod = parameters.intParameter("Sampling period");



   /** ************************************************************************
    * Initialise perturbation and priors
    ************************************************************************ */
    startingModel->initialiseMCMC( parameters("PERTURBATION")("PERTURBATION_MODEL") );
    startingTree->initialiseMCMC( parameters("PERTURBATION")("PERTURBATION_TREE") );
    if (rescueRead.is_open()){
        cout << "restoring transition kernel parameters..." << flush;
        vector<double> modelPerturbationRescue;
        unsigned int numberModelPerturbationParameters =
                startingModel->getNumberPerturbationParameters();
        modelPerturbationRescue.resize(numberModelPerturbationParameters);
        vector<double> treePerturbationRescue;
        unsigned int numberTreePerturbationParameters =
                startingTree->getNumberPerturbationParameters();
        treePerturbationRescue.resize(numberTreePerturbationParameters);

        for (unsigned int j = 0; j < numberModelPerturbationParameters; ++j){
            rescueRead.read( (char*)&modelPerturbationRescue[j], sizeof(double) );
        }
        startingModel->setAllPerturbationParameters( modelPerturbationRescue );
        for (unsigned int j = 0; j < numberTreePerturbationParameters; ++j){
            rescueRead.read( (char*)&treePerturbationRescue[j], sizeof(double) );
        }
        startingTree->setAllPerturbationParameters( treePerturbationRescue );
        cout << "done" << endl;
        cout << "restoring priors..." << flush;
        vector<double> priorRescue;
        unsigned int numberPriorParameters =
                startingModel->getNumberPriorParameters();
        priorRescue.resize(numberPriorParameters);
        for (unsigned int j = 0; j < numberPriorParameters; ++j){
            rescueRead.read( (char*)&priorRescue[j], sizeof(double) );
        }
        startingModel->setAllPriorParameters( priorRescue );

        numberPriorParameters =
                startingTree->getNumberPriorParameters();
        priorRescue.resize(numberPriorParameters);
        for (unsigned int j = 0; j < numberPriorParameters; ++j){
            rescueRead.read( (char*)&priorRescue[j], sizeof(double) );
        }
        startingTree->setAllPriorParameters( priorRescue );
    }

    //restoring the random box state
    if ( randomParams.size() ){
        randBox.setState( randomParams );
        parameters.touch( "Random seed" );
    }

    unsigned int perturbationTree =
        parameters("PERTURBATION").intParameter("Tree, proposal priority");
    unsigned int perturbationModel =
        parameters("PERTURBATION").intParameter("Model, proposal priority");

   /** ************************************************************************
    * Print out parameters to output file
    ************************************************************************ */

    if (!rescueRead.is_open()){
        cout << "Random seed = " << seed << endl;
        output << "Random seed = " << seed << endl;

        cout << "Data file          = " << dataFile << endl;
        cout << "Number of species  = " << table->getNumberSpecies() << endl;
        cout << "Model              = " << startingModel->getName() << endl;
        cout << "Burnin iterations      = " << numberBurninIterations << endl;
        cout << "Sampling iterations    = " << numberSamplingIterations << endl;
        cout << "Sampling period        = " << samplingPeriod << endl;

        output << "Data file          = " << dataFile << endl;
        output << "Number of species  = " << table->getNumberSpecies() << endl;
        output << "Model              = " << startingModel->getName() << endl;
        output << "Burnin iterations      = " << numberBurninIterations << endl;
        output << "Sampling iterations    = " << numberSamplingIterations << endl;
        output << "Sampling period        = " << samplingPeriod << endl;

        cout << "Random start model parameters             = "
             << ( randStartingModel ? "yes" : "no" ) << endl;
        output << "Random start model parameters             = "
               << ( randStartingModel ? "yes" : "no" ) << endl;
        if ( !randStartingModel ) {
            cout << "User's starting model parameters file = "
                 << startingModelFileName << endl;
            output << "User's starting model parameters file = "
                   << startingModelFileName << endl;
        }
        cout << "Random start tree                         = "
             << ( randStartingTree ? "yes" : "no" ) << endl;
        output << "Random start tree                         = "
               << ( randStartingTree ? "yes" : "no" ) << endl;
        if ( !randStartingTree ) {
            cout << "User's starting tree file             = "
                 << startingTreeFileName << endl;
            output << "User's starting tree file             = "
                   << startingTreeFileName << endl;
        }
        cout << "Output files root       = " << outputFileBaseName << endl;
        output << "Output files root       = " << outputFileBaseName << endl;
        // Output the initial tree and its loglikelihood
        cout << "Initial Tree" << endl;
        cout << startingTree->toString( false ) << endl;
        output << "Initial Tree" << endl;
        output << startingTree->toString( false ) << endl;
        cout << "Lnlikelihood = " << startingTree->loglikelihood() << endl;
        startingModel->printParameters( cout );
        output << "Lnlikelihood = " << startingTree->loglikelihood() << endl;
        startingModel->printParameters( output );
    }
    else{
        output << endl << "------------------restart-------------------" << endl;
    }

    // (re)start the timer and load the best likelihood and the current iteration
    if (rescueRead.is_open()){
        rescueRead.read( (char*)&currentIteration, sizeof(int) );
        ++currentIteration;
        rescueRead.read( (char*)&bestLikelihood, sizeof(double) );
        rescueRead.read( (char*)&acceptedModelPerturbation, sizeof(unsigned int) );
        rescueRead.read( (char*)&numberModelPerturbation, sizeof(unsigned int) );
        rescueRead.read( (char*)&acceptedTreePerturbation, sizeof(unsigned int) );
        rescueRead.read( (char*)&numberTreePerturbation, sizeof(unsigned int) );
        rescueRead.read( (char*)&start, sizeof(time_t) );
        start += time(0);
    }
    else{
        currentIteration = 0;
        acceptedModelPerturbation = 0;
        numberModelPerturbation = 0;
        acceptedTreePerturbation = 0;
        numberTreePerturbation = 0;
        start = time( 0 );
    }

    //initialise bestModelParams if beginning of the run or failure to restore the best state
    if (bestModelParams.size()==0){
        startingModel->getAllParameters(bestModelParams);
        bestLikelihood = startingTree->loglikelihood();
    }

    if (parameters.stringParameter( "Output format" ) == "phylip"){
        terminatorChar = ";";
    }
    else{
        if (parameters.stringParameter( "Output format" ) != "bambe"){
            cerr << "Output format must be \"phylip\" or \"bambe\"."
                 << endl;
            exit(EXIT_FAILURE);
        }
        terminatorChar = "";
    }

    startingTree->initSampling( parameters, overwrite );

    // to store the chain state regularly, open rescue for writing (in binary mode)
    filename = outputFileBaseName + string( ".rsc" );
    if ( rescueRead.is_open() ){
        // save rescueRead in memory and close it. Open as .rescue and copy back its content
        rescueRead.close();
    }
    rescue.open( filename.c_str(), ofstream::out | ofstream::binary | ofstream::ate );
    if ( !rescue.is_open() ){
        cerr << "error while creating the file " << filename << endl;
        exit(EXIT_FAILURE);
    }

    // to store the best state regularly, open rescueBest for writing (in binary mode)
    filename = outputFileBaseName + string( ".rsb" );
    // close rescueBestRead and open .rescue_best in append mode if recovering
    if ( rescueBestRead.is_open() ){
        rescueBestRead.close();
    }
    rescueBest.open( filename.c_str(), ofstream::out | ofstream::binary | ofstream::ate );
    if ( !rescueBest.is_open() ){
        cerr << "error while creating the file " << filename << endl;
        exit(EXIT_FAILURE);
    }

   /** ************************************************************************
    * Create/recover global output file
    ************************************************************************ */
    // Find out if the ".plot" file already exists and if it does
    // confirm overwrite
    filename = outputFileBaseName + string( ".plt" );
    FileParser::confirmOpenFile( plotfile, filename, overwrite  );
    //precision 4 for the outputs in plotfile
    plotfile.setf(ios::fixed);
    plotfile << setw(13) << setprecision(4);

    // Find out if the ".mp" file already exists and if it does
    // confirm overwrite (models are responsible for the precision and
    // the alignment of the numbers in this file)
    filename = outputFileBaseName + string( ".mp" );
    FileParser::confirmOpenFile( mpFile, filename, overwrite  );

    // Find out if the ".mp" file already exists and if it does
    // confirm overwrite (models are responsible for the precision and
    // the alignment of the figures in this file)
    filename = outputFileBaseName + string( ".hpm" );
    FileParser::confirmOpenFile( priorFile, filename, overwrite  );

    // Find out if the ".besttree" file already exists and if it does
    // confirm overwrite, the tree is responsible for the precision of
    // branch lengths
    filename = outputFileBaseName + string( "-best.tre" );
    FileParser::confirmOpenFile( bestTreeFile, filename, overwrite  );

    // Find out if the ".bestmodel" file already exists and if it does
    // confirm overwrite, models are responsible for the precision of the
    // figures
    filename = outputFileBaseName + string( "-best.mod" );
    FileParser::confirmOpenFile( bestModelParametersFile, filename, overwrite  );

    cout << "Chain launched" << endl;
    // launch the chain and wait....
    runMCMC( startingTree, startingModel, perturbationTree, perturbationModel,
            numberBurninIterations, numberSamplingIterations, samplingPeriod,
            currentIteration, terminatorChar );

    //delete the rescue files if the simulation finished
    filename = outputFileBaseName + string( ".rsc" );
    remove( filename.c_str() );
    filename = outputFileBaseName + string( ".rsb" );
    remove( filename.c_str() );

    // end the timer
    end = time( 0 );

    cout << "Time Taken = " << (int)( end - start ) << " seconds" << endl;
    output << "Time Taken = " << (int)( end - start ) << " seconds" << endl;
    cout << "Number of random calls = " << randBox.number_of_calls() << endl;

    // Check whether all parameters and categories have been used
    if (!parameters.checkAllUsed()){
        cerr << "(Double-check that these parameters were intended to be ignored.)" << endl;
    }

    delete startingTree;
    delete startingModel;
    delete table;
    return (EXIT_SUCCESS);
}



#include "Util/matrixmath.h"
#include <algorithm>
#include <numeric>

int main( int argc, char * argv[] ) {
    MCMCphase mcmcphase;

/*    void* modelModule;
    modelModule = dlopen("Gtr.dll", RTLD_NOW);
    if ( ! modelModule ){
        cerr << dlerror() << endl;
    exit(EXIT_FAILURE);
    }
    modelModule = dlopen("Rna7a.dll", RTLD_NOW);
    if ( ! modelModule ){
        cerr << dlerror() << endl;
    exit(EXIT_FAILURE);
    }
    modelModule = dlopen("Mixed.dll", RTLD_NOW);
    if ( ! modelModule ){
        cerr << dlerror() << endl;
    exit(EXIT_FAILURE);
    }
*/
    mcmcphase.run( argc, argv );
    return 0;
}
