#include "Tree/MCMCTreeBasic.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include "Tree/TreeMap.h"
#include "Models/Perturbator.h"

#include <iostream>
#include <iomanip>

using namespace std;

MCMCTreeBasic::MCMCTreeBasic( const string & registrationName ):
InferenceTree(), randBox(Singleton<randombox>::instance()) {
    Singleton < Factory<MCMCTree> > & treeFactory = Singleton < Factory<MCMCTree> >::instance();
    treeFactory.subscribe( this, registrationName );
}

MCMCTreeBasic::MCMCTreeBasic( ParametersSet & parameters ):
InferenceTree(), randBox(Singleton<randombox>::instance()){
    if (false) parameters = parameters;
    hyperPerturbator = new Perturbator();
}

MCMCTreeBasic::~MCMCTreeBasic(){
    delete hyperPerturbator;
}


void MCMCTreeBasic::initialiseMCMC( ParametersSet& parameters ){
    lastPerturbation = INVALID;
    treePriority = (double)
            parameters.intParameter("Topology changes, proposal priority");
    branchPriority = (double)
            parameters.intParameter("Branch lengths, proposal priority");
    assert(branchPriority);
    if ( hyperPerturbator->getLoad() ){
        hyperPriority = (double)
                parameters.intParameter("Hyperpriors, proposal priority");
    }
    else{
        hyperPriority = 0.0;
    }
    double totPriority = (double)(branchPriority+treePriority+hyperPriority);
    branchPriority /= totPriority;
    hyperPriority /= totPriority;
    treePriority = 1.0 - branchPriority - hyperPriority;
}


double MCMCTreeBasic::perturb(){
    Singleton < randombox > & rndBox = Singleton < randombox >::instance();
    assert(lastPerturbation == INVALID);
    double res;
    double rndValue = rndBox.ran();
    if (rndValue < treePriority){
        lastPerturbation = TOPOLOGY;
        res = globalTopologyChange();
    }
    else{
        rndValue -= treePriority;
        if ( rndValue < branchPriority ) {
            lastPerturbation = LOCAL;
            res = changeBranchLength();
        }
        else{
            assert(hyperPriority);
            lastPerturbation = HYPER_PRIORS;
            res = hyperPerturbator->perturb();
        }
    }
    return res;
}


void MCMCTreeBasic::initSampling(ParametersSet& parameters, bool overwrite){
    //check whether initialiseMCMC has been called before
    assert(lastPerturbation == INVALID);
    string outputFileBaseName( parameters.stringParameter( "Output file" ) );
    string outputFormat( parameters.stringParameter( "Output format" ) );

    FileParser::confirmOpenFile( samplesFile,
             outputFileBaseName + ".smp", overwrite );
    transform( outputFormat.begin(), outputFormat.end(), outputFormat.begin(), (int(*)(int))tolower );
    if ( outputFormat == "phylip" ){
        terminator = string( ";" );
    }
    else {
        assert(outputFormat == "bambe");
        terminator = string( "" );
    }

    FileParser::confirmOpenFile( branchLengthsFile,
             outputFileBaseName + ".bl", overwrite );
    branchLengthsFile.setf(ios::fixed);
    branchLengthsFile << setprecision(4);

    //hyperprior file
    if ( hyperPerturbator->getLoad() ){
        FileParser::confirmOpenFile( hyperParametersFile,
                 outputFileBaseName + ".hpt", overwrite );
        hyperParametersFile.setf(ios::fixed);
        hyperParametersFile << setprecision(3);
    }
}

void MCMCTreeBasic::sample(){
    //topology sampling
    samplesFile << toStringNumbered( true ) << terminator << endl;
    vector < BasicNode * >::const_iterator iter = nodeRefVector.begin();
    ++iter;
    //branch lengths sampling
    while ( iter != nodeRefVector.end() ){
        branchLengthsFile << setw(9) << (*iter)->getParentDistance() << ' ';
        ++iter;
    }
    branchLengthsFile << endl;
    if ( hyperPerturbator->getLoad() ){
        vector<double> params;
        getAllPriorParameters( params );
        for ( vector < double >::const_iterator iter = params.begin();
              iter != params.end(); ++iter ){
            hyperParametersFile << setw(9) << *iter << ' ';
        }
        //MCMCTreeBasic is responsible for the line skip
        hyperParametersFile << endl;
    }
}

void MCMCTreeBasic::printPerturbationParameters( ostream& outputStream ){
    hyperPerturbator->printPerturbationParameters( outputStream );
}


void MCMCTreeBasic::stopBurn(){
    hyperPerturbator->stopBurn();
}

bool MCMCTreeBasic::validatePerturbation( bool validation ){
    assert(lastPerturbation != INVALID);
    bool res = false;
    switch (lastPerturbation){
        case TOPOLOGY: res = validateTopologyChange( validation ); break;
        case LOCAL: res = validateBranchLength( validation ); break;
        case HYPER_PRIORS: res=hyperPerturbator->validatePerturbation( validation ); break;
        default : assert("ERROR: in validation perturbation variable lastPerturbation unknown"==0);
    }
    lastPerturbation = INVALID;
    return res;
}


void MCMCTreeBasic::update( UpdateMessage* message ){
    assert(message->hasType(UpdateMessage::MODEL_TYPE));
    if (message->hasFlag(UpdateMessage::MODEL_PRIOR_FLAG)){
        return;
    }
    int cat = -1;
    if (message->hasFlag(UpdateMessage::SYMBOL_CATEGORY)){
        cat = (int)message->modelMsg.symbolCategory;
    }
#ifdef DEBUG1
    cout << "MCMC tree received  update from symbol category : " << cat << endl;
#endif
    invalidateNodes( cat );
}


void MCMCTreeBasic::retrieveNodesAux(){
    retrieveNodes( -1 );
}

void MCMCTreeBasic::saveNodesAux(){
    saveNodes( -1 );
}


void MCMCTreeBasic::getAllParameters( vector<double>& params ) const{
    getBranchVector( params );
}

void MCMCTreeBasic::setAllParameters( const vector<double>& params){
    setBranchVector( params );
}

unsigned int MCMCTreeBasic::getNumberTreeParameters() const{
    return getNumberBranches();
}


void MCMCTreeBasic::getAllPerturbationParameters( vector<double>& params ) const{
    hyperPerturbator->getAllPerturbationParameters(params);
}

void MCMCTreeBasic::setAllPerturbationParameters( const vector<double>& params){
    hyperPerturbator->setAllPerturbationParameters(params);
}

unsigned int MCMCTreeBasic::getNumberPerturbationParameters() const{
    return hyperPerturbator->getNumberPerturbationParameters();
}

void MCMCTreeBasic::getAllPriorParameters( vector<double>& params ) const{
    hyperPerturbator->getAllPriorParameters( params );
}

void MCMCTreeBasic::setAllPriorParameters( const vector<double>& params){
    hyperPerturbator->setAllPriorParameters( params );
}

unsigned int MCMCTreeBasic::getNumberPriorParameters() const{
   return hyperPerturbator->getNumberPriorParameters();
}

double MCMCTreeBasic::getLnPrior() const{
   return hyperPerturbator->getLnPrior();
}
