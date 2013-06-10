#include "Models/Heterogeneous.h"

#include <algorithm>

#include "Models/MixedPerturbator.h"
#include "Models/Mixed.h"
#include "Models/MatrixModel.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/randombox.h"
#include "Util/ParametersSet.h"

#define SQR(n) (n*n)

using namespace std;

Heterogeneous Heterogeneous::prototype( "HETEROGENEOUS" );

Heterogeneous::Heterogeneous( const string & registrationName ) {
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model> >::instance();
    modelFactory.subscribe( this, registrationName );
}


Heterogeneous::Heterogeneous( ParametersSet & parameters ) {
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model> >::instance();

    Model* baseModel;
    string baseModelName;
    char str[50];

    model.clear();
    unsigned int numberModels = parameters.intParameter("Number of models");
    baseModelName = parameters("BASEMODEL").stringParameter("Model");
    for( unsigned int count = 0; count < numberModels; ++count ){
        baseModel = modelFactory.create( baseModelName,
                                     parameters("BASEMODEL") );
        model.push_back( baseModel );
    }
    sprintf( str, "ANCESTRAL_FREQUENCIES" );
    if (getNumberSymbolCategory()>1){
        char* numberPos = str;
        while(*numberPos) ++numberPos;
        ancestralFrequencies.resize(getNumberSymbolCategory());
        for( unsigned int count = 0; count < getNumberSymbolCategory(); ++count ){
            sprintf( numberPos, "%d", count+1 );
            ancestralFrequencies[count] =
                new AncestralFrequencies( parameters(str),
                            baseModel->getNumberGammaCategories(count),
                            baseModel->getInvariant(count) );
            parameters.touchCategory( str );
        }
    }
    else{
        ancestralFrequencies.resize(1);
        ancestralFrequencies[0] = new AncestralFrequencies( parameters(str),
                            baseModel->getNumberGammaCategories(0),
                            baseModel->getInvariant(0) );
        parameters.touchCategory( str );
    }
    perturbator = NULL;
}

Heterogeneous::~Heterogeneous() {
    //delete all the models
    for ( vector<Model*>::iterator iter = model.begin();
          iter != model.end(); ++iter ) {
        delete (*iter);
    }
}


/** ************************************************************************
 * getModelParameters
 * @semantics  create a ParametersSet and fill it with the parameters of
 *             the model
 ************************************************************************ */
ParametersSet Heterogeneous::getModelParameters() const{

    assert( averageSubstitutionRate == 1.0 );

    ParametersSet parameters("Heterogeneous model parameters");

    parameters["Model name"] = getName();
    char categoryName[30];
    char paramName[50];
    char paramValue[20];
    for ( unsigned int i = 0; i < getNumberSymbolCategory(); ++i ){
        sprintf( categoryName, "ANCESTRALFREQUENCIES%d" , i+1 );
        ancestralFrequencies[i]->getModelParameters( parameters(categoryName), this, i );
    }
    for ( unsigned int i = 1; i < model.size(); ++i ){
        sprintf( paramValue, "%.8f",
                 model[i]->getAverageSubstitutionRate(0)/averageSubstitutionRate );
        sprintf( paramName, "Model%d/model1 average substitution rate ratio" , i+1 );
        parameters[paramName] = paramValue;
    }

    for ( unsigned int i = 0; i < model.size(); ++i ){
        sprintf( categoryName, "MODEL%d", i+1);
        //we cannot give a name to newParametersSet yet
        ParametersSet& newParametersSet = parameters("");
        //because this copy will erase it at once
        newParametersSet = model[i]->getModelParameters();
        //set the name once newParametersSet is ready
        newParametersSet.setName( categoryName );
    }
    return parameters;
}


void Heterogeneous::setModelParameters( ParametersSet& parameters ){
    if (parameters.stringParameter("Model name")!=getName()){
        cerr << "You are trying to initialise the model " << getName() << endl;
        cerr << "with parameters for the model "
             << parameters.stringParameter("Model name") << endl;
        cerr << "Exit..." << endl;
        exit(EXIT_FAILURE);
    }
    char categoryName[60];
    char paramName[60];
    for ( unsigned int i = 0; i < getNumberSymbolCategory(); ++i ){
        sprintf( categoryName, "ANCESTRALFREQUENCIES%d" , i+1 );
        ancestralFrequencies[i]->setModelParameters( parameters(categoryName), this, i );
    }
    model[0]->setModelParameters(parameters("MODEL1"));
    model[0]->setAverageSubstitutionRate( averageSubstitutionRate, 0);
    for ( unsigned int i = 1; i < model.size(); ++i ){
        sprintf(categoryName, "MODEL%d",i+1);
        model[i]->setModelParameters( parameters(categoryName) );
        sprintf(paramName, "Model%d/model1 average substitution rate ratio", i+1);
        model[i]->setAverageSubstitutionRate(
                parameters.doubleParameter(paramName) * averageSubstitutionRate, 0 );
    }
}


void Heterogeneous::getAllParameters( vector < double > & parameters ) const {

    vector<double> partialParameters;

    parameters.clear();
    parameters.reserve( getNumberFreeParameters() );

    for (unsigned int i = 0; i < getNumberSymbolCategory(); ++i){
        ancestralFrequencies[i]->getAllParameters( partialParameters );
        parameters.insert( parameters.end(), partialParameters.begin(),
                           partialParameters.end() );
    }

    for ( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->getAllParameters( partialParameters );
        parameters.insert( parameters.end(), partialParameters.begin(),
                           partialParameters.end() );
    }
    for ( unsigned int i = 1; i < model.size(); ++i ){
        parameters.push_back( sqrt( model[i]->getAverageSubstitutionRate(0) ) );
    }
}

void Heterogeneous::setAllParameters( const vector < double > & newParameters ) {

    assert( newParameters.size() == getNumberFreeParameters() );

    vector < double >::const_iterator iter1;
    vector < double >::const_iterator iter2;

    iter1 = newParameters.begin();
    for (unsigned int i = 0; i < getNumberSymbolCategory(); ++i){
        iter2 = iter1 + ancestralFrequencies[i]->getNumberFreeParameters();
        ancestralFrequencies[i]->setAllParameters( vector<double>( iter1, iter2 ) );
        iter1=iter2;
    }
    for( unsigned int i = 0; i < model.size(); ++i ){
        iter2 = iter1 + model[i]->getNumberFreeParameters();
        model[i]->setAllParameters( vector<double>( iter1, iter2 ) );
        iter1=iter2;
    }
    for ( unsigned int i = 1; i < model.size(); ++i ){
        model[i]->setAverageSubstitutionRate(SQR((*iter1)),0);
        ++iter1;
    }
    assert( iter1 == newParameters.end() );
}

void Heterogeneous::validChange(){
    for ( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->validChange();
    }
}

unsigned int Heterogeneous::getNumberFreeParameters() const{
    unsigned int ret = model[0]->getNumberFreeParameters() * model.size();
    ret += model.size() - 1;
    for ( unsigned int i = 0; i < getNumberSymbolCategory(); ++i ){
        ret += ancestralFrequencies[i]->getNumberFreeParameters();
    }
    return ret;
}

unsigned int Heterogeneous::getNumberFreeFrequencyParameters() const{
	return 0;
}

void Heterogeneous::getOptimisableParameters( bool empiricalFreqs, vector<unsigned int> &optimisableParameters ) const{
	unsigned int nbModelParameters = getNumberFreeParameters();
	
	optimisableParameters.resize(nbModelParameters);
	for ( unsigned int i = 0; i < nbModelParameters; ++i ){
		optimisableParameters.at(i) = i;
	}
}

void Heterogeneous::initialiseMCMC( ParametersSet & parameters ) {

    perturbator = new MixedPerturbator();
    char param[30];
    char label[50];
    char* numberPos;
    sprintf( param, "PERTURBATION_BASEMODEL" );
    if (!parameters.findParameter("Average rates, prior")){
        parameters["Average rates, prior"] = "uniform(0.005,200.0)";
    }
    if (!parameters.findParameter("Average rates, proposal priority")){
        parameters["Average rates, proposal priority"] = "0";
    }
    //register the models first (this is compulsory)
    for( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->initialiseMCMC( parameters( param ) );
        perturbator->registerModel( model[i], "Model", parameters, "Average rates",
                                    &model[i]->averageSubstitutionRate,
                                    .2, ", initial step" );

        model[i]->attach( *this );
    }

    sprintf(label,"ANCESTRAL_FREQUENCIES");
    if (ancestralFrequencies.size()==1){
        ancestralFrequencies[0]->initialiseMCMC( parameters(label), perturbator );
        ancestralFrequencies[0]->attach( *this );
    }
    else{
        numberPos = label;
        while(*numberPos) ++numberPos;
        for( unsigned int i = 0; i < ancestralFrequencies.size(); ++i ){
            sprintf(numberPos,"%d", i+1);
            ancestralFrequencies[i]->initialiseMCMC( parameters(label), perturbator );
            ancestralFrequencies[i]->attach( *this );
        }
    }
}

void Heterogeneous::printPerturbationParameters( ostream & outputStream ) {
    perturbator->printPerturbationParameters( outputStream );
}

double Heterogeneous::perturb() {
    return perturbator->perturb();
}

void Heterogeneous::stopBurn(){
    perturbator->stopBurn();
    //there is no stopburn in MixedPerturbator, we have to transmit
    //the call ourself
    for ( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->stopBurn();
    }
}

void Heterogeneous::getAllPerturbationParameters( vector<double>& params ) const{
    perturbator->getAllPerturbationParameters( params );
}

void Heterogeneous::setAllPerturbationParameters( const vector<double>& params){
    perturbator->setAllPerturbationParameters( params );
}

unsigned int Heterogeneous::getNumberPerturbationParameters() const{
    return perturbator->getNumberPerturbationParameters();
}


void Heterogeneous::getAllPriorParameters( vector<double>& params ) const{
    perturbator->getAllPriorParameters( params );
}

void Heterogeneous::setAllPriorParameters( const vector<double>& params){
    perturbator->setAllPriorParameters( params );
}

unsigned int Heterogeneous::getNumberPriorParameters() const{
    return perturbator->getNumberPriorParameters();
}

double Heterogeneous::getLnPrior() const{
    return perturbator->getLnPrior();
}


void Heterogeneous::initialiseML( ParametersSet & parameters ){
    char param[70];

    cerr << "Time-heterogeneous models should not be used with ML methods" << endl;
    exit(EXIT_FAILURE);

    sprintf( param, "PENALTY_BASEMODEL" );
    for( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->initialiseML( parameters( param ) );
    }
    sprintf( param, "PENALTY_ANCESTRALFREQUENCIES" );
    if (ancestralFrequencies.size()==1){
        ancestralFrequencies[0]->initialiseML( parameters(param) );
    }
    else{
        char* numberPos = param;
        while(*numberPos) ++numberPos;
        for( unsigned int i = 0; i < ancestralFrequencies.size(); ++i ){
            sprintf(numberPos,"%d", i+1);
            ancestralFrequencies[i]->initialiseML( parameters(param) );
        }
    }
}

void Heterogeneous::diffLnPenalty( vector<double>& gradVector ) const{

    cerr << "Time-heterogeneous models should not be used with ML methods" << endl;
    exit(EXIT_FAILURE);

    gradVector.clear();
    for( unsigned int i = 0; i < model.size(); ++i ){
        vector<double> p;
        model[i]->diffLnPenalty(p);
        gradVector.insert(gradVector.end(),p.begin(),p.end());
    }
}

void Heterogeneous::getAllPenaltyParameters( vector<double>& params ) const{

    cerr << "Time-heterogeneous models should not be used with ML methods" << endl;
    exit(EXIT_FAILURE);

    params.resize(getNumberPenaltyParameters());
    vector<double>::iterator last = params.begin();
    for( unsigned int i = 0; i < model.size(); ++i ){
        vector<double> p;
        model[i]->getAllPenaltyParameters(p);
        last = copy(p.begin(),p.end(),last);
    }
    assert(last==params.end());
}

void Heterogeneous::setAllPenaltyParameters( const vector<double>& params){

    cerr << "Time-heterogeneous models should not be used with ML methods" << endl;
    exit(EXIT_FAILURE);

    assert(params.size() == getNumberPenaltyParameters());
    vector<double>::const_iterator iter1 = params.begin();
    vector<double>::const_iterator iter2;
    for( unsigned int i = 0; i < model.size(); ++i ){
        iter2=iter1+model[i]->getNumberPenaltyParameters();
        vector<double> p(iter1,iter2);
        model[i]->setAllPenaltyParameters(p);
        iter1=iter2;
    }
    assert(iter1==params.end());
}

unsigned int Heterogeneous::getNumberPenaltyParameters() const{

    cerr << "Time-heterogeneous models should not be used with ML methods" << endl;
    exit(EXIT_FAILURE);

    unsigned int nb = 0;
    for( unsigned int i = 0; i < model.size(); ++i ){
        nb += model[i]->getNumberPenaltyParameters();
    }
    return nb;
}

double Heterogeneous::getLnPenalty() const{

    cerr << "Time-heterogeneous models should not be used with ML methods" << endl;
    exit(EXIT_FAILURE);

    double penalty = 0.0;
    for( unsigned int i = 0; i < model.size(); ++i ){
        penalty += model[i]->getLnPenalty();
    }
    return penalty;
}


bool Heterogeneous::validatePerturbation( bool validation ) {
    return perturbator->validatePerturbation( validation );
}

void Heterogeneous::initialisation( SequenceTable * sequenceTable, int modelId ) {
    assert( modelId == 0);
    for ( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->initialisation( sequenceTable, modelId );
    }
    vector<double> freq;
    for (unsigned int i = 0; i < getNumberSymbolCategory(); ++i){
        //OK - technically we should use model[0]->numberFrequencies(i)
        //but this could cause huge design issue at the moment
        //(HETEROGENEOUS does not know how to find frequency of a state
        //from a frequency vector and use a 1<->1 matching
        ancestralFrequencies[i]->initialisation( model[0]->getNumberStates(i) );
    }
    averageSubstitutionRate = 1.0;
}


void Heterogeneous::printParameters( ostream & outputStream ) const {
    for (unsigned int i = 0; i < getNumberSymbolCategory(); ++i){
        outputStream << "ANCESTRAL FREQUENCIES " << (i+1) << ":   ";
        ancestralFrequencies[i]->printParameters(outputStream);
        outputStream << endl;
        outputStream.setf(ios::fixed);
        outputStream << setprecision( 5 );
        for ( unsigned int freqCat = 0; freqCat < ancestralFrequencies[i]->getNumberFrequenciesSets(); ++freqCat ){
            if (ancestralFrequencies[i]->getNumberFrequenciesSets()>=2){
                outputStream << "set of frequencies " << freqCat + 1 << ":" << endl;
            }
            else{
                outputStream << "frequencies:" << endl;
            }
            for ( unsigned int j = 0; j < getNumberStates(i); ++j ){
                outputStream << "f[" << getState( j, i ) << "] = " <<  getFrequency( j, freqCat, i ) << "  ";
            }
            outputStream << endl;
        }
        outputStream << endl;
    }
    outputStream << endl;
    for (unsigned int i = 0; i < model.size(); ++i ){
        outputStream << "MODEL" << (i+1) << endl;
        model[i]->printParameters( outputStream );
    }
}


void Heterogeneous::printLine( ostream & outputStream ) {
    outputStream.setf(ios::fixed);
    outputStream << setprecision( 4 );
    for (unsigned int i = 0; i < getNumberSymbolCategory(); ++i){
        ancestralFrequencies[i]->printLine(outputStream);
        outputStream << "  ";
    }
    for (unsigned int i = 0; i < model.size(); ++i ){
        if ( i!= 0 ){
            outputStream << setw( 9 ) << model[i]->getAverageSubstitutionRate( 0 ) << ' ';
        }
        model[i]->printLine( outputStream );
        outputStream << "    ";
    }
    outputStream.unsetf(ios::fixed);
}

void Heterogeneous::printPriorLine( ostream & outputStream ) {
    vector< double > priors;
    //prior for underlying models ar caught with that call
    getAllPriorParameters( priors );
    outputStream.setf(ios::fixed);
    outputStream << setprecision( 4 );
    for (unsigned int i = 0; i < priors.size(); ++i){
        outputStream << setw( 14 ) << priors[i] << ' ';
    }
    outputStream.unsetf(ios::fixed);
}


void Heterogeneous::fromLine( const vector<double>& parameters ){
    vector<double>::const_iterator iter1 = parameters.begin();
    vector<double>::const_iterator iter2;
    for (unsigned int i = 0; i < getNumberSymbolCategory(); ++i){
        iter2 = iter1 + ancestralFrequencies[i]->getNumberLineParameters();
        ancestralFrequencies[i]->fromLine( vector<double>(iter1,iter2));
        iter1=iter2;
    }
    for (unsigned int i = 0; i < model.size(); ++i ){
        if ( i!= 0 ){
            model[i]->setAverageSubstitutionRate( (*iter1), 0 );
            ++iter1;
        }
        iter2 = iter1 + model[i]->getNumberLineParameters();
        model[i]->fromLine( vector<double>(iter1,iter2) );
        iter1 = iter2;
    }
}


unsigned int Heterogeneous::getNumberLineParameters() const{
    //ratio between heterogeneous models
    unsigned int ret = model.size()-1;
    //ancestral frequencies
    for (unsigned int i = 0; i < getNumberSymbolCategory(); ++i){
        ret += ancestralFrequencies[i]->getNumberLineParameters();
    }
    //parameters
    ret += model[0]->getNumberLineParameters()*model.size();
    return ret;
}

void Heterogeneous::update(UpdateMessage* subject){
    clear();
    setType(MODEL_TYPE);
    //message from an internal model
    if (subject->hasType(MODEL_TYPE)){
        if (subject->hasFlag(MODEL_PRIOR_FLAG)){
            setFlag(MODEL_PRIOR_FLAG);
        }
        //if internal models are mixed then preserve the symbol category information
        if (subject->hasFlag(SYMBOL_CATEGORY)){
            setFlag(SYMBOL_CATEGORY);
            modelMsg.symbolCategory = subject->modelMsg.symbolCategory;
        }
        vector< Model* >::iterator iterModel = find( model.begin(), model.end(), (Model*)subject );
        assert(iterModel!=model.end());
        unsigned int modelId = iterModel - model.begin();
        setFlag(MODEL_FLAG);
        modelMsg.model = modelId;
    }
    //message from ancestral frequency
    else{
        assert(subject->hasType(PARAM_TYPE));
        assert(subject->hasFlag(FREQ));
        vector< AncestralFrequencies* >::iterator iterFreq = find( ancestralFrequencies.begin(), ancestralFrequencies.end(), (Frequencies*)subject );
        assert( iterFreq != ancestralFrequencies.end() );
        unsigned int symbolCat = iterFreq + 1 - ancestralFrequencies.begin();
        setFlag(SYMBOL_CATEGORY);
        modelMsg.symbolCategory = symbolCat;
        //assign an inexistant model
        setFlag(MODEL_FLAG);
        modelMsg.model = model.size();
    }
    notify();
}
