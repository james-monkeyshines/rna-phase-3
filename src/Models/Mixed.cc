#include "Models/Mixed.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"

#include "Sequence/SequenceTable.h"
#include "Models/MixedPerturbator.h"


#include <algorithm>


#define SQR(n) (n*n)

Mixed Mixed::prototype( "MIXED" );

Mixed::Mixed( const string & registrationName ) {
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model> >::instance();
    modelFactory.subscribe( this, registrationName );
}


Mixed::Mixed( ParametersSet & parameters ) {
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model> >::instance();

    Model* tempModel;
    string tempModelName;

    char categoryName[25];
    model.clear();
    int numberModels = parameters.intParameter("Number of models");
    for ( int i = 0; i < numberModels; ++i ){
        sprintf( categoryName, "MODEL%d", i+1 );
        tempModelName = parameters( categoryName ).stringParameter( "Model" );
        tempModel = modelFactory.create( tempModelName,
                                         parameters( categoryName ) );
        model.push_back(tempModel);
    }
    perturbator = NULL;
}

Mixed::~Mixed() {
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
ParametersSet Mixed::getModelParameters() const{
    ParametersSet parameters("Mixed model parameters");

    parameters["Model name"] = getName();
    char param[50];
    char averageSubstitutionRateRatio[20];
    char avRate[20];
    for ( unsigned int i = 1; i < model.size(); ++i ){
        sprintf( averageSubstitutionRateRatio, "%.8f",
                 model[i]->getAverageSubstitutionRate(0)/averageSubstitutionRate);
        sprintf( param, "Model%d/model1 average substitution rate ratio" , i+1 );
        parameters[param] = averageSubstitutionRateRatio;
    }
	sprintf(avRate, "%.8f", averageSubstitutionRate);
    parameters["Average substitution rate"] = avRate;

    char categoryName[15];
    for ( unsigned int i = 0; i < model.size(); ++i ){
        sprintf( categoryName, "MODEL%d", i+1);
        ParametersSet& newParametersSet = parameters("");
        newParametersSet = model[i]->getModelParameters();
        newParametersSet.setName( categoryName );
    }
    return parameters;
}


void Mixed::setModelParameters( ParametersSet& parameters ){
    if (parameters.stringParameter("Model name")!=getName()){
        cerr << "You are trying to initialise the model " << getName() << endl;
        cerr << "with parameters for the model "
             << parameters.stringParameter("Model name") << endl;
        cerr << "Exit..." << endl;
        exit(EXIT_FAILURE);
    }
    model[0]->setModelParameters(parameters("MODEL1"));
    model[0]->setAverageSubstitutionRate( averageSubstitutionRate, 0 );
    char param[60];
    for ( unsigned int i = 1; i < model.size(); ++i ){
        sprintf(param, "MODEL%d",i+1);
        model[i]->setModelParameters(parameters(param));
        sprintf(param, "Model%d/model1 average substitution rate ratio", i+1);
        model[i]->setAverageSubstitutionRate( parameters.doubleParameter(param)*averageSubstitutionRate, 0 );
    }
}


void Mixed::getAllParameters( vector < double > & parameters ) const {

    vector<double> partialParameters;

    parameters.clear();
    parameters.reserve( getNumberFreeParameters() );

    for ( unsigned int i = 1; i < model.size(); ++i ){
        parameters.push_back( sqrt( model[i]->getAverageSubstitutionRate(0) ) );
    }
    for ( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->getAllParameters( partialParameters );
        parameters.insert( parameters.end(), partialParameters.begin(),
                           partialParameters.end() );
    }
    assert(parameters.size()==getNumberFreeParameters());
}

void Mixed::setAllParameters( const vector < double > & newParameters ) {

    int current = 0;

    vector <double> partialParameters;

    for ( unsigned int i = 1; i < model.size(); ++i ){
        model[i]->setAverageSubstitutionRate( SQR( newParameters[current] ), 0 );
        ++current;
    }

    for( unsigned int i = 0; i < model.size(); ++i ){
        partialParameters.resize(model[i]->getNumberFreeParameters());
        for (unsigned int j = 0; j < partialParameters.size(); ++j){
            partialParameters[j] = newParameters[current];
            ++current;
        }
        model[i]->setAllParameters(partialParameters);
    }
    assert( current == (int)newParameters.size() );
}


void Mixed::setAverageSubstitutionRate( double newAverageRate,
                                     unsigned int symbolCategory ){
    //assert(symbolCategory==0);
    if (symbolCategory==0){
        averageSubstitutionRate = newAverageRate;
    }
    model[symbolCategory]->setAverageSubstitutionRate( newAverageRate, 0 );
}

void Mixed::validChange(){
    for ( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->validChange();
    }
}

unsigned int Mixed::getNumberFreeParameters() const{

    unsigned int ret = model.size() - 1;
    for ( unsigned int i = 0; i < model.size(); ++i ){
        ret += model[i]->getNumberFreeParameters();
    }
    return ret;
}

unsigned int Mixed::getNumberFreeFrequencyParameters() const{
    unsigned int ret = 0;
    for ( unsigned int i = 0; i < model.size(); ++i ){
        ret += model[i]->getNumberFreeFrequencyParameters();
    }
    return ret;
}

void Mixed::getOptimisableParameters( bool empiricalFreqs, vector<unsigned int> &optimisableParameters ) const{
    unsigned int nbModelParameters = getNumberFreeParameters();
	
	if ( empiricalFreqs ) {
		unsigned int nbFreqParameters = getNumberFreeFrequencyParameters();
		unsigned int nbRelativeRates = model.size() - 1;
		optimisableParameters.resize( nbModelParameters - nbFreqParameters );
		unsigned int cursor = 0;
		unsigned int k = 0;
		
		// All relative substitution rate parameters are together at the start of the vector
		while ( cursor < nbRelativeRates ){
			optimisableParameters.at(k++) = cursor++;
		}
		
		for ( unsigned int i = 0; i < model.size(); ++i ){
			unsigned int n = model[i]->getNumberFreeParameters();
			unsigned int m = model[i]->getNumberFreeFrequencyParameters();
			
			// Frequency parameters are the first ones within each model's set
        	for ( unsigned int j = cursor+m; j < cursor+n; ++j ){
				optimisableParameters.at(k++) = j;
			}
			cursor += n;
		}
	}
	else{
		optimisableParameters.resize(nbModelParameters);
		for ( unsigned int i = 0; i < nbModelParameters; ++i ){
			optimisableParameters.at(i) = i;
		}
	}
}

void Mixed::initialiseMCMC( ParametersSet & parameters ) {

    perturbator = new MixedPerturbator();
    char param[30];
    char label[30];

    sprintf( param, "PERTURBATION_MODEL" );
    char* numberPos = param;
    while(*numberPos) ++numberPos;
    if (!parameters.findParameter("Average rates, prior")){
        parameters["Average rates, prior"] = "uniform(0.005,200.0)";
    }
    for( unsigned int i = 0; i < model.size(); ++i ){
        sprintf( numberPos, "%d", i+1 );
        model[i]->initialiseMCMC( parameters( param ) );
        sprintf( label, "Model %d", i+1 );
        perturbator->registerModel( model[i], label,
                                    parameters, "Average rates",
                                    &model[i]->averageSubstitutionRate,
                                    .2, ", initial step" );
        model[i]->attach( *this );
    }
}

void Mixed::printPerturbationParameters( ostream & outputStream ) {
    perturbator->printPerturbationParameters( outputStream );
}

double Mixed::perturb() {
#ifdef DEBUG2
    cout << "Mixed model receives a perturb order transmit to perturbator" << endl;
#endif
    return perturbator->perturb();
}

void Mixed::stopBurn(){
    perturbator->stopBurn();
    for ( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->stopBurn();
    }
}

void Mixed::getAllPerturbationParameters( vector<double>& params ) const{
    perturbator->getAllPerturbationParameters( params );
}

void Mixed::setAllPerturbationParameters( const vector<double>& params){
    perturbator->setAllPerturbationParameters( params );
}

unsigned int Mixed::getNumberPerturbationParameters() const{
    return perturbator->getNumberPerturbationParameters();
}

void Mixed::getAllPriorParameters( vector<double>& params ) const{
    perturbator->getAllPriorParameters( params );
}

void Mixed::setAllPriorParameters( const vector<double>& params){
    perturbator->setAllPriorParameters( params );
}

unsigned int Mixed::getNumberPriorParameters() const{
    return perturbator->getNumberPriorParameters();
}

double Mixed::getLnPrior() const{
    return perturbator->getLnPrior();
}

void Mixed::initialiseML( ParametersSet & parameters ){
    char param[30];

    sprintf( param, "PENALTY_MODEL" );
    char* numberPos = param;
    while(*numberPos) ++numberPos;
    for( unsigned int i = 0; i < model.size(); ++i ){
        sprintf( numberPos, "%d", i+1 );
        model[i]->initialiseML( parameters( param ) );
    }
}

void Mixed::diffLnPenalty( vector<double>& gradVector ) const{
    gradVector.clear();
    for( unsigned int i = 0; i < model.size(); ++i ){
        vector<double> p;
        model[i]->diffLnPenalty(p);
        gradVector.insert(gradVector.end(),p.begin(),p.end());
    }
}

void Mixed::getAllPenaltyParameters( vector<double>& params ) const{
    params.resize(getNumberPenaltyParameters());
    vector<double>::iterator last = params.begin();
    for( unsigned int i = 0; i < model.size(); ++i ){
        vector<double> p;
        model[i]->getAllPenaltyParameters(p);
        last = copy(p.begin(),p.end(),last);
    }
    assert(last==params.end());
}

void Mixed::setAllPenaltyParameters( const vector<double>& params){
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

unsigned int Mixed::getNumberPenaltyParameters() const{
    unsigned int nb = 0;
    for( unsigned int i = 0; i < model.size(); ++i ){
        nb += model[i]->getNumberPenaltyParameters();
    }
    return nb;
}

double Mixed::getLnPenalty() const{
    double penalty = 0.0;
    for( unsigned int i = 0; i < model.size(); ++i ){
        penalty += model[i]->getLnPenalty();
    }
    return penalty;
}

bool Mixed::validatePerturbation( bool validation ) {
#ifdef DEBUG2
    cout << "Mixed model receives a validation " << validation << "transmit to perturbator" << endl;
#endif
    return perturbator->validatePerturbation( validation );
}

void Mixed::initialisation( SequenceTable * sequenceTable, int ) {
    if (sequenceTable){
        if ( sequenceTable->getNumberCategories() != model.size() ){
            cerr << "Error: the number of models declared in the MIXED model ("
                 << model.size() << ") does not match the number of categories ("
                 << sequenceTable->getNumberCategories() << ") in the datafile." << endl;
            exit(EXIT_FAILURE);
        }
    }
    for ( unsigned int i = 0; i < model.size(); ++i ){
        model[i]->initialisation( sequenceTable, i );
    }
    averageSubstitutionRate = 1.0;
}


void Mixed::printParameters( ostream & outputStream ) const {
    for (unsigned int i = 0; i < model.size(); ++i ){
        outputStream << "MODEL" << (i+1) << endl;
        model[i]->printParameters( outputStream );
    }
}



void Mixed::printLine( ostream & outputStream ) {
    outputStream.setf(ios::fixed);
    outputStream << setprecision( 4 );
    for (unsigned int i = 0; i < model.size(); ++i ){
        if ( i!= 0 ){
            outputStream << setw( 7 ) << model[i]->getAverageSubstitutionRate( 0 ) << ' ';
        }
        model[i]->printLine( outputStream );
        outputStream << "     ";
    }
    outputStream.unsetf(ios::fixed);
}

void Mixed::printPriorLine( ostream & outputStream ) {
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

void Mixed::fromLine( const vector<double>& parameters ){
    vector<double>::const_iterator iter = parameters.begin();
    for (unsigned int i = 0; i < model.size(); ++i ){
        if ( i!= 0 ){
            model[i]->setAverageSubstitutionRate( *iter, 0 ) ;
            ++iter;
        }
        model[i]->fromLine(
                vector<double>(iter,iter+model[i]->getNumberLineParameters()) );
       iter = iter + model[i]->getNumberLineParameters();
    }
}


unsigned int Mixed::getNumberLineParameters() const{
    unsigned int ret = model.size()-1;
    for (unsigned int i = 0; i < model.size(); ++i ){
        ret += model[i]->getNumberLineParameters();
    }
    return ret;
}

void Mixed::update(UpdateMessage* subject){
    clear();
    setType(MODEL_TYPE);
    setFlag(SYMBOL_CATEGORY);
    if (subject->hasType(PARAM_TYPE)){
        if (subject->hasFlag(PRIOR_FLAG)){
            setFlag(MODEL_PRIOR_FLAG);
#ifdef DEBUG2
            cout << "MIXED:received an update message from its average rate prior" << endl;
#endif
        }
        else{
#ifdef DEBUG2
            cout << "MIXED:received an update message from its average rate" << endl;
#endif
            //if the message does not come from any underlying model it must come from the average
            //substitution rate parameter of this mixed model.
            //update model[0]->averageSubstitutionRate with the new value of this->averageSubstitutionRate
            model[0]->setAverageSubstitutionRate( averageSubstitutionRate, 0 );
            model[0]->validChange();
        }
        modelMsg.symbolCategory = 0;
    }
    else{
        assert(subject->hasType(MODEL_TYPE));
        if (subject->hasFlag(MODEL_PRIOR_FLAG)){
            setFlag(MODEL_PRIOR_FLAG);
        }
        vector< Model* >::iterator iter = find( model.begin(), model.end(), subject );
        assert(iter!=model.end());
        unsigned int symbolCat = iter - model.begin();
#ifdef DEBUG2
        cout << "MIXED:received an update message, assigned to symbol " << symbolCat << endl;
        cout << "transmit" << endl;
#endif
        modelMsg.symbolCategory = symbolCat;
    }
    notify();
}
