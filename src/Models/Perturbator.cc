#include "Models/Perturbator.h"

#include "Util/ParametersSet.h"

#include "Models/PerturbatorParameter.h"

#include "Models/PerturbatorFrequenciesParameter.h"
#include "Models/PerturbatorGaussParameter.h"


#include <iostream>
#include <assert.h>

using namespace std;

Perturbator::Perturbator():
    rand(Singleton < randombox >::instance()){
    lastModified = -1;
    parameters.clear();
    rndTable.clear();
}

PerturbatorFrequenciesParameter* Perturbator::registerFrequencies( Observer<UpdateMessage>& obs, UpdateMessage::ParamFlags message,
    const string & name, vector<double> * frequencies,
    ParametersSet& parametersSet, const string & genericName, double defaultInitialFreqTuning,
    const string& initialTuningString ){
    PriorField priorField(parametersSet.stringParameter(genericName + ", prior")); 
    PerturbatorFrequenciesParameter* freqParam = new
            PerturbatorFrequenciesParameter( this, obs, message, name, frequencies,
                priorField, parametersSet, genericName, 
                defaultInitialFreqTuning, initialTuningString );
    this->parameters.insert(parameters.end(), freqParam);
    unsigned int freqPrio = 0 ;
    if (freqParam->getPriorType()!=PerturbatorFrequenciesParameter::FIXED){
        freqPrio = parametersSet.intParameter(genericName + ", proposal priority");
    }
    rndTable.insert( rndTable.end(), freqPrio, freqParam );
    return freqParam;
}

PerturbatorGaussParameter* Perturbator::registerGaussParameter( Observer<UpdateMessage>& obs, UpdateMessage::ParamFlags message,
        const string& name, double* element,
        const PriorField& prior, ParametersSet& parametersSet, const string & genericName, double initialStep, double minStep, double maxStep,
        const string& initialStepString ){
    PerturbatorGaussParameter* gaussParam = new
            PerturbatorGaussParameter( this, obs, message, name, element, prior, parametersSet, genericName,
            initialStep, minStep, maxStep, initialStepString );
    this->parameters.insert(parameters.end(), gaussParam);
    unsigned int ratesPrio = 0;
    if (gaussParam->getPriorType()!=PerturbatorGaussParameter::FIXED){
        ratesPrio = parametersSet.intParameter(genericName + ", proposal priority");
    }
    rndTable.insert( rndTable.end(), ratesPrio, gaussParam );
    return gaussParam;
}

PerturbatorGaussParameter* Perturbator::registerGaussParameter( Observer<UpdateMessage>& obs,  UpdateMessage::ParamFlags message,
        const string& name, double* element,
        ParametersSet& parametersSet, const string & genericName, double initialStep, double minStep, double maxStep,
        const string& initialStepString ){
    PriorField priorField(parametersSet.stringParameter(genericName + ", prior")); 
    PerturbatorGaussParameter* gaussParam = new
            PerturbatorGaussParameter( this, obs, message, name, element, priorField, parametersSet, genericName,
            initialStep, minStep, maxStep, initialStepString );
    this->parameters.insert(parameters.end(), gaussParam);
    unsigned int ratesPrio = 0;
    if (gaussParam->getPriorType()!=PerturbatorGaussParameter::FIXED){
        ratesPrio = parametersSet.intParameter(genericName + ", proposal priority");
    }
    rndTable.insert( rndTable.end(), ratesPrio, gaussParam );  
    return gaussParam;
}

double Perturbator::perturb(){
    if (rndTable.size()){
        assert(lastModified == -1);
        lastModified = (int)(rand.ran()*(double)rndTable.size());
        return rndTable[lastModified]->perturb();
    }
    return 0.0;
}

bool Perturbator::validatePerturbation(bool validation){
    if (rndTable.size()){
        assert(lastModified != -1);
        bool res = rndTable[lastModified]->validatePerturbation( validation );
        lastModified = -1;
        return res;
    }
    return true;
}

void Perturbator::stopBurn(){
    for ( vector<PerturbatorParameter*>::iterator iter = parameters.begin();
          iter != parameters.end(); ++iter ){
        (*iter)->stopBurn();
    }
}

void Perturbator::printPerturbationParameters(ostream& outputStream){
    for ( vector<PerturbatorParameter*>::iterator iter = parameters.begin();
          iter != parameters.end(); ++iter ){
        (*iter)->printPerturbationParameters( outputStream );
    }
}


void Perturbator::getAllPerturbationParameters( vector<double>& params ) const{
    params.clear();
    vector<double> p;
    for (unsigned int j=0; j < parameters.size(); ++j){
        parameters[j]->getAllPerturbationParameters( p );
        params.insert( params.end(), p.begin(), p.end() );
    }
}
    
void Perturbator::setAllPerturbationParameters( const vector<double>& params){
    assert( params.size() == getNumberPerturbationParameters() );
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp;
    for (unsigned int j=0; j < parameters.size(); ++j){
        nbp = parameters[j]->getNumberPerturbationParameters();
        p.resize(nbp);
        for ( unsigned int i = 0; i < nbp; ++i ){
            p[i] = *iter;
            ++iter;
        }
        parameters[j]->setAllPerturbationParameters( p );
    }
}
    
unsigned int Perturbator::getNumberPerturbationParameters() const{
    unsigned int numberPerturbationParameters = 0;
    for (unsigned int j=0; j < parameters.size(); ++j){
        numberPerturbationParameters += parameters[j]->getNumberPerturbationParameters();
    }
    return numberPerturbationParameters;
}

void Perturbator::getAllPriorParameters( vector<double>& params ) const{
    params.clear();
    vector<double> p;
    for (unsigned int j=0; j < parameters.size(); ++j){
        parameters[j]->getAllPriorParameters( p );
        params.insert( params.end(), p.begin(), p.end() );
    }
}
    
void Perturbator::setAllPriorParameters( const vector<double>& params){
    assert( params.size() == getNumberPriorParameters() );
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp;
    for (unsigned int j=0; j < parameters.size(); ++j){
        nbp = parameters[j]->getNumberPriorParameters();
        p.resize(nbp);
        for ( unsigned int i = 0; i < nbp; ++i ){
            p[i] = *iter;
            ++iter;
        }
        parameters[j]->setAllPriorParameters( p );
    }
}

unsigned int Perturbator::getNumberPriorParameters() const{
    unsigned int numberPriorParameters = 0;
    for (unsigned int j=0; j < parameters.size(); ++j){
        numberPriorParameters += parameters[j]->getNumberPriorParameters();
    }
    return numberPriorParameters;
}
    
double Perturbator::getLnPrior() const{
    double prior = 0.0;
    for (unsigned int j=0; j < parameters.size(); ++j){
        prior += parameters[j]->getLnPrior();
    }
    return prior;
}

