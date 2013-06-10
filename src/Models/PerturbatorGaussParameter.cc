#include "Models/PerturbatorGaussParameter.h"

#include "Util/ParametersSet.h"

#include <iostream>
#include <assert.h>


PerturbatorGaussParameter::PerturbatorGaussParameter( Perturbator* pertubator, Observer<UpdateMessage>& obs,
                      UpdateMessage::ParamFlags message, const string& name, double* pvalue,
                      const PriorField& prior, ParametersSet& parameters, const string& genericName,
                      double initialStep, double minStep, double maxStep,
                      const string& initialStepString ):
PerturbatorParameter(obs, message),
PerturbatorHelper( pertubator, obs, message, name, prior, parameters, genericName, initialStep, minStep, maxStep, initialStepString){
    this->pvalue = pvalue;
    if ( ((*pvalue)>extrMax) || ((*pvalue)<extrMin) ){
        *pvalue = mean;
    }
    //initialise lnPrior;
    lnPrior = computeLnPrior();
}
        
double PerturbatorGaussParameter::perturb(){
    //save oldPrior
    oldLnPrior = lnPrior;
    double lnHR = perturbValue( pvalue );
    PerturbatorParameter::perturb();
    //compute newPrior
    lnPrior = computeLnPrior();
    return lnHR;
}
        
bool PerturbatorGaussParameter::validatePerturbation( bool validation ){
    bool res = validatePerturbationValue( validation, pvalue );
    //restore prior if perturbation was refused
    if (!validation){
        lnPrior = oldLnPrior;
    }
    PerturbatorParameter::validatePerturbation( validation );
    return res;
}

void PerturbatorGaussParameter::update(UpdateMessage* subject){
    PerturbatorHelper::update(subject);
    if ( !subject->hasFlag(UpdateMessage::PRIOR_FLAG) ){
        lnPrior = computeLnPrior();
        //warn for a change in the prior
        setFlag(PRIOR_FLAG);
        notify();
        unsetFlag(PRIOR_FLAG);
    }
}

void PerturbatorGaussParameter::stopBurn(){
    PerturbatorHelper::stopBurn();
}


void PerturbatorGaussParameter::printPerturbationParameters( ostream& outputStream ) const{
    PerturbatorHelper::printPerturbationParameters(outputStream);
}

void PerturbatorGaussParameter::getAllPerturbationParameters( vector<double>& params ) const{
    PerturbatorHelper::getAllPerturbationParameters(params);
}
    
void PerturbatorGaussParameter::setAllPerturbationParameters( const vector<double>& params){
    PerturbatorHelper::setAllPerturbationParameters(params);
}
    
unsigned int PerturbatorGaussParameter::getNumberPerturbationParameters() const{
    return PerturbatorHelper::getNumberPerturbationParameters();
}

void PerturbatorGaussParameter::getAllPriorParameters( vector<double>& params ) const{    
    PerturbatorHelper::getAllPriorParameters(params);
}

void PerturbatorGaussParameter::setAllPriorParameters( const vector<double>& params){
    PerturbatorHelper::setAllPriorParameters(params);
    lnPrior=computeLnPrior();
}
    
unsigned int PerturbatorGaussParameter::getNumberPriorParameters() const{
   return  PerturbatorHelper::getNumberPriorParameters();
}

double PerturbatorGaussParameter::computeLnPrior() const{
    return PerturbatorHelper::getLnPriorValue(*pvalue);
}


double PerturbatorGaussParameter::getLnPrior() const{
//    double mom = computeLnPrior();
//    if (lnPrior != mom){
//        cout << "old ln prior = " << lnPrior << endl;
//        cout << "should be " << mom << endl;
//        exit(EXIT_FAILURE);
//    }
#ifdef DEBUG1
    double realPrior = computeLnPrior();
#ifdef DEBUG4    
    cout << "old ln prior = " << lnPrior << endl;
    cout << "should be " << realPrior << endl;
#endif
    assert(lnPrior==realPrior);
#endif
    return lnPrior;
}
