#include "Models/PerturbatorHelper.h"

#include "Models/PriorField.h"
#include "Models/PerturbatorGaussParameter.h"

#include <assert.h>
#include <string>
#include <iostream>
#include <float.h>
#include <cstdlib>

PerturbatorHelper::PerturbatorHelper( Perturbator* perturbator, Observer<UpdateMessage>& obs, 
                      UpdateMessage::ParamFlags message, const string& name, const PriorField& prior,
                      ParametersSet& parameters, const string& genericName,
                      double initialStep, double minStep, double maxStep,
                      const string& initialStepString ):
PerturbatorBase(name, parameters, genericName, initialStep, minStep, maxStep, initialStepString){
    priorType = NB_TYPE;
    savedLnValue = 0.0;
    double bound1, bound2;
    double mean1;
    lowBoundaryShockCounter = 0;
    highBoundaryShockCounter = 0;
    lowBoundaryShock = false;
    highBoundaryShock = false;
    
    hyperPriors.clear();
    
    if (prior.name() == "fixed"){
        if ((prior.getNumberParameters()!=1)||!prior.isConstant(0)){
            cerr << "invalid fixed prior: " << prior.toString() << ". 1 constant value expected."
                 << endl;
            exit(EXIT_FAILURE);
        }
        extrMin = prior.getValue(0);
        extrMax = prior.getValue(0);
        mean = prior.getValue(0);
        priorType = FIXED;
    }
    if (prior.name() == "uniform"){
        if ( (prior.getNumberParameters()!=2) || !prior.isConstant(0) || !prior.isConstant(1) ){
            cerr << "invalid uniform prior: " << prior.toString() << ". 2 fixed boundaries expected."
                 << endl;
            exit(EXIT_FAILURE);
        }
        priorParams.resize(2);
        //cheap trick to prevent non proper uniform prior
        extrMin = prior.getValue(0);
        extrMax = prior.getValue(1);
        if (extrMin>=extrMax){
            cerr << "invalid uniform prior (min>=max): " << prior.toString() << endl;
            exit(EXIT_FAILURE);
        }
        mean = (extrMin + extrMax) / 2.0;
        priorType = UNIFORM;
        savedLnValue = -log(extrMax-extrMin);
    }
    if (prior.name() == "normal"){
        if (prior.getNumberParameters()!=2){
            cerr << "invalid gaussian prior: " << prior.toString() << ". Mean and standard deviation expected."
                 << endl;
            exit(EXIT_FAILURE);
        }
        priorParams.resize(2);
        bound1 = -DBL_MAX;
        bound2 = DBL_MAX;
        create( perturbator, obs, message, name+" mean hyperparameter", prior, 0, bound1, bound2, mean,
                      parameters, genericName+" mean hyperparameter" );
        extrMin = -DBL_MAX;
        extrMax = DBL_MAX;
        bound1 = 0.0;
        bound2 = DBL_MAX;
        create( perturbator, obs, message, name+" standard deviation hyperparameter", prior, 1, bound1, bound2, mean1,
                      parameters, genericName+" standard deviation hyperparameter");
        priorType = NORMAL;
    }
    if (prior.name() == "exponential"){
        if (prior.getNumberParameters()!=1){
            cerr << "invalid exponential prior: " << prior.toString() << ". Mean expected."
                 << endl;
            exit(EXIT_FAILURE);
        }
        priorParams.resize(1);
        bound1 = 0.0;
        bound2 = DBL_MAX;
        create( perturbator, obs, message,name+" exponential hyperparameter", prior, 0, bound1, bound2, mean1,
                      parameters, genericName+" exponential hyperparameter" );
        extrMin = 0.0;
        extrMax = DBL_MAX;
        if(mean1==0.0){
            cerr << "invalid exponential hyperparameter "
                 << prior.getHyperPriorField(0).toString() << endl;
            exit(EXIT_FAILURE);
        }
        else{
            //if the "rate of change" parameter is equal to mean1 then the mean
            //of the exp distribution is 1/mean1
            mean = 1/mean1;
        }
        priorType = EXP;
    }
    if (prior.name() == "lognormal"){
        if ( (!prior.getNumberParameters()) || (prior.getNumberParameters()>3) ){
            cerr << "invalid logexponential prior: " << prior.toString()
                 << ". Shape, scale(default=1.0) and location(default=0.0) expected."
                 << endl;
            exit(EXIT_FAILURE);
        }
        priorParams.resize(3);
        bound1 = 0.0;
        bound2 = DBL_MAX;
        create( perturbator, obs, message, name+" shape hyperparameter", prior, 0, bound1, bound2, mean1,
                      parameters, genericName+" shape hyperparameter" );
        mean = exp(0.5*mean1*mean1);
        if (prior.getNumberParameters()<2){
            priorParams[1] = 1.0;
        }
        else{
            bound1 = 0.0;
            bound2 = DBL_MAX;
            create( perturbator, obs, message, name+" scale hyperparameter", prior, 1, bound1, bound2, mean1,
                          parameters, genericName+" scale hyperparameter" );
            mean = mean * mean1;
        }
        if (prior.getNumberParameters()<3){
            priorParams[2] = 0.0;
            bound1 = 0.0;
        }
        else{
            bound1 = -DBL_MAX;
            bound2 = DBL_MAX;
            create( perturbator, obs, message, name+" location hyperparameter", prior, 2, bound1, bound2, mean1,
                          parameters, genericName+" location hyperparameter" );
        }
        extrMin = bound1;
        mean += extrMin;
        extrMax = DBL_MAX;
        priorType = LOGNORMAL;
    }
    if (prior.name() == "gamma"){
        if ( (!prior.getNumberParameters()) || (prior.getNumberParameters()>3) ){
            cerr << "invalid logexponential prior: " << prior.toString()
                 << ". Shape, scale(default=1.0) and location(default=0.0) expected."
                 << endl;
            exit(EXIT_FAILURE);
        }
        priorParams.resize(3);
        bound1 = 0.0;
        bound2 = DBL_MAX;
        create( perturbator, obs, message, name+" shape hyperparameter", prior, 0, bound1, bound2, mean1,
                      parameters, genericName+" shape hyperparameter" );
        mean = mean1;
        if (prior.getNumberParameters()<2){
            priorParams[1] = 1.0;
        }
        else{
            bound1 = 0.0;
            bound2 = DBL_MAX;
            create( perturbator, obs, message, name+" scale hyperparameter", prior, 1, bound1, bound2, mean1,
                          parameters, genericName+" scale hyperparameter" );
            mean = mean * mean1;
        }
        if (prior.getNumberParameters()<3){
            priorParams[2] = 0.0;
            bound1 = 0.0;
        }
        else{
            bound1 = -DBL_MAX;
            bound2 = DBL_MAX;
            create( perturbator, obs, message, name+" location hyperparameter", prior, 2, bound1, bound2, mean1,
                          parameters, genericName+" location hyperparameter" );
        }
        extrMin = bound1;
        mean += extrMin;
        extrMax = DBL_MAX;
        priorType = GAMMA;
    }
    if ( priorType == NB_TYPE ){
        cerr << "unknown prior distibution: " << prior.toString() << endl;
        exit(EXIT_FAILURE);
    }
    if ( (priorType != FIXED) && (extrMin>=extrMax) ){
        cerr << "invalid boundaries, prior: " << prior.toString() << endl;
        exit(EXIT_FAILURE);
    }
    //initialise savedLnValue for EXP, NORMAL & LOGNORMAL 
    //this value is updated later if an update message from a hyperparameter
    //is received
    updateSavedLnValue();
}


void PerturbatorHelper::create( Perturbator* perturbator, Observer<UpdateMessage>& ,
                      UpdateMessage::ParamFlags message, const string& name, const PriorField& prior,
                      unsigned int priorField, double& min, double& max, double& mean,
                      ParametersSet& parameters, const string& genericName){
    if (prior.isConstant(priorField)){
        if ( ( prior.getValue(priorField) >= min ) &&
             ( prior.getValue(priorField) <= max ) ){
            min = prior.getValue(priorField);
            max = min;
            mean = min;
        }
        else{
            cerr << name << "invalid: " << prior.toString() << endl;
            exit(EXIT_FAILURE);
        }
        priorParams[priorField]=mean;
    }
    else{
        hyperPriors.push_back(priorField);
        PerturbatorGaussParameter* pert =
            perturbator->registerGaussParameter( *this, message,
                name, &(priorParams[priorField]),
                prior.getHyperPriorField(priorField), parameters, genericName,
                .3, .0001, 10.0, "initial step" );
        hyperPriorsPert.push_back(pert);
        pert->getExtrema( min, max );
        mean = pert->getMean();
        priorParams[priorField]=mean;
        //this call is vital to refresh the prior computed for the old value
        pert->invalidate();
    }
}

double PerturbatorHelper::perturbValue( double* pvalue ){
    
#ifdef DEBUG4  
    cout << "perturbator gauss received a perturb order" << endl;
#endif    
    //modify the value, gaussian perturbation with reflection
    //at the boundaries
    assert (pvalue);
    oldValue = *pvalue;
    assert(extrMin<extrMax);
    *pvalue += rand.gasdev() * step;
    //stay within boundaries
    double delta = extrMax - extrMin;
    //return in range [extrMin-delta; extrMax+delta]
    if ( *pvalue <= (extrMin-delta) ){
        int n = (int)( (extrMax-*pvalue)/(2.0*delta) );
        *pvalue += 2.0 * delta * (double)n;
        lowBoundaryShock = true;
        highBoundaryShock = true;
    }
    else{
        if ( *pvalue >= (extrMax+delta) ){
            int n = (int)( (*pvalue-extrMin)/(2.0*delta) );
            *pvalue -= 2.0 * delta * (double)n;
            lowBoundaryShock = true;
            highBoundaryShock = true;
        }
    }
    assert( *pvalue > (extrMin-delta) );
    assert( *pvalue < (extrMax+delta) );
    //reflection
    if (*pvalue<extrMin){
        *pvalue = 2.0*extrMin-*pvalue;
        lowBoundaryShock = true;
    }
    if (*pvalue>extrMax){
        *pvalue = 2.0*extrMax-*pvalue;
        highBoundaryShock = true;        
    }
    if (lowBoundaryShock) ++lowBoundaryShockCounter;
    if (highBoundaryShock) ++highBoundaryShockCounter;
    assert( *pvalue >= extrMin );
    assert( *pvalue <= extrMax );
    PerturbatorBase::addPerturbation();
#ifdef DEBUG4  
    cout << "value changed " << oldValue << " -> " << *pvalue << endl;
#endif    
    return 0.0;
}
        
bool PerturbatorHelper::validatePerturbationValue( bool validation, double* pvalue ){
    lowBoundaryShock = false;
    highBoundaryShock = false;
    //if the perturbation is refused
    if ( !validation ){
        //restore the old value
        *pvalue = oldValue;
#ifdef DEBUG4
        cout << "Gaussian perturbation refused, "
             << "return to:   " << (*pvalue) << endl;
#endif
    }
#ifdef DEBUG4
    else{
        cout << "Gaussian perturbation accepted" << endl;
    }
#endif    
    return PerturbatorBase::validatePerturbation( validation );
}

void PerturbatorHelper::getExtrema( double& min, double& max ){
    min = extrMin;
    max = extrMax;
}

double PerturbatorHelper::getMean(){
    return mean;
}

void PerturbatorHelper::stopBurn(){
    PerturbatorBase::stopBurn();
    lowBoundaryShockCounter = 0;
    highBoundaryShockCounter = 0;
}


void PerturbatorHelper::getAllPriorParameters( vector<double>& params ) const{
    params.resize(hyperPriors.size());
    for (unsigned int i = 0; i < hyperPriors.size(); ++i){
        params[i] = priorParams[hyperPriors[i]];
    }
}
    
void PerturbatorHelper::setAllPriorParameters( const vector<double>& params){
    assert(params.size()==hyperPriors.size());
    for (unsigned int i = 0; i < hyperPriors.size(); ++i){
        priorParams[hyperPriors[i]] = params[i];
        hyperPriorsPert[i]->invalidate();
    }
    updateSavedLnValue();
}
    
unsigned int PerturbatorHelper::getNumberPriorParameters() const {
    return hyperPriors.size();
}


void PerturbatorHelper::updateSavedLnValue(){
    switch (priorType){
        case EXP:
            savedLnValue = log(priorParams[0]);
        break;
        case NORMAL:
            savedLnValue = -log(priorParams[1]) - M_LN_SQRT_2PI;
        break;
        case LOGNORMAL:
            savedLnValue = -log(priorParams[0]) - M_LN_SQRT_2PI;
        break;
        case GAMMA:
            savedLnValue = -log(priorParams[1]) - lgamma(priorParams[0]);
        break;
        default: break;
    }
}

        
double PerturbatorHelper::getLnPriorValue( double value ) const{
#ifdef DEBUG2
    double test = savedLnValue;
    switch (priorType){
        case EXP:
            savedLnValue = log(priorParams[0]);
        break;
        case NORMAL:
            savedLnValue = -log(priorParams[1]) - M_LN_SQRT_2PI;
        break;
        case LOGNORMAL:
            savedLnValue = -log(priorParams[0]) - M_LN_SQRT_2PI;
        break;
        case GAMMA:
            savedLnValue = -log(priorParams[1]) - lgamma(priorParams[0]);
        break;
        default: break;
    }
#ifdef DEBUG4
   cout << "PerturbatorHelper::savedLnValue is " << test << endl;
   cout << "PerturbatorHelper::savedLnValue should be " << savedLnValue << endl;   
#endif
    assert(test == savedLnValue);    
#endif
    double calcVar1;
    double calcVar2;    
    switch (priorType){
        case FIXED: assert (value==extrMin); return 0.0;
        case UNIFORM: if ( (value>=extrMin) && (value<=extrMax) ){
                          return savedLnValue;
                      }
                      else {
                          cerr << "Warning, proposing value out of the range "
                               << "defined by the uniform prior on " << name
                               << ": " << value << endl;
                          return -DBL_MAX;
                      }
        case EXP:
            return savedLnValue - (priorParams[0]*value);
        case NORMAL:
            calcVar1 = (value - priorParams[0]);
            return  -(calcVar1*calcVar1)/(2.0*priorParams[1]*priorParams[1])
                    + savedLnValue;
        case LOGNORMAL: 
            calcVar1 = log(value-priorParams[2]);
            calcVar2 = calcVar1-log(priorParams[1]);
            calcVar2 = calcVar2*calcVar2;
            return -(calcVar2/(2.0*(priorParams[0]*priorParams[0]))) - calcVar1 + savedLnValue;
        case GAMMA:
            calcVar1 = (value-priorParams[2])/priorParams[1];
            return ((priorParams[0]-1.0)*(log(calcVar1))) - calcVar1 +savedLnValue;
        default:
            assert(0);
    }
    assert(0);
    return 0.0;
}
