#include "Models/PerturbatorBase.h"

#include "Util/ParametersSet.h"

#include <assert.h>
#include <iostream>
#include <cstdlib>

PerturbatorBase::PerturbatorBase( const string& name,
    double minAcceptance, double maxAcceptance,
    double initialStep, double minStep, double maxStep):
    rand(Singleton < randombox >::instance()){

    assert(0);
    this->name = name;
    this->minAcceptance = minAcceptance;
    this->maxAcceptance = maxAcceptance;
    assert(minAcceptance >= 0.0);
    assert(maxAcceptance <= 1.0);
    assert(minAcceptance<maxAcceptance);
    if ( (minAcceptance == 0.0) && (maxAcceptance == 1.0) ){
        adaptativePerturbation = false;
    }
    else{
        adaptativePerturbation = true;
    }
    numberPerturbation = 0;
    numberAcceptedPerturbation = 0;
    this->step = initialStep;
    this->minStep = minStep;
    this->maxStep = maxStep;
    assert(initialStep<=maxStep);
    assert(minStep<=initialStep);
    stepModifier = 1.6;
    lastModification = 0;
}

PerturbatorBase::PerturbatorBase( const string& name,
                      ParametersSet& parameters, const string& genericName,
                      double initialStep, double minStep, double maxStep,
                      const string& initialStepString):
    rand(Singleton < randombox >::instance()){
    
    this->name = name;
    numberPerturbation = 0;
    numberAcceptedPerturbation = 0;
    this->step = initialStep;
    this->minStep = minStep;
    this->maxStep = maxStep;
    assert(initialStep<=maxStep);
    assert(minStep<=initialStep);
    stepModifier = 1.6;
    lastModification = 0;
    minAcceptance = .20;
    maxAcceptance = .25;
    //dynamic adaptation, acceptance rate range provided...
    if ( parameters.findParameter( genericName + ", proposal minimum acceptance rate") ){
        minAcceptance = parameters.doubleParameter
                    ( genericName + ", proposal minimum acceptance rate");
        maxAcceptance = parameters.doubleParameter
                    ( genericName + ", proposal maximum acceptance rate");
        if (minAcceptance < 0.0){
            cerr << "Proposal minimum acceptance rate lower than 0.0 for "
                 << genericName << ": " << minAcceptance << endl;
            exit(EXIT_FAILURE);
        }
        if (maxAcceptance > 1.0){
            cerr << "Proposal maximum acceptance rate greater than 1.0 for "
                 << genericName << ": " << maxAcceptance << endl;
            exit(EXIT_FAILURE);
        }
        if ( (minAcceptance>=maxAcceptance) ){
            cerr << "Error with the proposal acceptance rates for "
                 << genericName << ':' << endl
                 << "minimum acceptance rate = " << minAcceptance << endl
                 << "maximum acceptance rate = " << maxAcceptance << endl
                 << minAcceptance << " >= " << maxAcceptance << endl;
            exit(EXIT_FAILURE);
        }            
    }
    // if the dynamical adaptation is turned off, the tuning parameter
    // is compulsory (otherwise it is optional)
    if ( (minAcceptance == 0.0) && (maxAcceptance==1.0) ){
        adaptativePerturbation = false;
        step = parameters.doubleParameter
                    (genericName + ", " + initialStepString);
    }
    else{
        adaptativePerturbation = true;
        if ( parameters.findParameter(genericName + ", " + initialStepString) ){
            step = parameters.doubleParameter
                    (genericName + ", " + initialStepString);
        }
    }
    if ( ( step < minStep ) || ( step > maxStep ) ){
        cerr << genericName << " tuning/step parameter = " << step
             << "out of range [" << minStep << ", " << maxStep << ']' << endl;
        exit(EXIT_FAILURE);
    }
}


bool PerturbatorBase::validatePerturbation( bool validation ){
    if (validation){
        ++numberAcceptedPerturbation;
    }
    //if the dynamical step modification is on
    if( adaptativePerturbation ){
        if ( numberPerturbation == 200 ){
            double rate = (double)numberAcceptedPerturbation/
                          (double)numberPerturbation;
            if ( rate > maxAcceptance ){
                if (lastModification == -1){
                    //reduce the tuning modifier if close to the good tuning
                    //(when it jumps around the targeted acceptance rate)
//                    cout << "jump around : " << stepModifier;
                    stepModifier -= (stepModifier-1.0)*0.3;
//                    cout << "->" << stepModifier << endl;
                }
//                cout << name << " acceptance too high, step "
//                     << step << "->";
                step = step * stepModifier;
//                cout << step << endl;
                lastModification = 1;
            }
            else{
                if ( rate < minAcceptance ){
                    if (lastModification == 1){
                        //reduce the tuning modifier if close to the good tuning
                        //(when it jumps around the targeted acceptance rate)
//                        cout << "jump around : " << stepModifier;
                        stepModifier -= (stepModifier-1.0)*0.3;
//                        cout << "->" << stepModifier << endl;
                    }
//                    cout << name << " acceptance too low, step "
//                         << step << "->";
                    step = step / stepModifier;
//                    cout << step << endl;
                    lastModification = -1;
                }
            }
            numberAcceptedPerturbation = 0;
            numberPerturbation = 0;
            if (step>maxStep){
                step = maxStep;
//                cout << "max step reached ->" << step << endl;
            }
            if (step<minStep){
                step = minStep;
//                cout << "min step reached ->" << step << endl;
            }
        }
    }
    return true;    
}

double PerturbatorBase::addPerturbation(){
    ++numberPerturbation;
    return 0.0;
}

void PerturbatorBase::stopBurn(){
    adaptativePerturbation = false;
    numberPerturbation = 0;
    numberAcceptedPerturbation = 0;
}
        
        
PerturbatorBase::~PerturbatorBase(){
}


void PerturbatorBase::printPerturbationParameters( ostream& outputStream ) const{
   //redefined in PerturbatorFrequenciesParameter (to have a better sentence)
    outputStream << name << ", acceptance rate : "
                 << numberAcceptedPerturbation << '/' << numberPerturbation << "="
                 << (double)numberAcceptedPerturbation*100.0/(double)numberPerturbation
                 << "%, final step : "<< step << endl;
}

void PerturbatorBase::getAllPerturbationParameters( vector<double>& params ) const{
    params.resize(6);
    params[0] = (double)adaptativePerturbation;
    params[1] = (double)numberPerturbation;
    params[2] = (double)numberAcceptedPerturbation;
    params[3] = (double)lastModification;
    params[4] = (double)step;
    params[5] = (double)stepModifier;
}
    
void PerturbatorBase::setAllPerturbationParameters( const vector<double>& params){
    assert(params.size()==6);
    adaptativePerturbation = (bool)params[0];
    numberPerturbation = (int)params[1];
    numberAcceptedPerturbation = (int)params[2];
    lastModification = (int)params[3];
    step = (double)params[4];
    stepModifier = (double)params[5];
}
    
unsigned int PerturbatorBase::getNumberPerturbationParameters() const{
    return 6;
}

