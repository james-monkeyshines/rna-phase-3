#include "Models/PerturbatorFrequenciesParameter.h"

#include <iostream>
#include <cstdlib>

#include "Util/statlib.h"

#include "Models/Perturbator.h"
#include "Models/PriorField.h"
#include "Models/PerturbatorGaussParameter.h"



PerturbatorFrequenciesParameter::PerturbatorFrequenciesParameter( Perturbator* perturbator, Observer<UpdateMessage>& obs,
      UpdateMessage::ParamFlags message, const string & name, vector<double>* pvalues,
      const PriorField& prior, ParametersSet& parameters, const string& genericName,
      double initialTuning, const string& initialTuningString ):
PerturbatorParameter(obs, message),        
PerturbatorBase(name, parameters, genericName, initialTuning, 100.0, 1e10, initialTuningString){

    this->pvalues = pvalues;
    //stepModifier < 1.0 this is a trick because for frequencies, the bigger
    //the tuning parameter the lower the change. With a stepModifier < 1.0 we
    //can use the PerturbatorBase adaptative step function for both
    //frequencies and gauss perturbation
    stepModifier = 0.6;
    priorType = NB_TYPE;
    hyperPriors.clear();
    savedLnValue = 0.0;
    if (prior.name() == "flat"){
        if (prior.getNumberParameters()){
            cerr << "invalid prior for frequencies: " << prior.toString() << endl;
            exit(EXIT_FAILURE);
        }
        priorType = FLAT;
        //p.d.f for a flat dirichlet prior is constant; D(O|1,1,...1) = (numberFrequencies-1)!
        //=>ln(D(O|1,1,1,...,1) =  sigma ln(x)        x=2..numberFrequencies-1)  
        for (double i = 2.0; i < pvalues->size() ; ++i ){
            savedLnValue += log(i);
        }
    }
    if (prior.name() == "dirichlet"){
        if (prior.getNumberParameters()==1){
            priorType = DIRICHLET_EQUAL;
            dirichletValues.resize(1);
            if (prior.isConstant(0)){
                //p.d.f for this dirichlet prior is function
                //of a constant; D(O|1,1,...1) = (numberFrequencies-1)!
                savedLnValue = lgamma(pvalues->size()*prior.getValue(0))-
                               pvalues->size()*lgamma(prior.getValue(0));
                dirichletValues[0] = prior.getValue(0);
            }
            else{
                double min, max;
                hyperPriors.push_back(0);
                dirichletValues[0] = 1.0;
                PerturbatorGaussParameter* pert =
                    perturbator->registerGaussParameter( *this, UpdateMessage::HYPER_PARAM, name+" hyperprior", &(dirichletValues[0]),
                        prior.getHyperPriorField(0), parameters, genericName+" hyperprior",
                        .5, .001, 5.0, "initial step" );
                hyperPriorsPert.push_back(pert);
                pert->getExtrema(min,max);
                if (min < 0.0){
                     cerr << "invalid range [" << min << ',' << max
                         << "] in a dirichlet prior... min < 0.0" << endl;
                    exit(EXIT_FAILURE);
                }
                dirichletValues[0] = pert->getMean();
            }
        }
        if (prior.getNumberParameters()==pvalues->size()){
            priorType = DIRICHLET;
            dirichletValues.resize(pvalues->size());
            hyperPriors.reserve(pvalues->size());
            hyperPriorsPert.reserve(pvalues->size());
            for ( unsigned int i = 0; i < pvalues->size(); ++i ){
                if (prior.isConstant(i)){
                    dirichletValues[i] = prior.getValue(i);
                    savedLnValue -= lgamma(dirichletValues[i]);
                    if (dirichletValues[i]<0.0){
                        cerr << "invalid negative parameter in a dirichlet prior"
                             << endl;
                        exit(EXIT_FAILURE);
                    }
                }
                else{
                    char paramName[35];
                    sprintf( paramName, " hyperprior for freq ");
                    char* endName = paramName + 21;
                    double min, max;
                    hyperPriors.push_back(i);
                    dirichletValues[i] = 1.0;
                    sprintf( endName, "%d", i+1);
                    string genName;
                    //set up the name used for the proposal probability
                    if (parameters.findParameter(genericName+paramName+", proposal priority")){
                        genName = genericName+paramName;
                    }
                    else{
                        genName =  genericName+" hyperprior";
                    }
                    PerturbatorGaussParameter* pert =
                        perturbator->registerGaussParameter( *this, UpdateMessage::HYPER_PARAM, name+paramName, &(dirichletValues[i]),
                            prior.getHyperPriorField(i), parameters, genName,
                            .5, .001, 5.0, "initial step" );
                    hyperPriorsPert.push_back(pert);
                    pert->getExtrema(min,max);
                    if (min < 0.0){
                        cerr << "invalid range [" << min << ',' << max
                             << "] in a dirichlet prior... min < 0.0" << endl;
                        exit(EXIT_FAILURE);
                    }
                    dirichletValues[i] = pert->getMean();
                }
            }
        }
        //unrecognized prior, wrong number of parameter
        if (priorType == NB_TYPE){
            cerr << "error, expecting 1 or " << pvalues->size()
                 << " parameters in the dirichlet prior on frequencies " << endl;
            exit(EXIT_FAILURE);
        }
    }
    //fixed prior
    if (prior.name() == "fixed"){
        priorType = FIXED;
        if (prior.getNumberParameters()!=0){
            if ( prior.getNumberParameters()!= pvalues->size() ){
                cerr << "error, expecting " << pvalues->size() << " parameters in the fixed prior on frequencies " << endl;
                exit(EXIT_FAILURE);
            }
            double tot = 0.0;
            for (unsigned int i = 0; i < pvalues->size(); ++i){
                if(!prior.isConstant(i)){
                    cerr << "error in the fixed prior on frequencies, element must be constant" << endl;
                    exit(EXIT_FAILURE);
                }
                else{
                    (*pvalues)[i] = prior.getValue(i);
                    tot += (*pvalues)[i];
                }
            }
            if (tot != 1.0){
                cerr << "error: fixed prior on frequencies does not sum to one" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    //unrecognized prior
    if (priorType == NB_TYPE){
        cerr << "unrecognized prior for frequencies: " << prior.name() << endl;
        exit(EXIT_FAILURE);
    }
    //initialise lnPrior;
    lnPrior = computeLnPrior();
}

bool PerturbatorFrequenciesParameter::validatePerturbation( bool validation ){
    //if the perturbation is refused
    if (!validation){
        //restore the old values and prior
        *pvalues = oldValues;
        lnPrior = oldLnPrior;
    }
#ifdef DEBUG4
    if (!validation){
        cout << "Frequencies, perturbation refused" << endl
             << "return to               : ";
        for( unsigned int i = 0; i < pvalues->size(); ++i ){
            cout << (*pvalues)[i] << "  ";
        }
        cout << endl;
    }
    else{
        cout << "Frequencies, perturbation accepted" << endl
             << "staying at              : ";
        for( unsigned int i = 0; i < pvalues->size(); ++i ){
            cout << (*pvalues)[i] << "  ";
        }
        cout << endl;
    }
#endif
    PerturbatorParameter::validatePerturbation( validation );
    return PerturbatorBase::validatePerturbation( validation );
}


double PerturbatorFrequenciesParameter::perturb(){
    
#ifdef DEBUG4  
    cout << "perturbator freq received a perturb order" << endl;
#endif    
    //save oldPrior
    oldLnPrior = lnPrior;
    
    if ( priorType == FIXED ){
        cerr << "Error, trying to modify the frequencies with a fixed prior... did you forget to set the proposal probability to 0?" << endl;
        exit(EXIT_FAILURE);
    }
    //store the old value in case the perturbation is not accepted
    oldValues = *pvalues;

    //modify the value   
    double numer = 0.0;
    double denom = 0.0;

    vector < double > dirichp1( oldValues.size() );
    vector < double > dirichp2( oldValues.size() );
    for ( unsigned int i = 0; i < oldValues.size(); ++i ) {
        dirichp1[i] = oldValues[i] * step;
    }
    rand.sdirichlet( dirichp1, *pvalues );
    double sum = 0.0;
    for ( unsigned int j = 0; j < (*pvalues).size(); ++j ){
        //to prevent instability...
        if ((*pvalues)[j]<1e-12){
            (*pvalues)[j] = 1e-12;
        }
        sum += (*pvalues)[j];
    }
    for ( unsigned int j = 0; j < (*pvalues).size(); ++j ){
        (*pvalues)[j] /= sum;
        dirichp2[j] = (*pvalues)[j] * step;
    }
    numer = statlib::ldirichlet( oldValues, dirichp2 );
    denom = statlib::ldirichlet( *pvalues, dirichp1 );
//#ifdef DEBUG4  
    cout << "frequencies changed from: ";
    for( unsigned int i = 0; i < oldValues.size(); ++i ){
        cout << oldValues[i] << "  ";
    }
    cout << endl;
    cout << "to                      : ";
    for( unsigned int i = 0; i < oldValues.size(); ++i ){
        cout << (*pvalues)[i] << "  ";
    }
    cout << endl;
    cout << "ln HR = " << numer-denom << endl;
//#endif    
    PerturbatorBase::addPerturbation();
    PerturbatorParameter::perturb();
    
    //compute newPrior
    lnPrior = computeLnPrior();
    return ( numer - denom );
}


void PerturbatorFrequenciesParameter::stopBurn(){
    PerturbatorBase::stopBurn();
}

void PerturbatorFrequenciesParameter::update(UpdateMessage* subject){
    //if hyperparam modification
    if ( !subject->hasFlag(UpdateMessage::PRIOR_FLAG) ){
        assert(subject->hasFlag(UpdateMessage::HYPER_PARAM));
        lnPrior = computeLnPrior();
        //warn for a change in the prior
        setFlag(PRIOR_FLAG);
        notify();
        unsetFlag(PRIOR_FLAG);
    }
}


void PerturbatorFrequenciesParameter::printPerturbationParameters( ostream& outputStream ) const{
    //overloaded to have a nice output
    outputStream << name << ", acceptance rate : "
                 << numberAcceptedPerturbation << '/' << numberPerturbation << "="
                 << (double)numberAcceptedPerturbation*100.0/(double)numberPerturbation
                 << "%, final dirichlet tuning parameter : "<< step << endl;
}


void PerturbatorFrequenciesParameter::getAllPerturbationParameters( vector<double>& params ) const{
    PerturbatorBase::getAllPerturbationParameters(params);
}

void PerturbatorFrequenciesParameter::setAllPerturbationParameters( const vector<double>& params){
    PerturbatorBase::setAllPerturbationParameters(params);
}
   
unsigned int PerturbatorFrequenciesParameter::getNumberPerturbationParameters() const{
    return PerturbatorBase::getNumberPerturbationParameters();
}

void PerturbatorFrequenciesParameter::getAllPriorParameters( vector<double>& params ) const{
    params.resize(hyperPriors.size());
    assert( !(hyperPriors.size()) || (priorType == DIRICHLET) || (priorType ==DIRICHLET_EQUAL) );
    for (unsigned int i = 0; i < hyperPriors.size(); ++i){
        params[i] = dirichletValues[hyperPriors[i]];
    }
}
    
void PerturbatorFrequenciesParameter::setAllPriorParameters( const vector<double>& params){
    assert(params.size()==hyperPriors.size());
    assert( !(hyperPriors.size()) || (priorType == DIRICHLET) || (priorType ==DIRICHLET_EQUAL) );
    for (unsigned int i = 0; i < hyperPriors.size(); ++i){
        dirichletValues[hyperPriors[i]] = params[i];
        hyperPriorsPert[i]->invalidate();
    }
    lnPrior=computeLnPrior();
}
    
unsigned int PerturbatorFrequenciesParameter::getNumberPriorParameters() const {
    return hyperPriors.size();
}

double PerturbatorFrequenciesParameter::getLnPrior() const{
#ifdef DEBUG1
    assert((float)lnPrior==(float)computeLnPrior());
#endif    
    return lnPrior;
}

double PerturbatorFrequenciesParameter::computeLnPrior() const{
    double partialLnDirichlet = 0.0;
    double sumDi = 0.0;
    
    if ( (priorType==DIRICHLET_EQUAL) || (priorType==DIRICHLET) ){
        sumDi = 0.0;
        partialLnDirichlet = 0.0;
        for (unsigned int i = 0; i < pvalues->size(); ++i){
            partialLnDirichlet +=  ( (priorType == DIRICHLET_EQUAL ? dirichletValues[0]-1 : dirichletValues[i]-1 ) * log((*pvalues)[i]) );
        }
    }
    switch (priorType){
        case FIXED: return 1.0;
        case FLAT: return savedLnValue;
        // otherwise, compute the term prod( freq_i ^ D_i-1 )
        case DIRICHLET_EQUAL:
            if (hyperPriors.size()==0){
                return savedLnValue + partialLnDirichlet;
            }
            else{
                return lgamma(pvalues->size()*dirichletValues[0])-
                               pvalues->size()*lgamma(dirichletValues[0]) + partialLnDirichlet;
            }
        case DIRICHLET:
            //first term: lgamma( sum (D_i) )
            for (unsigned int i = 0; i < pvalues->size(); ++i){
                sumDi += dirichletValues[i];
            }
            // some constant values in the term 1/Prod( gamma(D_i) ) have been included in priorLnValue but some might be missing
            for (unsigned int i = 0; i < hyperPriors.size(); ++i){
                partialLnDirichlet -= lgamma(dirichletValues[hyperPriors[i]]);
            }
            return (lgamma(sumDi) + savedLnValue + partialLnDirichlet);
        default: assert (0);
    }
    assert(0);
    return 1.0;
}

