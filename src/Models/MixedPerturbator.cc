#include "Models/MixedPerturbator.h"

#include "Models/PerturbatorParameter.h"
#include "Models/PerturbatorGaussParameter.h"

#include "Util/ParametersSet.h"

#include <iostream>
#include <assert.h>

using namespace std;

MixedPerturbator::MixedPerturbator():Perturbator(){
    lastModifiedModel = -1;
    models.clear();
    rndModels.clear();
}

MixedPerturbator::~MixedPerturbator(){
}

void MixedPerturbator::registerModel( Model* model, const string& genericModelName,
                       ParametersSet& parameters, const string& genericAverageRateName,
                       double* averageRate, double averageRateInitialStep,
                       const string& initialStepString ){
    unsigned int modelPrio = parameters.intParameter(genericModelName + ", proposal priority");
    if (models.size()){
        double extrMin, extrMax;
        char label[60];
        sprintf( label, "Model%d/model1 average substitution rate ratio",
                        models.size()+1 );
        //when its average rate is modified, the underlying model will receive a message AVERAGE_RATE
        //it has to treat it
        PerturbatorGaussParameter* pert =
                registerGaussParameter( *model, UpdateMessage::AVERAGE_RATE,
                          label, averageRate,
                          parameters, genericAverageRateName,
                          averageRateInitialStep, 0.001, 5.0,
                          initialStepString );
        pert->getExtrema( extrMin, extrMax );
        if (extrMin<0.0){
            cerr << "invalid prior for the " << genericAverageRateName
                 << parameters.stringParameter( genericAverageRateName + ", prior") << endl;
            exit(EXIT_FAILURE);
        }
    }
    modelToModel.insert(modelToModel.end(), modelPrio, models.size() );
    models.push_back( model );
    rndModels.insert( rndModels.end(), modelPrio, model );
}
        
double MixedPerturbator::perturb(){
    double res = 0.0;
    if( rndTable.size() || rndModels.size() ){
        assert(lastModified == -1);
        assert(lastModifiedModel == -1);
        int nb = (int)(rand.ran()*(double)(rndTable.size()+rndModels.size()));
        if ( nb >= (int)rndTable.size() ){
            lastModifiedModel = nb - (int)rndTable.size();
            res = rndModels[lastModifiedModel]->perturb();
        }
        else{
            lastModified = nb;
            //if we change an average substitution rate
            res = rndTable[lastModified]->perturb();
        }
    }
    return res;
}

bool MixedPerturbator::validatePerturbation(bool validation){
    if( rndTable.size() || rndModels.size() ){
        if (lastModified != -1){
            assert(lastModifiedModel==-1);
            bool res = rndTable[lastModified]->
                                 validatePerturbation(validation);
            lastModified = -1;
            return res;
        }
        else{
            assert(lastModifiedModel!=-1);
            bool res = rndModels[lastModifiedModel]->validatePerturbation(validation);
            lastModifiedModel = -1;
            return res;
        }
    }
    return true;
}

void MixedPerturbator::printPerturbationParameters(ostream& outputStream) {
    Perturbator::printPerturbationParameters( outputStream );
    for ( unsigned int i = 0; i < models.size(); ++i ){
        outputStream << "MODEL " << i+1 << endl;
        models[i]->printPerturbationParameters( outputStream );
    }
}


void MixedPerturbator::getAllPerturbationParameters( vector<double>& params ) const{
    //first add the basic parameters
    Perturbator::getAllPerturbationParameters( params );
    //then add the parameters from each model
    vector<double> p;
    for (unsigned int j=0; j < models.size(); ++j){
        models[j]->getAllPerturbationParameters( p );
        params.insert( params.end(), p.begin(), p.end() );
    }
}
    
void MixedPerturbator::setAllPerturbationParameters( const vector<double>& params){
    assert( params.size() == getNumberPerturbationParameters() );
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp;
    //first set basic parameters
    nbp = Perturbator::getNumberPerturbationParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    Perturbator::setAllPerturbationParameters( p );
    //then for each model
    for (unsigned int j=0; j < models.size(); ++j){
        nbp = models[j]->getNumberPerturbationParameters();
        p.resize(nbp);
        for ( unsigned int i = 0; i < nbp; ++i ){
            p[i] = *iter;
            ++iter;
        }
        models[j]->setAllPerturbationParameters( p );
    }
}
    
unsigned int MixedPerturbator::getNumberPerturbationParameters() const{
    //count the basic number of parameters first (use ancestor method)
    unsigned int numberPerturbationParameters =
            Perturbator::getNumberPerturbationParameters();
    //then add the parameters from each model perturbator
    for (unsigned int j=0; j < models.size(); ++j){
        numberPerturbationParameters += models[j]->getNumberPerturbationParameters();
    }
    return numberPerturbationParameters;
}

void MixedPerturbator::getAllPriorParameters( vector<double>& params ) const{
    //first add the basic parameters
    Perturbator::getAllPriorParameters( params );
    //then add the parameters from each model
    vector<double> p;
    for (unsigned int j=0; j < models.size(); ++j){
        models[j]->getAllPriorParameters( p );
        params.insert( params.end(), p.begin(), p.end() );
    }
    assert( params.size() == getNumberPriorParameters() );
}
    
void MixedPerturbator::setAllPriorParameters( const vector<double>& params){
    assert( params.size() == getNumberPriorParameters() );
    vector<double>::const_iterator iter = params.begin();
    vector<double> p;
    unsigned int nbp;
    //first set basic parameters
    nbp = Perturbator::getNumberPriorParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    Perturbator::setAllPriorParameters( p );
    //then for each model
    for (unsigned int j=0; j < models.size(); ++j){
        nbp = models[j]->getNumberPriorParameters();
        p.resize(nbp);
        for ( unsigned int i = 0; i < nbp; ++i ){
            p[i] = *iter;
            ++iter;
        }
        models[j]->setAllPriorParameters( p );
    }
}
    
unsigned int MixedPerturbator::getNumberPriorParameters() const{
    //count the basic number of parameters first (use ancestor method)
    unsigned int numberPriorParameters =
            Perturbator::getNumberPriorParameters();
    //then add the parameters from each model perturbator
    for (unsigned int j=0; j < models.size(); ++j){
        numberPriorParameters += models[j]->getNumberPriorParameters();
    }
    return numberPriorParameters;
}
    
double MixedPerturbator::getLnPrior() const{
    double lnPrior = Perturbator::getLnPrior();
    for (unsigned int j=0; j < models.size(); ++j){
        lnPrior += models[j]->getLnPrior();
    }
    return lnPrior;
}

