#ifndef PERTURBATORPARAMETER_H
#define PERTURBATORPARAMETER_H

#include <string>
#include <iostream>
using namespace std;

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Observer.h"
#include "Models/UpdateMessage.h"
#include "Util/randombox.h"

class ParametersSet;



/**
 * An ancestor class used to ease the perturbation of parameters during a
 * MCMC run. Define an interface for the parameter registered in a perturbator
 * This class is used mainly to handle the dynamic variation of the step
 * during the burnin period and to output perturbation values at the end
 * of the run. It cannot be instanciated (use PerturbatorHelper instead)
 */
class PerturbatorParameter: public Subject< UpdateMessage > {

protected:

public:
    /** ***********************************************************************
     * PerturbatorParameter
     * @semantics  empty default constructor
     *********************************************************************** */
    PerturbatorParameter(){};

    /** ***********************************************************************
     * PerturbatorParameter
     * @semantics register an observer for the parameter
     *********************************************************************** */
    PerturbatorParameter(Observer<UpdateMessage>& obs, UpdateMessage::ParamFlags message){
        attach(obs);
        clear();
        setType( UpdateMessage::PARAM_TYPE );
        setFlag( message );
    };

    /** ***********************************************************************
     * stopBurn
     * @semantics  end of the burnin period:
     *             -stop the dynamic tuning of the step parameter
     *             -reset counters to give acceptance rate in sampling period
     ************************************************************************ */
    virtual void stopBurn() = 0;


    /** ***********************************************************************
     * perturb
     * @return    ln Hasting Ratio
     * @semantics perturb the parameter
     ************************************************************************ */
    virtual double perturb(){
        notify();
        //the returned value is not meant to be used there
        return 0.0;
    }

    /** ***********************************************************************
     * validatePerturbation
     * @input      validation, true if the proposal was accepted
     ************************************************************************ */
    virtual bool validatePerturbation(bool validation){
        if (!validation){
            notify();
        }
        return true;
    }

    /** ***********************************************************************
     * printParameter
     * @input       the output stream
     * @semantics   print the final step and final acceptance rate of 'name'
     ************************************************************************ */
    virtual void printPerturbationParameters( ostream& outputStream ) const = 0;

    /** ***********************************************************************
     * ~PerturbatorParameter
     * @semantics   destructor
     ************************************************************************ */
    virtual ~PerturbatorParameter(){};

    /** ************************************************************************
     * getAllPerturbationParameters
     * @semantics      get all the parameters used in the perturbator
     *                 to save its state
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    virtual void getAllPerturbationParameters( vector<double>& params ) const = 0;

    /** ************************************************************************
     * setAllPerturbationParameters
     * @semantics      restore the state of the perturbator
     ************************************************************************ */
    virtual void setAllPerturbationParameters( const vector<double>& params) = 0;

    /** ************************************************************************
     * getNumberPerturbationParameters
     * @semantics      return the number of parameters used in the perturbator
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    virtual unsigned int getNumberPerturbationParameters() const = 0;

    /** ************************************************************************
     * getAllPriorParameters
     * @semantics      get all the parameters used in the priors
     *                 to save its state
     * @preconditions  the prior should have been initialised
     ************************************************************************ */
    virtual void getAllPriorParameters( vector<double>& params ) const = 0;

    /** ************************************************************************
     * setAllPriorParameters
     * @semantics      restore the state of the prior
     ************************************************************************ */
    virtual void setAllPriorParameters( const vector<double>& params) = 0;

    /** ************************************************************************
     * invalidate
     * @semantics      invalidate the previous computation of the prior
     ************************************************************************ */
    virtual void invalidate() = 0;

    /** ************************************************************************
     * getNumberPriorParameters
     * @semantics      return the number of free parameters used in the prior
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    virtual unsigned int getNumberPriorParameters() const = 0;

    /** ************************************************************************
     * getLnPrior
     * @input      none
     * @semantics  this function returns ln(prior) according to the user prior
     *             specification
     ************************************************************************ */
    virtual double getLnPrior() const = 0;

};

#endif //PERTURBATORPARAMETER_H
