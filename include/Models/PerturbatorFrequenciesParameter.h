#ifndef PERTURBATORFREQUENCIESPARAMETER_H
#define PERTURBATORFREQUENCIESPARAMETER_H

#include "Util/randombox.h"
#include "PatternDesign/Singleton.h"
#include "Models/PerturbatorParameter.h"
#include "Models/PerturbatorBase.h"


class Perturbator;
class PriorField;
/*
 * a class to perturb a set of frequencies (which sums to 1), new values are
 * drawn from a dirichlet distribution centred at the current value
 */
class PerturbatorFrequenciesParameter : public PerturbatorParameter, public PerturbatorBase, Observer<UpdateMessage> {        
public:

    typedef enum{
        FLAT,
        DIRICHLET_EQUAL,
        DIRICHLET,
        FIXED,
        NB_TYPE
    } PriorType;
        
    /** ************************************************************************
     * getPriorType
     * @input      none
     * @semantics  this function returns priorType, (ie FLAT, DIRICHLET, ...)
     ************************************************************************ */
    inline PriorType getPriorType() const { return priorType; }
        

    /** ***********************************************************************
     * PerturbatorFrequenciesParameter
     * @input     name, the name of the frequencies set, used when we print
     *            the final perturbation step
     * @input     pvalues, a pointer to the vector of frequencies
     * @input     parameters, a ParametersSet which contains rules for
     *            the perturbation
     * @input     genericName the beggining of the label for the 'rules'
     *            in 'parameters'. Eg, "Frequencies", "Ancestral frequencies",
     *            ...
     * @input     initialTuning, the DEFAULT initial tuning parameter (see
     *            previous constructor)
     * @input     initialStepString
     *            genericName + ", " + initialStepString is the label used
     *            in 'parameters' for user defined initialTuning value
     *            if a value is specified in parameters, it replace the
     *            default one.
     * @semantics usual constructor
     ************************************************************************ */
    PerturbatorFrequenciesParameter( Perturbator* perturbator,  Observer<UpdateMessage>& obs,
                      UpdateMessage::ParamFlags message, const string & name, vector<double>* pvalues,
                      const PriorField& prior, ParametersSet& parameters, const string& genericName,
                      double initialTuning, const string& initialTuningString );

    /** ***********************************************************************
     * perturb
     * @return     ln(Hasting ratio) of the proposal (non-trivial value)
     * @semantics  propose a new set of frequencies for the parameter *pvalues
     *             new vector = DirichletProp(old vector, tuning)
     ************************************************************************ */
    virtual double perturb();

    /** ***********************************************************************
     * stopBurn
     * @semantics  end of the burnin period:
     *             -stop the dynamic tuning of the step parameter
     *             -reset counters to give acceptance rate in sampling period
     *             This call from PerturbatorParameter if sent to the base
     *             PerturbatorBase::stopBurn() method
     ************************************************************************ */
    virtual void stopBurn();
    
    /** ***********************************************************************
     * validatePerturbation
     * @input      validation, true if the change is accepted, false if it is
     *             rejected
     * @return     true
     * @semantics  restore the old frequencies or valid the new ones.
     *             During burnin time, if the user did not prevent change,
     *             the Dirichlet tuning parameter can grow or shrink to
     *             increase or decrease the acceptance rate. However, it
     *             remains in the hard-coded interval [100.0, 100000.0]
     ************************************************************************ */
    virtual bool validatePerturbation( bool validation );

    /** ***********************************************************************
     * printParameter
     * @input       the output stream
     * @semantics   print the final step and final acceptance rate of 'name'
     ************************************************************************ */     
    virtual void printPerturbationParameters( ostream& outputStream ) const;

    /** ************************************************************************
     * getAllPerturbationParameters
     * @semantics      get all the parameters used in the perturbator
     *                 to save its state
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    virtual void getAllPerturbationParameters( vector<double>& params ) const;
    
    /** ************************************************************************
     * setAllPerturbationParameters
     * @semantics      restore the state of the perturbator
     ************************************************************************ */
    virtual void setAllPerturbationParameters( const vector<double>& params);
   
    /** ************************************************************************
     * getNumberPerturbationParameters
     * @semantics      return the number of parameters used in the perturbator
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    virtual unsigned int getNumberPerturbationParameters() const;
    
    /** ************************************************************************
     * getAllPriorParameters
     * @semantics      get all the parameters used in the priors
     *                 to save its state
     * @preconditions  the prior should have been initialised
     ************************************************************************ */
    virtual void getAllPriorParameters( vector<double>& params ) const;
    
    /** ************************************************************************
     * setAllPriorParameters
     * @semantics      restore the state of the prior
     ************************************************************************ */
    virtual void setAllPriorParameters( const vector<double>& params);
    
    /** ************************************************************************
     * getNumberPriorParameters
     * @semantics      return the number of free parameters used in the prior
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    virtual unsigned int getNumberPriorParameters() const;
    
    /** ************************************************************************
     * getLnPrior
     * @input      none
     * @semantics  this function returns ln(prior) according to the user prior
     *             specification
     ************************************************************************ */
    virtual double getLnPrior() const;
    
    /** ************************************************************************
     * update
     * @input      The message
     * @semantics  this function is called by hyperpriors when modified
     ************************************************************************ */
    virtual void update(UpdateMessage* subject);
    
     /** ************************************************************************
     * invalidate
     * @input      none
     * @semantics  this function MUST BE called when the parameter is modified
     *             with something else than perturb (to refresh the prior)
     ************************************************************************ */
    inline void invalidate(){
        lnPrior = computeLnPrior();
    }
   
private:
    /** ************************************************************************
     * computeLnPrior
     * @input      none
     * @semantics  called when hyperpriors are changing or after perturbation
     ************************************************************************ */
    virtual double computeLnPrior() const;
    
    vector<double> * pvalues;
    vector<double> oldValues;
    
    //store hyperPriors parameters indices and perturbator
    vector< unsigned int > hyperPriors;
    vector< PerturbatorParameter* > hyperPriorsPert;
    vector< double > dirichletValues;
    PriorType priorType;
    //to save a constant value when a flat prior is used
    double savedLnValue;
    
    //to save the prior
    double lnPrior;
    double oldLnPrior;
};

#endif //PERTURBATORFREQUENCIESPARAMETER_H
