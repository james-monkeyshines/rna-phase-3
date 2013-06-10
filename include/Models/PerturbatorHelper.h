#ifndef PERTURBATORHELPER_H
#define PERTURBATORHELPER_H

#include "Util/randombox.h"
#include "PatternDesign/Singleton.h"
#include "Models/PerturbatorBase.h"
#include "Models/Perturbator.h"
#include "Util/ParametersSet.h"

#include <string>

class PriorField;

/**
 * A class used to ease the perturbation of parameters during a MCMC run
 * This class uses PerturbatorBase and adds some functionality for the prior
 * The important difference with PerturbatorGaussParameter is that no
 * parameter is assigned to this object. It can be used with any of the
 * branch lengths for example (if they share the pertubation behaviour and
 * the prior)
 */
class PerturbatorHelper : public PerturbatorBase, public Observer<UpdateMessage>{
public:
        
    typedef enum{
        UNIFORM,
        NORMAL,
        LOGNORMAL,
        EXP,
        GAMMA,
        FIXED,
        NB_TYPE
    } PriorType;
        
    /** ************************************************************************
     * getPriorType
     * @input      none
     * @semantics  this function returns priorType, (ie FLAT, DIRICHLET, ...)
     ************************************************************************ */
    inline PriorType getPriorType(){ return priorType; }

    /** ***********************************************************************
     * PerturbatorHelper
     * @input     name, the name of the perturbator parameter,
     *            used when we print the final perturbation step
     * @input     parameters, a ParametersSet which contains rules for
     *            the perturbation
     * @input     genericName the beggining of the label for the 'rules'
     *            in 'parameters'. Eg, "Branch length", "Swap probability",
     * @input     initialStep, the DEFAULT initial step (see
     *            previous constructor)
     * @input     minStep, maxStep, the extrema value for the step
     * @input     initialStepString
     *            genericName + ", " + initialStepString is the label used
     *            in 'parameters' for user defined initialStep
     *            if a value is specified in parameters, it replaces the
     *            default one
     * @semantics usual constructor
     ************************************************************************ */
    PerturbatorHelper( Perturbator* perturbator, Observer<UpdateMessage>& obs, 
                      UpdateMessage::ParamFlags message, const string& name, const PriorField& prior,
                      ParametersSet& parameters, const string& genericName,
                      double initialStep, double minStep, double maxStep,
                      const string& initialStepString );
                      
                      
                      
    /** ***********************************************************************
     * ~PerturbatorHelper
     * @semantics virtual destructor
     ************************************************************************ */
     virtual ~PerturbatorHelper(){}
     
    /** ************************************************************************
     * perturbValue
     * @input      a pointer to the value to perturb
     * @input      [min, max] boundaries for that value
     * @input      stepMultiplier, a factor used to multiply the step
     * @semantics  *pvalue is replaced by a new value drawn from a normal
     *             distribution N(old value, step * stepMultiplier)
     *             PerturbatorParameter::perturb() is called in order to deal
     *             with the dynamic modification of the step during burnin
     *             One does not have to use such a perturbation mechanism
     *             and the perturb function defined in PerturbatorParameter
     *             can still be called directly if only interested by the
     *             adaptative step functionality.
     *             This function is used by PerturbatorGaussParameter
     ************************************************************************ */
    double perturbValue( double* pvalue );

    /** ************************************************************************
     * validatePerturbationValue
     * @input      validation, true if the proposal is accepted, false otherwise
     * @input      a pointer to the value which was perturbed (warning, no
     *             check done)
     * @semantics  the old value was stored during the call to perturbValue.
     *             if the proposal is rejected then pvalue is restored.
     *             PerturbatorParameter::validatePerturbation() is called
     *             to perform the dynamic modification of step
     *             One does not have to use such a validation mechanism
     *             and the validatePerturbation function defined in
     *             PerturbatorParameter must be called directly if only
     *             interested by the adaptative step functionality.
     *             This function is used by PerturbatorGaussParameter
     ************************************************************************ */
    bool validatePerturbationValue( bool validation, double* pvalue );
    
    
    /** ************************************************************************
     * getAllPriorParameters
     * @semantics      get all the parameters used in the priors
     *                 to save its state
     * @preconditions  the prior should have been initialised
     ************************************************************************ */
    void getAllPriorParameters( vector<double>& params ) const;
    
    /** ************************************************************************
     * setAllPriorParameters
     * @semantics      restore the state of the prior
     ************************************************************************ */
    void setAllPriorParameters( const vector<double>& params);
    
    /** ************************************************************************
     * getNumberPriorParameters
     * @semantics      return the number of free parameters used in the prior
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    unsigned int getNumberPriorParameters() const;
    
    /** ************************************************************************
     * getLnPriorValue
     * @input      the value to compute p(x)
     * @semantics  this function returns ln(prior) according to the user prior
     *             specification
     ************************************************************************ */
    double getLnPriorValue( double value ) const;
    
     /** ************************************************************************
     * updateSavedLnValue
     * @input      none
     * @semantics  when a priorParams is modified, part of the new prior can
     *             be computed without an actual reference to the value.
     *             this function update this constant part when necessary
     ************************************************************************ */
    void updateSavedLnValue();
   
    
    /** *************************************************************************
     * getExtrema
     * @return    the boundaries for the parameters
     ************************************************************************** */
    void getExtrema( double& min, double& max );

    
    /** *************************************************************************
     * getMean
     * @return    an "acceptable" mean value for the parameter
     ************************************************************************** */
    double getMean();
     
    
    /** *************************************************************************
     * stopBurn
     * @semantics    reinitialise some counters to have a count for the
     *               sampling period, stop the adaptive step modification
     *               ( with a call to PerturbatorBase::stopBurn() )
     ************************************************************************** */
    void stopBurn();
    
    /** ************************************************************************
     * update
     * @input      The message
     * @semantics  This function is called by hyperpriors when modified
     *             PerturbatorHelper react by updating savedLnValue (since
     *             it is expensive to compute it for each call to
     *             getLnPriorValue. PerturbatorGaussParameter is overriding it
     *             to update the saved lnPrior value simultaneously
     ************************************************************************ */
    virtual void update(UpdateMessage* subject){
        /* if this is more than a message to warn that a prior changed */
        /* ie, if priorParams has changed                              */
        if ( !subject->hasFlag(UpdateMessage::PRIOR_FLAG) ){
            updateSavedLnValue();
        }
        /* PeturbatorHelper does not send update message when its prior change */
    }
   
    /** *************************************************************************
     * getLowShock
     * @return    an "acceptable" mean value for the parameter
     ************************************************************************** */
    inline bool getLowShock(){
        return lowBoundaryShock;
    } 
    
    /** *************************************************************************
     * getHighShock
     * @return    an "acceptable" mean value for the parameter
     ************************************************************************** */
    inline bool getHighShock(){
        return highBoundaryShock;
    }
    
    /** *************************************************************************
     * getLowShockCounter
     * @return    an "acceptable" mean value for the parameter
     ************************************************************************** */
    inline unsigned int getLowShockCounter(){
        return lowBoundaryShockCounter;
    }
    
    
    /** *************************************************************************
     * getHighShockCounter
     * @return    an "acceptable" mean value for the parameter
     ************************************************************************** */
    inline unsigned int getHighShockCounter(){
        return highBoundaryShockCounter;
    }
    
       
protected:
    /** ************************************************************************
     * create
     * @input      name, the name for a parameter
     * @input      prior, the prior which contains info for that parameter
     * @input      priorField, index to use in prior
     * @inout      min/max, input:  min/max allowed values
     *                      output: min/max possible values
     * @output     mean, the mean value for the parameter, according to the prior
     * @input      parametersSet, necessary parameter for this prior and its
     *             hyperpriors
     * @input      generic name, name used to look in the parametersSet
     * @semantics  a generic function to add a new parameter and its prior/
     *             hyperpriors
     ************************************************************************ */
    void create( Perturbator* perturbator, Observer<UpdateMessage>& obs,
                      UpdateMessage::ParamFlags message, const string& name, const PriorField& prior,
                      unsigned int priorField, double& min, double& max, double& mean,
                      ParametersSet& parameters, const string& genericName);
                      
    double oldValue;

    double extrMin, extrMax, mean;
    
    bool lowBoundaryShock;
    bool highBoundaryShock;    
    
    unsigned int lowBoundaryShockCounter;
    unsigned int highBoundaryShockCounter;
    
#ifdef DEBUG2
    mutable
#endif
    double savedLnValue;
    
    
    //store hyperPriors parameters indices and perturbator
    vector< unsigned int > hyperPriors;
    vector< PerturbatorParameter* > hyperPriorsPert;
    vector< double > priorParams;
    PriorType priorType;
};

#endif //PERTURBATORHELPER_H
