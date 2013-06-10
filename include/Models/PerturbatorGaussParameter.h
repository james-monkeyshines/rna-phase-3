#ifndef PERTURBATORGAUSSPARAMETER_H
#define PERTURBATORGAUSSPARAMETER_H

#include "Util/randombox.h"
#include "PatternDesign/Singleton.h"
#include "Models/PerturbatorParameter.h"
#include "Models/PerturbatorHelper.h"

class Perturbator;

/*
 * a class to perturb a parameter, new values are drawn from a
 * gaussian distribution centred at the current values
 * this class can be registered in a perturbator
 */
 
class PerturbatorGaussParameter : public PerturbatorParameter, public PerturbatorHelper {
    
public:

    /** ***********************************************************************
     * PerturbatorGaussParameter
     * @input     name, the name of the parameter, used when we print the
     *            final perturbation step
     * @input     pvalue, a pointer to the parameter
     * @input     parameters, a ParametersSet which contains rules for
     *            the perturbation
     * @input     genericName the beggining of the label for the 'rules'
     *            in 'parameters'. Eg, "Rate ratios", "Gamma parameter",
     *            "Average Rates", ...
     * @input     initialStep, the DEFAULT initial standard deviation (see
     *            previous constructor)
     * @input     minStep, maxStep, the extrema value for the standard deviation
     * @input     min, max, the DEFAULT boundaries for the uniform
     *            prior distribution (see previous constructor)
     * @input     initialStepString, upperBoundString
     *            genericName + ", " + initialStepString
     *            genericName + ", " + upperBoundString are the label used
     *            in 'parameters' for user defined initialStep and max values
     *            if values are specified in parameters, they replace default
     *            ones
     * @semantics usual constructor
     ************************************************************************ */
    PerturbatorGaussParameter( Perturbator* pertubator, Observer<UpdateMessage>& obs,
                      UpdateMessage::ParamFlags message, const string& name, double* pvalue,
                      const PriorField& prior, ParametersSet& parameters, const string& genericName,
                      double initialStep, double minStep, double maxStep,
                      const string& initialStepString );
    
    /** ***********************************************************************
     * perturb
     * @return     ln(Hasting ratio) of the proposal (should be 0.0)
     * @semantics  propose a new value for the parameter *pvalue
     *             newValue = oldValue + N(0,standard deviation)
     *             newValue is reflected to make sure it remains in
     *             [min, max]
     ************************************************************************ */
    virtual double perturb();

    /** ***********************************************************************
     * stopBurn
     * @semantics  end of the burnin period:
     *             -stop the dynamic tuning of the step parameter
     *             -reset counters to give acceptance rate in sampling period
     *             This call from PerturbatorParameter if sent to the base
     *             PerturbatorHelper::stopBurn() method
     ************************************************************************ */
    virtual void stopBurn();
    
    /** ***********************************************************************
     * validatePerturbation
     * @input      validation, true if the change is accepted, false if it is
     *             rejected
     * @return     true
     * @semantics  restore the old value or valid the new one.
     *             During burnin time, if the user did not prevent change in
     *             the step, the standard deviation can grow or shrink to
     *             decrease or increase the acceptance rate. The SD remains
     *             in [minStep, maxStep] in any case.
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
    
    /** *************************************************************************
     * getExtrema
     * @return    the boundaries for the parameters
     ************************************************************************** */
    inline void getExtrema( double& min, double& max ){
        PerturbatorHelper::getExtrema( min, max );
    }
    
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
     * updateLnPrior
     * @input      none
     * @semantics  called when hyperpriors are changing or after perturbation
     ************************************************************************ */
    virtual double computeLnPrior() const;

    double* pvalue;

    //to save the prior
    double lnPrior;
    double oldLnPrior;
};

#endif //PERTURBATORGAUSSPARAMETER_H
