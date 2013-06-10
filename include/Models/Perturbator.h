#ifndef PERTURBATOR_H
#define PERTURBATOR_H

#include <vector>
#include <fstream>

#include "Util/randombox.h"
#include "PatternDesign/Singleton.h"
#include "Models/PerturbatorParameter.h"
#include "Models/PriorField.h"

#include <stdarg.h>

using namespace std;

class ParametersSet;
class PerturbatorFrequenciesParameter;
class PerturbatorGaussParameter;

/******************************************************************************
 * a class used to register frequencies and parameters of
 * a model in order to perturb them during a MCMC run
 * frequencies are perturbed according to a Dirichlet proposal
 * new values for other parameters are drawn from a gaussian
 * distribution centred at their current value with a given
 * standard deviation (the standard deviation changes during the
 * burnin of the markov chain in order to reach a "good" acceptance
 * rate)
 **************************************************************************** */
class Perturbator {

protected:
    vector<PerturbatorParameter*> parameters;
    vector<PerturbatorParameter*> rndTable;
    Singleton < randombox > & rand;
    int lastModified;
    
public:
    /** ************************************************************************
     * Perturbator
     * @semantics   constructor, initialise the multiple model perturbator
     ************************************************************************ */
    Perturbator();

    /** ************************************************************************
     * registerFrequencies
     * @input       name, the name given to the set of frequencies (used when
     *              printing the values of perturbation parameters in
     *              printParameters)
     * @input       frequencies, a pointer to the frequencies of the model
     * @input       parameters, a parameters set which contains the "rules"
     *              for the perturbation of the frequencies
     * @input       genericName, the beginning for the labels in this ParametersSet.
     *              It should be something like "Frequencies".
     *              For instance, a corresponding label in 'parameters'
     *              could be "Frequencies, minimum acceptance rate"
     * @input       freqPrio, the priority for the perturbation of
     *              the frequencies wrt. to other frequencies and gaussian parameters
     *              registered
     * @input       defaultInitialFreqTuning, the initial tuning parameter for the
     *              Dirichlet perturbation of frequencies (the lower the value,
     *              the bigger the step). It may change during the burnin.
     * @input       the user can specify the initial tuning parameter, a label
     *              genericName + ", " + initialTuningString can be put in
     *              'parameters' for that.
     *              eg, "Frequencies, initial Dirichlet tuning parameter"
     * @semantics   register a set of frequencies (which sums to 1) in the
     *              perturbator
     ************************************************************************ */
    PerturbatorFrequenciesParameter* registerFrequencies( Observer<UpdateMessage>& obs, UpdateMessage::ParamFlags message,
        const string & name, vector<double> * frequencies,
        ParametersSet& parameters, const string & genericName,
        double defaultInitialFreqTuning, const string& initialTuningString );
        
        
    /** ************************************************************************
     * registerGaussParameter (two versions: normal parameter/prior parameter)
     * @input       name, the name given to the parameter (used when
     *              printing the values of perturbation parameters in
     *              printParameters)
     * @input       element, a pointer to the parameter of the model
     * @input       parameters, a ParametersSet which contains the "rules"
     *              for the perturbation.
     * @input(opt)  prior, a PriorField when using a hyperprior
     *              the prior field is read from the ParametersSet in the
     *              normal case
     * @input       genericName, the beginning for the labels in this ParametersSet.
     *              It should be something like "Rate ratios" or "Gamma parameter".
     *              For instance, a corresponding label in 'parameters'
     *              could be "Gamma parameter, minimum acceptance rate"
     * @input       elementPrio, the priority for the perturbation of
     *              this parameter wrt. to other frequencies and gaussian parameters
     *              registered
     * @input       initialStep, the initial standard deviation for the
     *              perturbation of the parameter. It may change during the burnin
     *              to get closer to a 'good' acceptance rate
     * @input       minStep,maxStep, the extrema values for the step.
     *              the step will not pass these boundaries during the run even if
     *              the targeted acceptance rate is not reached.
     * @input       the user can specify an initial standard deviation to replace
     *              initialStep the default one, a label
     *              genericName + ", " + initialStepString can be put in
     *              'parameters' for that.
     *              eg, "Invariant parameter, initial step"
     * @input       the user can specify an upper bound value to replace
     *              max the default one, a label
     *              genericName + ", " + upperBoundString can be put in
     *              'parameters' for that.
     *              eg, "Rate ratios, upper bound"
     * @input       defaultPrior, the default prior for the frequencies
     * @input       ... variable list of arguments for the prior
     * @semantics   register a set of frequencies (which sums to 1) in the
     *              perturbator
     ************************************************************************ */
    PerturbatorGaussParameter* registerGaussParameter( Observer<UpdateMessage>& obs, UpdateMessage::ParamFlags message,
        const string& name, double* element,
        ParametersSet& parameters, const string & genericName,
        double initialStep, double minStep, double maxStep,
        const string& initialStepString );
        
    PerturbatorGaussParameter* registerGaussParameter( Observer<UpdateMessage>& obs, UpdateMessage::ParamFlags message,
            const string& name, double* element,
        const PriorField& prior, ParametersSet& parameters, const string & genericName,
        double initialStep, double minStep, double maxStep,
        const string& initialStepString );

    /** ************************************************************************
     * perturb
     * @semantics   first phase of a MCMC cycle: perturb some parameters
     *              we only perturb one parameter at a time usually
     ************************************************************************ */
    double perturb();

    /** ************************************************************************
     * validatePerturbation
     * @input      validation, true if the change is accepted
     * @semantics  phase 2 of a MCMC cycle: accept or reject the change
     ************************************************************************ */
    bool validatePerturbation(bool validation);
    
    /** ************************************************************************
     * stopBurn
     * @semantics  call stopBurn for each parameter, it stops their adaptative
     *             perturbation and initialise the counter to compute the
     *             acceptance rate during the sampling period
     ************************************************************************ */
    void stopBurn();

    /** ************************************************************************
     * printPerturbationParameters
     * @input      outputStream, the stream used to print the parameters
     * @semantics  print the final values for the pertubartion parameters
     *             (not the parameters themselves), eg, final standard
     *             deviation, acceptance rate during the sampling period, ...
     ************************************************************************ */
    void printPerturbationParameters( ostream & outputStream );
        
    /** ************************************************************************
     * getAllPerturbationParameters
     * @semantics      get all the parameters used in the perturbator
     *                 to save its state
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    void getAllPerturbationParameters( vector<double>& params ) const;
    
    /** ************************************************************************
     * setAllPerturbationParameters
     * @semantics      restore the state of the perturbator
     ************************************************************************ */
    void setAllPerturbationParameters( const vector<double>& params);
    
    /** ************************************************************************
     * getNumberPerturbationParameters
     * @semantics      return the number of parameters used in the perturbator
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    unsigned int getNumberPerturbationParameters() const;

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
     * getLnPrior
     * @input      none
     * @semantics  this function returns ln(prior) according to the user prior
     *             specification
     ************************************************************************ */
    double getLnPrior() const;
    
    inline unsigned int getLoad(){
        return rndTable.size();
    }
};

#endif //PERTURBATOR_H




