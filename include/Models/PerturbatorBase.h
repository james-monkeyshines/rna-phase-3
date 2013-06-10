#ifndef PERTURBATORBASE_H
#define PERTURBATORBASE_H

#include <string>
#include "Util/randombox.h"
#include "PatternDesign/Singleton.h"

class ParametersSet;

/**
 * A class used to ease the perturbation of parameters during a MCMC run
 * This class is used mainly to handle the dynamic variation of the step
 * during the burnin period and to output perturbation values at the end
 * of the run. These are the basic functionalities of the
 * PerturbatorGaussParameter
 */
class PerturbatorBase{
protected:  
    string name;
    Singleton<randombox>& rand;
    bool adaptativePerturbation;
    int numberPerturbation;
    int numberAcceptedPerturbation;
    
    double minAcceptance;
    double maxAcceptance;
    int lastModification;
    double step;
    double minStep;
    double maxStep;
    double stepModifier;
    
public:
    /** ***********************************************************************
     * PerturbatorBase
     * @input     name, the name of the parameter, used when we print the
     *            final perturbation step
     * @input     minAcceptance, maxAcceptance targeted acceptance range
     * @input     initialStep, the initial step. The step can change
     *            during the burnin time to get closer to the targeted
     *            acceptance rate
     * @input     minStep, maxStep, the extrema value for the standard deviation
     * @semantics simple constructor
     ************************************************************************ */
    PerturbatorBase( const string & name,
                       double minAcceptance, double maxAcceptance,
                       double initialStep, double minStep, double maxStep );

    /** ***********************************************************************
     * PerturbatorBase
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
    PerturbatorBase( const string& name,
                      ParametersSet& parameters, const string& genericName,
                      double initialStep, double minStep, double maxStep,
                      const string& initialStepString );
     
    /** ************************************************************************
     * ~PerturbatorBase
     * @return     0.0
     * @semantics  useless destructor
     ************************************************************************ */
     virtual ~PerturbatorBase();
             
    /** ************************************************************************
     * addPerturbation
     * @return     0.0
     * @semantics  add one to the counter of proposed perturbation
     ************************************************************************ */
    double addPerturbation();

    /** ************************************************************************
     * validatePerturbation
     * @input      validation, true if the proposal was accepted
     * @semantics  add one to the counter of accepted perturbation if
     *             validation == true. During the burnin period,
     *             if the adaptative step modification is not turned off,
     *             The step might be adapted each time the counter of proposed
     *             perturbation reaches 200, in order to have a correct
     *             accepted/proposed ratio.
     *             step is multiplied or divised by stepModifier if the
     *             acceptance is out of the specified range
     *             stepModifier converges to one when oscillations are
     *             spotted ( <minAcceptance then >maxAcceptance and vice-versa)
     ************************************************************************ */
    bool validatePerturbation( bool validation );
        
    /** ************************************************************************
     * getStep
     * @return     the step parameter
     * @semantics  since the PerturbatorHelper can be used just to have a
     *             dynamic step, we need a direct access to the step if we
     *             do not want to perturb a parameter with a normal
     *             distribution. When using the step to perturb something
     *             manually, one must use perturb and validatePerturbation
     *             in order to let the step adapt itself dynamically
     ************************************************************************ */
    inline double getStep(){
        return step;
    }
    
    /** ***********************************************************************
     * stopBurn
     * @semantics  end of the burnin period:
     *             -stop the dynamic tuning of the step parameter
     *             -reset counters to give acceptance rate in sampling period
     ************************************************************************ */
    void stopBurn();

    /** ***********************************************************************
     * printPerturbationParameter
     * @input       the output stream
     * @semantics   print the final step and final acceptance rate of 'name'
     ************************************************************************ */     
    void printPerturbationParameters( ostream& outputStream ) const;

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
    
};


#endif //PERTURBATORBASE_H
