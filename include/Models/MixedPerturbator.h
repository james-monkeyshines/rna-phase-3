#ifndef MIXEDPERTURBATOR_H
#define MIXEDPERTURBATOR_H

#include <ostream>

#include <list>

#include "Util/randombox.h"

#include "PatternDesign/Singleton.h"

#include "Models/Perturbator.h"

#include "Models/Model.h"

class ParametersSet;

using namespace std;

/******************************************************************************
 * a class used to register frequencies and parameters
 * in order to perturb them during a MCMC run
 * the specificity of the mixed perturbator over the normal one
 * is that models can be parameter of the mixed perturbator, a perturbation
 * of the model calls model->perturb(), the mixed perturbator can perturb
 * the average substitution rate of the model
 * NB : ONE MUST REGISTER THE MODELS FIRST (then other parameters)
 **************************************************************************** */
class MixedPerturbator : public Perturbator {
protected:
    vector<Model*> models;
    vector<Model*> rndModels;
        
    int lastModifiedModel;
    
    vector<int> modelToModel;
    
public:
    /** ************************************************************************
     * MixedPerturbator
     * @semantics   constructor, initialise the multiple model perturbator
     ************************************************************************ */
    MixedPerturbator();

    /** ************************************************************************
     * ~MixedPerturbator
     * @semantics   destructor, free the memory
     ************************************************************************ */
    ~MixedPerturbator();
        
    /** ************************************************************************
     * registerModel
     * @input       model, a pointer to the model
     * @input       genericModelName, the beginning for the model label in
     *              parameters.
     * @input       parameters, a parameters set which contains the "rules"
     *              for the priors/perturbation of the average substitution
     *              rate and the model
     * @input       genericAverageRateName, the beginning for the labels in
     *              parameters should be something like "Average rates".
     *              For instance, a corresponding parameter could be
     *              "Average rates, minimum acceptance rate"
     * @input       averageRateInitialStep, the initial perturbation step for the
     *              average rate (the new value is drawn from a normal
     *              distribution centred at the current value, the
     *              initial standard deviation is that step and it may change
     *              during the burnin)
     * @input       the user can specify the initial step, a label
     *              genericName + ", " + initialStepString can be put in
     *              'parameters' for that.
     *              eg, "Average rates, initial step"
     * @semantics   register fully a new model in the mixed perturbator.
     *              NB: when a model is perturbed model->perturb() is
     *              called, therefore each model must have prepared its own
     *              prior/perturbation system
     ************************************************************************ */    
    void registerModel( Model* model, const string& genericModelName,
                       ParametersSet& parameters, const string& genericAverageRateName,
                       double* averageRate, double averageRateInitialStep,
                       const string& initialStepString );
            
    /** ************************************************************************
     * perturb
     * @semantics   first phase of a MCMC cycle: perturb some parameters
     *              we only perturb one parameter at a time usually
     *              when a model is perturbed, 
     * @output      int modelPerturbed, the model pertrubed
     *              int inModelPretrubed, the inside model perturbed
     *              heterogeneous::perturb( category, model )  -->
     *              mixedPerturbed( model, category ) -->
     *              mixed::perturb( category, junk ) -->
     *              mixedPerturbed( category, junk ) -->
     ************************************************************************ */
    double perturb();

    /** ************************************************************************
     * validatePerturbation
     * @input      validation, true if the change is accepted
     * @semantics  phase 2 of a MCMC cycle: accept or reject the change
     ************************************************************************ */
    bool validatePerturbation(bool validation);

    /** ************************************************************************
     * printParameters
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

};

#endif //MIXEDPERTURBATOR_H




