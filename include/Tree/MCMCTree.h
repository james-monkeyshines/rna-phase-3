#ifndef MCMCTREE_H
#define MCMCTREE_H

#include "Tree/InferenceTree.h"
#include "PatternDesign/Observer.h"

#include "Models/UpdateMessage.h"

#include <ostream>

class ParametersSet;

class MCMCTree : virtual public InferenceTree,
        public Observer< UpdateMessage >{
public:

    /** ************************************************************************
     * ~MCMCTree
     * @semantics  virtual destructor needed for inheritance
     ************************************************************************ */
    virtual ~MCMCTree(){};
                        
    /** ***********************************************************************
     * clone
     * @semantics  function used by the Factory<MCMCTree> to clone the prototype
     *             and create a real model.
    ************************************************************************ */
    virtual MCMCTree * clone( ParametersSet & parameters ) const = 0;    
    
    /** ************************************************************************
     * initialiseMCMC
     * @input      perturbation parameters
     * @semantics   initialise the perturbation parameter for this tree
     ************************************************************************ */
    virtual void initialiseMCMC( ParametersSet & parameters ) = 0;
            
    /** ************************************************************************
     * printPerturbationParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the perturbation parameters names and values of the
     *             model in a readable format
     ************************************************************************ */
    virtual void printPerturbationParameters( ostream & outputStream ) = 0;
    
    /** ************************************************************************
     * perturb
     * @return      ln(Hasting ratio) of the perturbation
     * @semantics   perturb the tree
     ************************************************************************ */
    virtual double perturb() = 0;

    /** ************************************************************************
     * initSampling
     * @input       a set of parameters
     * @semantics   initialise what is necessary to do the sample
     ************************************************************************ */
    virtual void initSampling(ParametersSet& parameters, bool overwrite) = 0;
    
    /** ************************************************************************
     * sample
     * @semantics   store the current state on files
     ************************************************************************ */
    virtual void sample() = 0;
    
    
     /** ************************************************************************
     * retrieveNodes
     * @semantics   come back to the last saved computation for the last
     *              category cat modified
     ************************************************************************ */
    virtual void retrieveNodesAux() = 0;
    
      /** ************************************************************************
     * saveNodes
     * @semantics   save the computation for the last category modified.
     ************************************************************************ */
    virtual void saveNodesAux() = 0;
  
    /** ************************************************************************
     * validatePerturbation
     * @input      boolean, yes if the perturbation is accepted
     * @semantics  return the tree to its initial state if the perturbation
     *             was not accepted
     ************************************************************************ */
    virtual bool validatePerturbation( bool validation ) = 0;
    
    /** ************************************************************************
     * stopBurn
     * @input      none
     * @semantics  this function is called to inform the tree that the burnin
     *             period is finished while performing a MCMC run
     ************************************************************************ */
    virtual void stopBurn() = 0;
    
    /** ************************************************************************
     * getAllParameters
     * @semantics    return all the parameters of the tree (in order to
     *               save/restore the state)
     ************************************************************************ */
    virtual void getAllParameters( vector<double>& params ) const = 0;
    
    /** ************************************************************************
     * setAllParameters
     * @input        the parameters to restore the state
     * @semantics    restore all the parameters of the tree with
     *               the given parameters
     ************************************************************************ */
    virtual void setAllParameters( const vector<double>& params) = 0;
    
    /** ************************************************************************
     * getNumberTreeParameters
     * @input      none
     * @semantics  return the number of free parameters in the tree (branch
     *             lengths, model for each branch, ...)
     ************************************************************************ */
     virtual unsigned int getNumberTreeParameters() const = 0;
     
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
     * getNumberPriorParameters
     * @semantics      return the number of free parameters used in the prior
     * @preconditions  the prior should have been initialised
     ************************************************************************ */
    virtual unsigned int getNumberPriorParameters() const = 0;

    /** ************************************************************************
     * getLnPrior
     * @input      none
     * @semantics  this function returns ln(prior) according to the user prior
     *             specification (for MCMC runs)
     ************************************************************************ */
    virtual double getLnPrior() const = 0;
    
};

#endif //MCMCTREE_H
