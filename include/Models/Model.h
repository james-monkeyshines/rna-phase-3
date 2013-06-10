#ifndef MODEL_H
#define MODEL_H

#include "configfix.h"

#include <string>
#include <fstream>
#include <vector>
#include <assert.h>

#include "Util/array2D.h"
#include "Util/ParametersSet.h"

#include "PatternDesign/Observer.h"
#include "Models/UpdateMessage.h"

class Perturbator;
class SequenceTable;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

class Model :  public Observer< UpdateMessage >, public Subject< UpdateMessage >{
public:

    /** ************************************************************************
     * ~Model
     * @semantics A virtual destructor needed for inheritance
     ************************************************************************ */
    virtual ~Model() {
    };

    /** ************************************************************************
     * initialisation
     * @semantics  initialise the model and its parameters according to
     *             default value or the content of the sequences
     ************************************************************************ */
    virtual void initialisation( SequenceTable * sequenceTable = NULL, int modelId = 0 ) = 0;

    /* ***********************************************************************
     * clone
     * @semantics  function used by the Factory<Model> to clone the prototype
     *             and create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const = 0;

    /** ************************************************************************
     * getModelParameters
     * @semantics  create a ParametersSet and fill it with the parameters of
     *             the model
     ************************************************************************ */
    virtual ParametersSet getModelParameters() const = 0;

    /** ************************************************************************
     * setModelParameters
     * @semantics      set the parameters of the model according to the content
     *                 of the parameters set
     * @preconditions  the parameters in the set are consistent
     * @postconditions issue a call to validChange() after. DO NOT FORGET
     ************************************************************************ */
    virtual void setModelParameters( ParametersSet & parameters ) = 0;

    /** ************************************************************************
     * validChange
     * @input  none
     * @semantics  to be called before using probability or prior if 
     *             the average substitution rate or setAllParameters has been
     *             called
     ************************************************************************ */
    virtual void validChange() = 0;
    
    /** ************************************************************************
     * getNumberStates
     * @return  The number of states a model has
     ************************************************************************ */
    virtual unsigned int getNumberStates( unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getState
     * @input   statenumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int statenumber, unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getFrequencyState
     * @input   freqStateNumber, the index of a frequency of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getFrequencyState( unsigned int freqStateNumber, unsigned int = 0 ) const = 0;
      
    /** ************************************************************************
     * getNumberSymbolCategory
     * @return  The number of categories of symbol the model deals with
     ************************************************************************ */
    virtual unsigned int getNumberSymbolCategory() const = 0;
    
        
    /** ************************************************************************
     * getNumberSymbols
     * @return  The number of distinct symbols recognized by the model
     ************************************************************************ */
    virtual unsigned int getNumberSymbols( unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getSymbolNumber
     * @input   base, a string representing a state of the model
     * @return  The basenumber of the state
     *           e.g. "A" = 1 , "C"="2" , "G"=3 and U="4" for the DNA models
     ************************************************************************ */
    virtual int getSymbolNumber( const string & base, unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getSymbol
     * @input   basenumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getSymbol( unsigned int basenumber, unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName() const = 0;

    /** ************************************************************************
     * getFrequency
     * @intput  residue, the index of a state of the model
     * @return  The frequencies of the given state of the model
     ************************************************************************ */
    virtual double getFrequency( unsigned int residue, unsigned int rateCategory,
                                 unsigned int symbolCategory ) const = 0;


    /** ************************************************************************
     * getFrequencies
     * @input   freq, an instance of vector <double>
     * @output  freq, an instance of vector <double> filled with the frequencies
     ************************************************************************ */
    virtual void getFrequencies( vector < double > & freq,
             unsigned int rateCategory, unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getNumberGammaCategories
     * @return  The number of gamma categories a model has
     ************************************************************************ */
    virtual unsigned int getNumberGammaCategories(  unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getInvariant
     * @return  invariantCategory?
     ************************************************************************ */
    virtual unsigned int getInvariant(  unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getNumberRatesCategories
     * @return  The number of rate categories a model has
     ************************************************************************ */
    virtual unsigned int getNumberRatesCategories(  unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * getRateCategoryProbability
     * @input   category, the index of a category
     * @return  The probability of a particular rate category
     ************************************************************************ */
    virtual double getRateCategoryProbability(  unsigned int category, unsigned  int symbolCategory ) const = 0;

    /** ************************************************************************
     * getAverageSubstitutionRate
     * @return  The average rate of substitution
     ************************************************************************ */
    virtual double getAverageSubstitutionRate(  unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * SetAverageSubstitutionRate
     * @semantics  Change the average rate of substitution for the given category
     * @postconditions issue a call to validChange() after
     ************************************************************************ */
    virtual void setAverageSubstitutionRate( double newAverageRate,  unsigned int symbolCategory ) = 0;


    /** ************************************************************************
     * getAllParameters
     * @input   pars, an instance of vector <double>
     * @output  pars, an instance of vector <double> filled with all the
     *          parameters of the model reparametrized so as to be independant
     ************************************************************************ */
    virtual void getAllParameters( vector < double > & pars ) const = 0;

    /** ************************************************************************
     * setAllParameters
     * @input  pars, an instance of vector <double>, filled with all the new
     *         parameters of the model
     * @postconditions issue a call to validChange() after
     ************************************************************************ */
    virtual void setAllParameters( const vector < double > & pars ) = 0;

    /** ************************************************************************
     * getNumberFreeParameters
     * @input   The number of free parameters in the model (and
     *          the size of the pars vector above
     ************************************************************************ */
    virtual unsigned int getNumberFreeParameters() const = 0;
    
    /** ************************************************************************
     * getNumberFreeFrequencyParameters
     * @return   The number of free frequency parameters in the model
     ************************************************************************ */
    virtual unsigned int getNumberFreeFrequencyParameters() const = 0;
	
    /** ************************************************************************
     * getOptimisableParameters
     * @input   empiricalFreqs, a boolean
	 * @inout	optimisableParameters, an instance of vector <unsigned int>
				that will be filled with all the indices of the parameters
				that need to be optimised.
     ************************************************************************ */
    virtual void getOptimisableParameters( bool empiricalFreqs, vector<unsigned int> &optimisableParameters ) const = 0;
	
    /** ************************************************************************
     * probability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    virtual double probability( unsigned int oldState, unsigned int newState, double time,
        unsigned  int category, unsigned int symbolCategory ) const = 0 ;

    /** ************************************************************************
     * diffProbability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    virtual double diffProbability( unsigned int oldState, unsigned int newState,
    double time, unsigned int category, unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * secondDiffProbability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    virtual double secondDiffProbability( unsigned int oldState,
    unsigned int newState, double time, unsigned int category, unsigned int symbolCategory ) const = 0;

    virtual const array2D< double > &
    getEquivalencyTable( unsigned int symbolCategory ) const = 0;

    /** ************************************************************************
     * printParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the parameters names and values of the model in a
     *             readable format
     ************************************************************************ */
    virtual void printParameters( ostream & outputStream ) const = 0;

    /** ************************************************************************
     * printLine
     * @input      outputStream  The stream used to output the parameters
     * @semantics  output the parameters values of the model in a single line
     ************************************************************************ */
    virtual void printLine( ostream & outputStream ) = 0;

    /** ************************************************************************
     * printPriorLine
     * @input      outputStream  The stream used to output the parameters
     * @semantics  output the parameters values of the priors in a single line
     ************************************************************************ */
    virtual void printPriorLine( ostream & outputStream ) = 0;
    
    /** ************************************************************************
     * fromLine
     * @input the new parameters of the model as they are printed with
     *        printLine (ie, same order)
     ************************************************************************ */
    virtual void fromLine( const vector<double>& parameters ) = 0;
            
    /** ************************************************************************
     * getNumberLineParameters
     * @input the number of parameters printed/read with printLine and fromLine
     ************************************************************************ */
    virtual unsigned int getNumberLineParameters() const = 0;
        
    /** ************************************************************************
     * initialiseMCMC
     * @input      perturbation/prior parameters
     * @semantics  Warn the model that it will have to perturb itself and
     *             provide some useful perturbation parameters
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
     * getAllPenaltyParameters
     * @output         a vector with the square root of the parameters
     * @semantics      get all the free parameters used in the parametric
     *                 penalized approach
     * @preconditions  initialiseML performed
     ************************************************************************ */
    virtual void getAllPenaltyParameters( vector<double>& params ) const = 0;
    
    /** ************************************************************************
     * setAllPenaltyParameters
     * @input          a vector with +/- square root of the new parameters
     * @semantics      change the parameters in the penalized approach
     * @preconditions  initialiseML performed
     ************************************************************************ */
    virtual void setAllPenaltyParameters( const vector<double>& params ) = 0;
    
    /** ************************************************************************
     * getNumberPenaltyParameters
     * @output         size of params in get/setAllPenaltyParameters
     * @preconditions  initialiseML performed
     ************************************************************************ */
    virtual unsigned int getNumberPenaltyParameters() const = 0; 
     
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
     *             specification (for MCMC runs)
     ************************************************************************ */
    virtual double getLnPrior() const = 0;
    
    
    
    /** ************************************************************************
     * initialiseML
     * @input      penalty parameters for the model
     * @semantics  add penalty terms in the likelihood for a penalized approach
     ************************************************************************ */
    virtual void initialiseML( ParametersSet & parameters ) = 0;
    
    /** ************************************************************************
     * diffLnPenalty
     * @input          none
     * @output         d(Lpenalty)/dparam, partial derivative vector wrt.
     *                 penalty free parameters
     * @semantics      this function returns diff(ln(penalty)) according to the user
     *                 prior specification (for ML runs)
     * @preconditions  initialiseML previously called
     ************************************************************************ */
    virtual void diffLnPenalty( vector<double>& gradVector) const = 0;
             
    /** ************************************************************************
     * getLnPenalty
     * @input      none
     * @semantics  this function returns ln(penalty) according to the user
     *             prior specification (for ML runs)
     ************************************************************************ */
    virtual double getLnPenalty() const = 0;
    
    
    /** ************************************************************************
     * perturb
     * @semantics  ask the model to perturb some of its parameters
     *             the flag cat is changed if a particular category is affected
     *             the flag node is changed if not all the nodes are affected
     ************************************************************************ */
    virtual double perturb() = 0;

    /** ************************************************************************
     * validatePerturbation
     * @input      validation status : true to validate the change from the
     *             last perturbation, false otherwise
     * @semantics  according to the validation flag, restore the model
     *             parameters as they were before the last call to perturb
     *             or validate the changes.
     ************************************************************************ */
    virtual bool validatePerturbation( bool validation ) = 0;
    
    /** ************************************************************************
     * stopBurn
     * @input      none
     * @semantics  this function is called to inform the model that the burnin
     *             period is finished while performing a MCMC run
     ************************************************************************ */
    virtual void stopBurn() = 0;    
    
        
    /** ************************************************************************
     * update
     * @semantics  Observer method called by the Subject
     ************************************************************************ */
    virtual void update(UpdateMessage* subject) = 0;


	virtual double getlnLAdjustmentEq() {
		return lnL_adj_eq_freq;
	}

	virtual double getlnLAdjustmentEmp() {
		return lnL_adj_emp_freq;
	}

    // average substitution rate do not access directly
    // in general case (use with perturbator)
    double averageSubstitutionRate;

	// 7-state models need a lnL adjustment value to be added to the final lnL,
	// in order to be compared to 4- and 16-state models.
	double lnL_adj_emp_freq;
	double lnL_adj_eq_freq;

};
#endif

