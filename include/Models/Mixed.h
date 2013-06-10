#ifndef MIXED_H
#define MIXED_H

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Util/array2D.h"
#include "PatternDesign/Observer.h"
#include "Models/Model.h"
#include "Models/UpdateMessage.h"

class ParametersSet;
class SequenceTable;
class MixedPerturbator;

using namespace std;

class Mixed : public Model{
private:

    //array2D< double > equivalencyTable;

    vector< Model* > model;

    static Mixed prototype;
    
protected:
    MixedPerturbator * perturbator;


private:
    /** ***********************************************************************
     * Mixed
     * @semantics  default constructor of the mixed substitution
     *             model object, this construtor must be used only once
     *             to register the prototype to the ModelFactory
     ************************************************************************ */
    Mixed( const string & registrationName );

public:
    /** .**********************************************************************
     * clone
     * @semantics  function used by the ModelFactory to clone the prototype and
     *             create a real model.
     ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new Mixed( parameters );
    }

    /** ************************************************************************
     * getModelParameters
     * @semantics  create a ParametersSet and fill it with the parameters of
     *             the model
     ************************************************************************ */
    virtual ParametersSet getModelParameters() const;


    /** ************************************************************************
     * setModelParameters
     * @semantics      set the parameters of the model according to the content
     *                 of the parameters set
     * @preconditions  the parameters in the set are consistent
     ************************************************************************ */
    virtual void setModelParameters( ParametersSet & parameters );

protected:
    /** ************************************************************************
     * Mixed
     * @pre        parameters coherent
     * @semantics  Constructor of the Mixed substitution model object
     ************************************************************************ */
    Mixed( ParametersSet & parameters );


public :
    void initialisation( SequenceTable * SequenceTable, int modelId = 0 );

    /** ************************************************************************
     * getNumberSymbolCategory
     * @return  The number of categories of symbol the model deals with
     ************************************************************************ */
    virtual inline unsigned int getNumberSymbolCategory() const {
        return (unsigned int)model.size();
    }

    /** ************************************************************************
     * getModel
     * @input   the model id
     * @return  the model
     ************************************************************************ */
    inline virtual Model* getModel( unsigned int modelId ) const {
        return model[modelId];
    }

    /** ************************************************************************
     * ~Mixed
     * @semantics    the destructor
     ************************************************************************ */
     virtual ~Mixed();

    /** ************************************************************************
     * getNumberStates
     * @return  The number of states a model has
     ************************************************************************ */
    inline virtual unsigned int getNumberStates( unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getNumberStates( 0 );
    }
    
    /** ************************************************************************
     * getState
     * @input   stateNumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int stateNumber, unsigned int symbolCategory ) const{
        return model[symbolCategory]->getState( stateNumber, 0 );
    }

    /** ************************************************************************
     * getFrequencyState
     * @input       an index of frequency
     * @return      a string
     * @semantics   return a name to give to the frequency
     ************************************************************************ */
    virtual string getFrequencyState( unsigned int freqState, unsigned int symbolCategory ) const{
        return model[symbolCategory]->getFrequencyState( freqState, 0 );
    }

    /** ************************************************************************
     * getNumberSymbols
     * @return  The number of distinct symbols recognized by the model
     ************************************************************************ */
    virtual unsigned int getNumberSymbols( unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getNumberSymbols( 0 );
    }

    /** ************************************************************************
     * getSymbolNumber
     * @input   base, a string representing a state of the model
     * @return  The basenumber of the state
     *           e.g. "A" = 1 , "C"="2" , "G"=3 and U="4" for the DNA models
     ************************************************************************ */
    inline virtual int getSymbolNumber( const string & base, unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getSymbolNumber( base, 0 );
    }



    /** ************************************************************************
     * getSymbol
     * @input   basenumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    inline virtual string getSymbol( unsigned int baseNumber, unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getSymbol( baseNumber, 0 );
    }

    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    inline virtual string getName( void ) const {
        string name( "Mixed model : " + model[0]->getName() );
        for (unsigned int i = 1; i < model.size(); ++i ){
             name = name + " & " + model[i]->getName();
        }
        return name;
    }

    /** ************************************************************************
     * getFrequency
     * @intput  residue, the index of a state of the model
     * @return  The frequencies of the given state of the model
     ************************************************************************ */
    virtual double getFrequency( unsigned int residue, unsigned int rateCategory,
                                        unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getFrequency( residue, rateCategory, 0 );
    }

    /** ************************************************************************
     * getFrequency
     * @intput  sbase, a state of the model
     * @return  The frequencies of the given state of the model
     ************************************************************************ */
    virtual double getFrequency( const string & sbase, int rateCategory,
                                        unsigned int symbolCategory ) const {
        return getFrequency( getSymbolNumber( sbase, symbolCategory ),
                             rateCategory, symbolCategory );
    }
    
    /** ************************************************************************
     * getFrequencies
     * @input   freq, an instance of vector <double>
     * @output  freq, an instance of vector <double> filled with the frequencies
     ************************************************************************ */
    virtual void getFrequencies( vector < double > & freq,  unsigned int rateCategory, unsigned int symbolCategory ) const{
        model[symbolCategory]->getFrequencies( freq, rateCategory, 0 );
    }

    /** ************************************************************************
     * getNumberGammaCategories
     * @return  The number of gamma categories a model has
     ************************************************************************ */
    virtual unsigned int getNumberGammaCategories(  unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory() );
        return model[symbolCategory]->getNumberGammaCategories( 0 );
    }

    /** ************************************************************************
     * getInvariant
     * @return  invariantCategory?
     ************************************************************************ */
    virtual unsigned int getInvariant(  unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory() );
        return model[symbolCategory]->getInvariant( 0 );
    }
    
    /** ************************************************************************
     * getNumberRatesCategories
     * @return  The number of rate categories a model has
     ************************************************************************ */
    virtual unsigned int getNumberRatesCategories(  unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getNumberRatesCategories( 0 );
    }

    /** ************************************************************************
     * getRateCategoryProbability
     * @input   category, the index of a category
     * @return  The probability of a particular rate category
     ************************************************************************ */
    inline virtual double getRateCategoryProbability(  unsigned int category,
                                                    unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getRateCategoryProbability( category, 0 );
    }

    /** ************************************************************************
     * getAverageSubstitutionRate
     * @return  The average rate of substitution
     ************************************************************************ */
    inline virtual double getAverageSubstitutionRate(  unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getAverageSubstitutionRate( 0 );
    }

    /** ************************************************************************
     * setAverageSubstitutionRate
     * @semantics  Provided to comply to the model interface. Must not be used
     *             with symbolCategory != 0
     ************************************************************************ */
    virtual void setAverageSubstitutionRate( double newAverageRate,
                                                 unsigned int symbolCategory );
    
    /** ************************************************************************
     * validChange
     * @input  none
     * @semantics  to be called before using probability or prior if 
     *             the average substitution rate or setAllParameters has been
     *             called
     ************************************************************************ */
    virtual void validChange();

    /** ************************************************************************
     * getAllParameters
     * @input   pars, an instance of vector <double>
     * @output  pars, an instance of vector <double> filled with all the
     *          parameters of the model reparametrized so as to be independant
     ************************************************************************ */
    virtual void getAllParameters( vector < double > & pars ) const;

    /** ************************************************************************
     * setAllParameters
     * @input  pars, an instance of vector <double>, filled with all the new
     *         parameters of the model
     * @postconditions issue a call to validChange() after. DO NOT FORGET
     ************************************************************************ */
    virtual void setAllParameters( const vector < double > & pars );

    /** ************************************************************************
     * getNumberFreeParameters
     * @return   The number of free parameters in the model (and
     *          the size of the pars vector above
     ************************************************************************ */
    virtual unsigned int getNumberFreeParameters() const;
    
    /** ************************************************************************
     * getNumberFreeFrequencyParameters
     * @return   The number of free frequency parameters in the model
     ************************************************************************ */
    virtual unsigned int getNumberFreeFrequencyParameters() const;
    
    /** ************************************************************************
     * getOptimisableParameters
     * @input   empiricalFreqs, a boolean
	 * @inout	optimisableParameters, an instance of vector <unsigned int>
				that will be filled with all the indices of the parameters
				that need to be optimised.
     ************************************************************************ */
    virtual void getOptimisableParameters( bool empiricalFreqs, vector<unsigned int> &optimisableParameters ) const;
	
    /** ************************************************************************
     * probability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    inline virtual double probability( unsigned int oldState, unsigned int newState,
                       double time, unsigned int category,unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->probability( oldState, newState,
                                                           time, category, 0 );
    }

    /** ************************************************************************
     * diffProbability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    virtual double diffProbability( unsigned int oldState, unsigned int newState,
                       double time, unsigned int category, unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->diffProbability( oldState, newState,
                                                           time, category, 0 );
    }

    /** ************************************************************************
     * secondDiffProbability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    virtual double secondDiffProbability( unsigned int oldState, unsigned int newState,
                      double time, unsigned int category, unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->secondDiffProbability( oldState, newState,
                                                           time, category, 0 );
    }


    inline virtual const array2D< double > &
    getEquivalencyTable( unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[symbolCategory]->getEquivalencyTable( 0 );
    }

    /** ************************************************************************
     * printParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the parameters names and values of the model in a
     *             readable format
     ************************************************************************ */
    virtual void printParameters( ostream & outputStream ) const;

    /** ************************************************************************
     * printLine
     * @input      outputStream  The stream used to output the parameters
     * @semantics  output the parameters values of the model in a single line
     ************************************************************************ */
    virtual void printLine( ostream & outputStream );
    
    /** ************************************************************************
     * printPriorLine
     * @input      outputStream  The stream used to output the parameters
     * @semantics  output the parameters values of the priors in a single line
     ************************************************************************ */
    virtual void printPriorLine( ostream & outputStream );
    
    /** ************************************************************************
     * fromLine
     * @input the new parameters of the model as they are printed with
     *          printLine (ie, same order)
     ************************************************************************ */
    virtual void fromLine( const vector<double>& parameters );
            
    /** ************************************************************************
     * getNumberLineParameters
     * @return the number of parameters printed/read with printLine and
     *         fromLine
     ************************************************************************ */
    inline virtual unsigned int getNumberLineParameters() const;

    /** ************************************************************************
     * initialiseMCMC
     * @input      perturbation/prior parameters
     * @semantics  Warn the model that it will have to perturb itself and
     *             provide some useful perturbation parameters, call
     *             transmitted to internal models
     ************************************************************************ */
    virtual void initialiseMCMC( ParametersSet & parameters );
    

    /** ************************************************************************
     * printPerturbationParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the perturbation parameters names and values of the
     *             model in a readable format
     ************************************************************************ */
    virtual void printPerturbationParameters( ostream & outputStream );


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
     *             specification (for MCMC runs)
     ************************************************************************ */
    virtual double getLnPrior() const;
    
    /** ************************************************************************
     * getAllPenaltyParameters
     * @output         a vector with the square root of the parameters
     * @semantics      get all the free parameters used in the parametric
     *                 penalized approach
     * @preconditions  initialiseML performed
     ************************************************************************ */
    virtual void getAllPenaltyParameters( vector<double>& params ) const;
    
    /** ************************************************************************
     * setAllPenaltyParameters
     * @input          a vector with +/- square root of the new parameters
     * @semantics      change the parameters in the penalized approach
     * @preconditions  initialiseML performed
     ************************************************************************ */
    virtual void setAllPenaltyParameters( const vector<double>& params );
    
    /** ************************************************************************
     * getNumberPenaltyParameters
     * @output         size of params in get/setAllPenaltyParameters
     * @preconditions  initialiseML performed
     ************************************************************************ */
    virtual unsigned int getNumberPenaltyParameters() const;
    
    /** ************************************************************************
     * initialiseML
     * @input      penalty parameters for the model
     * @semantics  add penalty terms in the likelihood for a penalized approach
     ************************************************************************ */
    virtual void initialiseML( ParametersSet & parameters );
    
    /** ************************************************************************
     * diffLnPenalty
     * @input          none
     * @output         d(Lpenalty)/dparam, partial derivative vector wrt.
     *                 penalty free parameters
     * @semantics      this function returns diff(ln(penalty)) according to the user
     *                 prior specification (for ML runs)
     * @preconditions  initialiseML previously called
     ************************************************************************ */
    virtual void diffLnPenalty( vector<double>& gradVector ) const;
    
    /** ************************************************************************
     * getLnPenalty
     * @input      none
     * @semantics  this function returns -ln(prior) according to the user
     *             prior specification (for ML runs)
     ************************************************************************ */
    virtual double getLnPenalty() const;
    
    
    /** ************************************************************************
     * perturb
     * @semantics  ask the model to perturb some of its parameters
     *             the flag cat is changed if a particular category is affected
     *             the flag node is changed if not all the nodes are affected
     ************************************************************************ */
    virtual double perturb();
    
    /** ************************************************************************
     * validatePerturbation
     * @input      validation status : true to validate the change from the
     *             last perturbation, false otherwise
     * @semantics  according to the validation flag, restore the model
     *             parameters as they were before the last call to perturb
     *             or validate the changes.
     ************************************************************************ */
    virtual bool validatePerturbation( bool validation );

    /** ************************************************************************
     * stopBurn
     * @input      none
     * @semantics  this function is called to inform the model that the burnin
     *             period is finished while performing a MCMC run
     ************************************************************************ */
    virtual void stopBurn();

    /** ************************************************************************
     * update
     * @input      the message sent by an underlying model or average
     *             substitution rates
     * @semantics  when an underlying Mixed model is updated we need to
     *             invalidate its branches and their ancestors. relay this
     *             message to the (heterogeneous) MCMC tree so that it updates
     *             the given symbolCategory.
     ************************************************************************ */
    virtual void update( UpdateMessage* subject );


	virtual double getlnLAdjustmentEq() {
		double total = 0.0;
		bool required = false;
		for (unsigned int i = 0; i < model.size(); ++i ) {
			if (getNumberStates(i) == 7) {
				total += model[i]->getlnLAdjustmentEq();
				required = true;
			}
		}
		return required ? total : 1;
	}

	virtual double getlnLAdjustmentEmp() {
		double total = 0.0;
		bool required = false;
		for (unsigned int i = 0; i < model.size(); ++i ) {
			if (getNumberStates(i) == 7) {
				total += model[i]->getlnLAdjustmentEmp();
				required = true;
			}
		}
		return required ? total : 1;
	}

};

#endif //MIXED_H

