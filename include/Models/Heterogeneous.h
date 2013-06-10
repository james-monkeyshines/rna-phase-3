#ifndef HETEROGENEOUS_H
#define HETEROGENEOUS_H

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <vector>

#include "Util/array2D.h"
#include "PatternDesign/Observer.h"
#include "Models/Model.h"
#include "Models/UpdateMessage.h"
#include "Models/Frequencies.h"
#include "Models/AncestralFrequencies.h"


class ParametersSet;
class MixedPerturbator;
using namespace std;

class Heterogeneous : public Model{
  
private:
    //array2D< double > equivalencyTable;
        
    vector< Model* > model;

    static Heterogeneous prototype;
    
protected:
    MixedPerturbator * perturbator;


private:
    /** ***********************************************************************
     * Heterogeneous
     * @semantics  default constructor of the mixed substitution
     *             model object, this construtor must be used only once
     *             to register the prototype to the ModelFactory
     ************************************************************************ */
    Heterogeneous( const string & registrationName );

public:
        
    vector<AncestralFrequencies*> ancestralFrequencies;

    /** ***********************************************************************
     * clone
     * @semantics  function used by the ModelFactory to clone the prototype and
     *             create a real model.
     ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new Heterogeneous( parameters );
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
    /** .***********************************************************************
     * Heterogeneous
     * @pre        parameters coherent
     * @semantics  Constructor of the Mixed substitution model object
     ************************************************************************ */
    Heterogeneous( ParametersSet & parameters );


public :
    void initialisation( SequenceTable * SequenceTable, int modelId = 0 );

    /** ************************************************************************
     * getNumberSymbolCategory
     * @return  The number of categories of symbol the model deals with
     ************************************************************************ */
    virtual inline unsigned int getNumberSymbolCategory() const {
        return model[0]->getNumberSymbolCategory();
    }

    /** ************************************************************************
     * ~Heterogeneous
     * @semantics    the destructor
     ************************************************************************ */
    virtual ~Heterogeneous();

    /** ************************************************************************
     * getNumberStates
     * @return  The number of states a model has
     ************************************************************************ */
    inline virtual unsigned int getNumberStates( unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory() );
        return model[0]->getNumberStates( symbolCategory );
    }
    
    /** ************************************************************************
     * getState
     * @input   statenumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int statenumber, unsigned int symbolCategory ) const {
        assert( symbolCategory < getNumberSymbolCategory() );
        return model[0]->getState( statenumber, symbolCategory );
    }
    
    /** ************************************************************************
     * getFrequencyState
     * @input       an index of frequency
     * @return      a string
     * @semantics   return a name to give to the frequency
     ************************************************************************ */
    virtual string getFrequencyState( unsigned int freqState, unsigned int symbolCategory ) const{
        assert( symbolCategory < getNumberSymbolCategory() );
        return model[0]->getFrequencyState( freqState, symbolCategory );
    }

    /** ************************************************************************
     * getNumberSymbols
     * @return  The number of distinct symbols recognized by the model
     ************************************************************************ */
    virtual unsigned int getNumberSymbols( unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory());
        return model[0]->getNumberSymbols( symbolCategory );
    }

    /** ************************************************************************
     * getSymbolNumber
     * @input   base, a string representing a state of the model
     * @return  The basenumber of the state
     *           e.g. "A" = 1 , "C"="2" , "G"=3 and U="4" for the DNA models
     ************************************************************************ */
    inline virtual int getSymbolNumber( const string & base, unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory() );
        return model[0]->getSymbolNumber( base, symbolCategory );
    }


    /** ************************************************************************
     * getSymbol
     * @input   basenumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    inline virtual string getSymbol( unsigned int baseNumber, unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory());
        return model[0]->getSymbol( baseNumber, symbolCategory );
    }
    
    /** ************************************************************************
     * getNumberModels
     * @return  The number of models
     ************************************************************************ */
    inline virtual unsigned int getNumberModels() const {
        return model.size();
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
     * getName
     * @return  The name of the model
     ************************************************************************ */
    inline virtual string getName( void ) const {
        char str[50];
        sprintf( str, "Heterogeneous model : %d * ", getNumberModels() );
        return string( string(str) + model[0]->getName() );
    }

    /** ************************************************************************
     * getFrequency
     * @intput  residue, the index of a state of the model
     * @return  The ancestral frequency of the given state of the model
     ************************************************************************ */
    virtual double getFrequency( unsigned int residue, unsigned int rateCategory,
                                        unsigned int symbolCategory ) const {
        return (*ancestralFrequencies[symbolCategory])[rateCategory][residue];
    }

    /** ************************************************************************
     * getFrequency
     * @intput  sbase, a state of the model
     * @return  The frequencies of the given state of the model
     ************************************************************************ */
    virtual double getFrequency( const string & sbase, unsigned int rateCategory,
                                        unsigned int symbolCategory ) const {
        return getFrequency( getSymbolNumber( sbase, symbolCategory ),
                             rateCategory, symbolCategory );
    }
    
    /** ************************************************************************
     * getFrequencies
     * @input   freq, an instance of vector <double>
     * @output  freq, an instance of vector <double> filled with the frequencies
     ************************************************************************ */
    virtual void getFrequencies( vector < double > & freq, unsigned int rateCategory, unsigned int symbolCategory ) const{
        cerr << "Bug, Heterogeneous::getFrequencies called" << endl;
        exit(1);
        freq = freq;
        rateCategory = symbolCategory;
    }

    /** ************************************************************************
     * getNumberGammaCategories
     * @return  The number of gamma categories a model has
     ************************************************************************ */
    virtual unsigned int getNumberGammaCategories(  unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory() );
        return model[0]->getNumberGammaCategories( symbolCategory );
    }

    /** ************************************************************************
     * getInvariant
     * @return  invariantCategory?
     ************************************************************************ */
    virtual unsigned int getInvariant(  unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory() );
        return model[0]->getInvariant( symbolCategory );
    }
    
    /** ************************************************************************
     * getNumberRatesCategories
     * @return  The number of rate categories a model has
     ************************************************************************ */
    virtual unsigned int getNumberRatesCategories( unsigned int symbolCategory ) const {
        assert(symbolCategory < getNumberSymbolCategory() );
        return model[0]->getNumberRatesCategories( symbolCategory );
    }

    /** ************************************************************************
     * getRateCategoryProbability
     * @input   category, the index of a category
     * @return  The probability of a particular rate category
     ************************************************************************ */
    inline virtual double getRateCategoryProbability( unsigned int category,
                                                   unsigned int symbolCategory ) const {
        return model[0]->getRateCategoryProbability( category, symbolCategory );
    }

    /** ************************************************************************
     * getAverageSubstitutionRate
     * @return  The average rate of substitution
     ************************************************************************ */
    inline virtual double getAverageSubstitutionRate( unsigned int symbolCategory ) const {
        cerr << "Bug, Heterogeneous::getAverageSubstitutionRate called" << endl;
        exit(1);
        return model[0]->getAverageSubstitutionRate( symbolCategory );
    }

    /** ************************************************************************
     * setAverageSubstitutionRate
     * @semantics  Provided to comply to the model interface. Must not be used
     ************************************************************************ */
    inline virtual void setAverageSubstitutionRate( double ,
                                                unsigned int  ){
        cerr << "Bug, Heterogeneous::setAverageSubstitutionRate called" << endl;
        exit(EXIT_FAILURE);
    }
 
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
     * @semantics forbiden call...
     ************************************************************************ */
    inline virtual double probability( unsigned int oldState, unsigned int newState,
                       double time, unsigned int category,unsigned int symbolCategory ) const {
        cerr << "You must use the Heterogeneous MCMC tree with Heterogeneous models" << endl;
        exit(1);
        oldState = newState;
        category=symbolCategory;
        time = 0.0;
        return -1.0;
    }

    /** ************************************************************************
     * diffProbability
     * @semantics forbiden call...
     ************************************************************************ */
    virtual double diffProbability( unsigned int oldState, unsigned int newState,
                       double time, unsigned int category, unsigned int symbolCategory ) const {
        cerr << "Bug, Heterogeneous::diffProbability called, check your control file" << endl;
        exit(1);
        oldState = newState;
        category=symbolCategory;
        time = 0.0;
        return -1.0;
    }

    /** ************************************************************************
     * secondDiffProbability
     * @semantics forbiden call...
     ************************************************************************ */
    virtual double secondDiffProbability( unsigned int oldState, unsigned int newState,
                      double time, unsigned int category, unsigned int symbolCategory ) const {
        cerr << "Bug, Heterogeneous::secondDiffProbability called, check your control file" << endl;
        exit(1);
        oldState = newState;
        category=symbolCategory;
        time = 0.0;
        return -1.0;
    }


    inline virtual const array2D< double > &
    getEquivalencyTable( unsigned int symbolCategory ) const {
        assert( symbolCategory < model.size() );
        return model[0]->getEquivalencyTable( symbolCategory );
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
     * @input      the set of prior/perturbation parameters for the
     *             perturbator (model specific)
     * @semantics  initialise the perturbation
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
     * @input      the message sent by a matrix or mixed model or the
     *             perturbator of ancestral frequencies
     * @semantics  when an underlying MatrixModel is updated we need to
     *             invalidate its branches and their ancestors. relay this
     *             message to the (heterogeneous) MCMC tree.
     ************************************************************************ */
    virtual void update( UpdateMessage* subject );
    
    /** ************************************************************************
     * validChange
     * @input  none
     * @semantics  to be called before using probability or prior if 
     *             the average substitution rate or setAllParameters has been
     *             called
     ************************************************************************ */
    virtual void validChange();
};

#endif //HETEROGENEOUS_H




