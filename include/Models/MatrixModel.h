#ifndef MATRIXMODEL_H
#define MATRIXMODEL_H

#include "Models/Model.h"
#include "Models/RatesRatios.h"
#include "Models/Frequencies.h"

using namespace std;

class Perturbator;
class ParametersSet;
class SequenceTable;

class MatrixModel : public Model{

protected :
        
    MatrixModel();

    MatrixModel( const string & registrationName );

    MatrixModel( ParametersSet & parameters );
    
    void initInv( ParametersSet & parameters );
    void initGamma( ParametersSet & parameters );
    void initLink( ParametersSet & parameters );
    void initMatrix( unsigned int matrixSize );
    
    virtual void setEigenMatrix() = 0;
    
    /** ************************************************************************
     * outputFreqVector/outputMatrix
     * @input       labels, some labels (eg "GU", "A", "G+C")
     * @input       values, values to output
     * @input       outputStream, where to output
     * @input       pad, (padding) space to use (see code for more)
     * @input       prec, precision for the double (see code for more)
     * @semantics   simple functions to help output parameters in a readable
     *              format on the screen
     ************************************************************************ */
    void outputFreqVector( const vector<string>& labels, const vector<double>& values, ostream & outputStream, unsigned int pad, unsigned int prec ) const;
    void outputMatrix( const vector<string>& labels, const array2D<double>& values, ostream & outputStream, unsigned int pad, unsigned int prec ) const;
    
public :
    virtual ~MatrixModel();

    void initialisation( SequenceTable * sequenceTable = NULL, int modelId = 0 );
    
    /** ************************************************************************
     * getNumberSymbolCategory
     * @return  The number of categories of symbol the model deals with
     ************************************************************************ */
    virtual inline unsigned int getNumberSymbolCategory() const {
        return 1;
    }
    
    /** ************************************************************************
     * getSymbol
     * @input   symbolNumber, the index of a symbol of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getSymbol( unsigned int symbolNumber, unsigned int symbolCategory = 0 ) const = 0;

    /** ************************************************************************
     * getState
     * @input   stateNumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int stateNumber, unsigned int = 0 ) const = 0;
    
    /** ************************************************************************
     * getFrequencyState
     * @input   freqStateNumber, the index of a frequency of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getFrequencyState( unsigned int freqStateNumber, unsigned int = 0 ) const{
        assert(freqStateNumber<numberFrequencies);
        assert(numberFrequencies==getNumberStates());
        //by default, name of the state
        return getState( freqStateNumber );
    }

    /** ************************************************************************
     * getNumberStates
     * @input   stateNumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual unsigned int getNumberStates( unsigned int symbolCategory = 0 ) const = 0;
    
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


    
    /** ************************************************************************
     * getExchangeability
     * @input  residue1, the index of a state of the model
     * @input  residue2, the index of a state of the model
     * @input  gammaCategory, a discrete gamma category
     * @return  The residue1<->residue2 exchangeability parameter for the given
     *          gamma category
     ************************************************************************ */
    virtual double getExchangeability( unsigned int residue1,
        unsigned int residue2, unsigned int gammaCategory ) const;
    
    
    /** ************************************************************************
     * getFrequency
     * @intput  residue, the index of a state of the model
     * @intput  rateCategory, the discrete rate category
     * @intput  symbolCategory, must be 0
     * @return  The frequencies of the given state of the model
     ************************************************************************ */
    virtual double getFrequency( unsigned int residue,
        unsigned int rateCategory, unsigned int = 0 ) const {
        assert ( residue < (*frequencies)[rateCategory].size() );
        return (*frequencies)[rateCategory][residue];
    }


    /** ************************************************************************
     * getFrequencies
     * @input   freq, an instance of vector <double>
     * @output  freq, an instance of vector <double> filled with the frequencies
     ************************************************************************ */
    virtual void getFrequencies( vector < double > & freq, unsigned int rateCategory, unsigned int ) const{
        if( frequencies ){
            freq = (*frequencies)[rateCategory];
        }
        else{
            freq.clear();
        }
    }

    
    /** ************************************************************************
     * getNumberGammaCategories
     * @return  The number of gamma categories a model has
     ************************************************************************ */
    inline unsigned int getNumberGammaCategories( unsigned int = 0 ) const {
        return numberGammaCategories;
    }
    
    /** ************************************************************************
     * getInvariant
     * @return  invariantCategory?
     ************************************************************************ */
    inline unsigned int getInvariant( unsigned int = 0 ) const {
        return invariantCategory;
    }
    
    /** ************************************************************************
     * getNumberRatesCategories
     * @return  The number of rate categories a model has
     ************************************************************************ */
    inline unsigned int getNumberRatesCategories( unsigned int = 0 ) const {
        return numberRatesCategories;
    }
        
    /** ************************************************************************
     * getRateCategoryProbability
     * @input   category, the index of a category
     * @return  The probability of a particular rate category
     ************************************************************************ */
    double getRateCategoryProbability( unsigned int category, unsigned int symbolCategory = 0 ) const;

    /** ************************************************************************
     * getAverageSubstitutionRate
     * @return  The average rate of substitution
     ************************************************************************ */
    virtual double getAverageSubstitutionRate( unsigned int symbolCategory = 0 ) const;

    /** ************************************************************************
     * setAverageSubstitutionRate
     * @input  newAverageRate the new average rate of substitution
     ************************************************************************ */
    virtual void setAverageSubstitutionRate( double newAverageRate, unsigned int symbolCategory = 0);

    /** ************************************************************************
     * updateAverageRateVector
     * @input  none
     * @semantics  compute averageRate for each gamma category with the actual
     *             values of alpha and average substitution rate
     ************************************************************************ */
    virtual void updateAverageRateVector();
     
    /** ************************************************************************
     * updateEigenMatrix
     * @input  none
     * @semantics  compute the multiplier in front of each rate matrix to have
     *             the right averageRate for each of them. Call setEigenMatrix
     *             once it is done.
     ************************************************************************ */
    virtual void updateEigenMatrix();
    
    /** ************************************************************************
     * validChange
     * @input  none
     * @semantics  to be called before using probability or prior if 
     *             the average substitution rate or setAllParameters has been
     *             called
     ************************************************************************ */
    virtual void validChange();
                          
    /** ************************************************************************
     * initialiseMCMC
     * @input      perturbation parameters
     * @semantics  Warn the model that it will have to perturb itself and
     *             provide some useful perturbation parameters
     ************************************************************************ */
    virtual void initialiseMCMC( ParametersSet & parameters );

    /** ************************************************************************
     * perturb
     * @semantics  ask the model to perturb some of its parameters
     ************************************************************************ */
    virtual double perturb();

    /** ************************************************************************
     * validatePerturbation
     * @input      validation, true if the perturbation is accepted
     * @semantics  restore the model to the state before the perturbation
     *             if validation==false
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
     * getLnPenalty
     * @input      none
     * @semantics  this function returns -ln(prior) according to the user
     *             prior specification (for ML runs)
     ************************************************************************ */
    virtual double getLnPenalty() const;

    /** ************************************************************************
     * getAllParameters
     * @input   pars, an instance of vector <double>
     * @output  pars, an instance of vector <double> filled with all the
     *          parameters of the model reparametrized so as to be independant
     ************************************************************************ */
    virtual void getAllParameters(vector < double > & pars ) const;

    /** ************************************************************************
     * setAllParameters
     * @input  gamma, an instance of vector <double>, filled with all the new
     *         parameters of the model
     * @postconditions issue a call to validChange() after. DO NOT FORGET
     ************************************************************************ */
    virtual void setAllParameters(const vector < double > & pars );

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
     * printLine
     * @input outputStream  The stream used to output the parameters
     * @semantics output the parameters values of the model in a single line
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
     *        printLine (ie, same order)
     ************************************************************************ */
    virtual void fromLine( const vector<double>& parameters );
            
    /** ************************************************************************
     * getNumberLineParameters
     * @input the number of parameters printed/read with printLine and fromLine
     ************************************************************************ */
    virtual unsigned int getNumberLineParameters() const;
    
    /** ************************************************************************
     * printParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the parameters names and values of the model in a
     *             readable format. (the frequencies and matrices cannot be
     *             handled at that level
     *             -gamma shape parameter (if relevant)
     *             -proportion of invariant sites (if relevant)
     *             -substitution rate for each category
     ************************************************************************ */
    virtual void printParameters( ostream & outputStream ) const;

    /** ************************************************************************
     * printParametersAux
     * @input      outputStream, the stream used to output the parameters
     * @semantics  General method to print the exchangeability matrice
     *             and the Markov model matrix on the screen.
     *             Duplicated in RNA16 => update conjointly...
     ************************************************************************ */
    virtual void printParametersAux( ostream & outputStream ) const;
    
    
    /** ************************************************************************
     * getEquivalencyTable
     * @return  The equivalency table (symbol->state)
     ************************************************************************ */
    inline const array2D< double > & getEquivalencyTable( unsigned int = 0 ) const {
        return equivalencyTable;
    }
    
    /** ************************************************************************
     * errorSymbol
     * @input      the sequenceTable used to initialise the model
     * @input      the modelId used during initialisation
     * @input      speciesId, the species where the problem was found
     * @input      site, the position were the problem occured
     * @semantics  function called when the model find an unrecognized
     *             symbol when initialising from the sequence table.
     *             this function print an error message and abort.
     *             The symbol which caused the problem (modelId,speciesId,site)
     *             is converted back into a corresponding position in the data
     *             file to let the user know what is going wrong (hopefully)
     ************************************************************************ */
    void errorSymbol( SequenceTable * sequenceTable, unsigned int modelId,
                      unsigned int speciesId, int site);
                         
    /** ************************************************************************
     * getName
     * @return     The name of the model
     * @semantics  MatrixModel is not supposed to be instanciated but
     *             implements the basic mechanisms for gamma categories, ..
     *             This function is called by descendants to give the postfix
     *             of the model name (+dG4+I+...)
     ************************************************************************ */
    virtual string getName() const;
    
    /** ************************************************************************
     * update
     * @semantics  Observer method called by the Subject
     ************************************************************************ */
    virtual void update(UpdateMessage* subject);    
 
    

protected :

    Frequencies * frequencies;
    RatesRatios * ratesRatios;
    
    /* the matrixIndex is used to retrieve the exchangeability
     * parameter i<->j in the ratesRatios vector
     * some specific model do not use it and must redefine the function
     * getExchangeability in that case. Others have to fill this array2D*/
    array2D<int> matrixIndex;
    
    Perturbator * perturbator;

    unsigned int numberRatesCategories;

    //invariant category flag
    unsigned int invariantCategory;
    //proportion of invariant sites
    double proportionInvariantSites;    
    
    // discrete gamma model flag
    bool discreteGamma;
    // gamma shape parameter
    double alpha;
    // number of gamma categories
    unsigned int numberGammaCategories;
            
        
    array2D< double > equivalencyTable;
    
    // eigen system
    vector < double > * eigenValues;
    array2D< double > * eigenMatrix;
    array2D< double > * ieigenMatrix;
    //transition rate matrix
    array2D< double > * rateMatrix;
    /* average substitution rate for each transition matrice */
    vector < double > averageRate;
    /* multiplication factor used in front of each transition matrice */
    vector < double > substitutionRate;
    /* probability for each gamma categories
     * uniform if percentile gamma approximation */
    vector < double > categoriesWeights; 
    bool percentileGammaApproximation;

    unsigned int numberFrequencies;
    unsigned int numberRatesRatios;
};

#endif //MATRIXMODEL_H
