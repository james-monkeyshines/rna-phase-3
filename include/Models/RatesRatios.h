#ifndef RATESRATIOS_H
#define RATESRATIOS_H

#include <vector>
#include <assert.h>
#include <iostream>

using namespace std;

#include "PatternDesign/Observer.h"
#include "Models/UpdateMessage.h"
#include "Util/array2D.h"

class Model;
class MatrixModel;
class Perturbator;
class PerturbatorParameter;
class PriorField;
class ParametersSet;

class RatesRatios : public Observer< UpdateMessage >, public Subject< UpdateMessage >{
protected:
    RatesRatios() {}

    /* a helper function to initialise prior and pertubation parameters
     * for S or X when used. Called in initialiseMCMC                      */
    void registerVariationParameter( ParametersSet& parameters, Perturbator* perturbator );
    
public:
    RatesRatios( ParametersSet& parameters, MatrixModel* model );
    
    void initPrimitive( ParametersSet& parameters,
          unsigned int numberGammaCategories, unsigned int invariant );
    
    inline unsigned char getGammaCatAdjust(){
        return gammaCatAdjust;
    }
    
    inline unsigned int getNumberRatesRatiosSets(){
        return numberRatesRatiosSets;
    }
    
    inline unsigned int getNumberRatesRatiosXSets(){
        return numberRatesRatiosXSets;
    }
    
    inline const vector< double > & operator[](unsigned int gammaCat ){
        return ratesRatios[ratesCategory[gammaCat]];
    }
    
    inline const vector< double > & operator()(unsigned int rateRatiosCat ){
        assert(rateRatiosCat<numberRatesRatiosSets);
        return ratesRatios[rateRatiosCat];
    }
    
    inline unsigned int ratesRatiosCat( unsigned int gammaCat ){
        return ratesCategory[gammaCat];
    }
    
    void initialisation( unsigned int numberRatesRatios );
    
    virtual void initialisation( const vector<double>& rates);

    double* linkS();    
    void assignS( double* sR );    
    vector<double> * linkX();    
    void assignX( vector<double>* xR );
    
    
    /** ************************************************************************
     * prepare
     * @input     the new vector of substitution rates of the model
     * @semantics called by the underlying MatrixModel when either the
     *            gamma shape parameter or the average substitution rate
     *            has changed. Update priors if type=var
     *            and values if L, LS or LX
     ************************************************************************ */
    virtual void prepare( const vector<double> & rates );
            
    /** ************************************************************************
     * updateLinear
     * @input     none
     * @semantics called by the previous function once savedRates has been
     *            updated when a L/LS or LX model is used. Update the
     *            values of the frequencies for intermediate gamma categories
     *            (all the values but the first and the last frequencies vector
     *            are invalid when the averageRate of the model has changed)
     ************************************************************************ */
    virtual void updateLinear();
    
    
    /** ************************************************************************
     * getModelParameters
     ************************************************************************ */
    virtual void getModelParameters( ParametersSet & parameters ) const;
    
    /** ************************************************************************
     * setModelParameters
     ************************************************************************ */
    virtual void setModelParameters( ParametersSet & parameters );
    
    /** ************************************************************************
     * printParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the parameters names and values of the model in a
     *             readable format.
     *             The frequencies and matrices cannot be handled at that level
     *             since their corresponding state is unknowm
     *             sF and xF values can be printed however
     ************************************************************************ */
    void printParameters( ostream & outputStream ) const;
    
    /** ************************************************************************
     * initialiseMCMC
     * @input      prior/perturbation parameters
     * @input      the perturbator to register
     * @semantics  prepare the rates ratios for a MCMC run
     ************************************************************************ */
    virtual void initialiseMCMC( ParametersSet & parameters, Perturbator* perturbator );
   
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
     ************************************************************************ */
    virtual void setAllParameters(const vector < double > & pars );

    /** ************************************************************************
     * getNumberFreeParameters
     * @return   The number of free parameters in the model (and
     *          the size of the pars vector above
     ************************************************************************ */
    virtual unsigned int getNumberFreeParameters() const;
    
    /** ************************************************************************
     * printLine
     * @input outputStream  The stream used to output the parameters
     * @semantics output the parameters values of the model in a single line
     ************************************************************************ */
    virtual void printLine( ostream & outputStream );

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
     * @return     lnPrior
     * @semantics  please note that basic prior on rates (ie uniform) are
     *             handled by the perturbator given in parameter to 
     *             initialiseMCMC. This function only deals with priors which
     *             are related to variation over rate categories
     ************************************************************************ */
    virtual double getLnPrior() const;


    /** ************************************************************************
     * update
     * @semantics  Observer method called by the Subject
     ************************************************************************ */
    virtual void update(UpdateMessage* subject);
    
    enum{
        L = 0x01,
        LS = 0x02,
        LX = 0x04,
        LINK = 0x08,
        LINKY = 0x10,
        LINKN = 0x20,
    };
        

protected: 
    /** ************************************************************************
     * create
     * @semantics  helper method, read in prior the parameter for the 
     *             hyperparameter index priorField, create a hyperparameter
     *             (and register it to the perturbator) or store a constant value 
     ************************************************************************ */
//    void create( Perturbator* perturbator,
//                      const string& name, const PriorField& prior,
//                      unsigned int priorField, double& min, double& max, double& mean,
//                      ParametersSet& parameters);
                      
    /** ************************************************************************
     * computeCovMatrix
     * @semantics  set up the covariance matrices for the numberNucleotides
     *             gaussian process. Should be updated each time
     *             hyperpriors are changing
     ************************************************************************ */
//    void computeCovMatrix();
        
    /** ************************************************************************
     * computeLnPrior
     * @return     the new prior
     * @semantics  used to update the prior when hyperparameters or frequencies
     *             are changing
     ************************************************************************ */
    double computeLnPrior();
    
    // vector used to store the sets of rateRatios 
    vector < double > * ratesRatios;

    // each gamma category can be assigned its own private rate matrix.
    // Usually numberRatesRatiosSets == 1 (the usual case)
    // each different rate matrix will be assigned its own rate ratios
    unsigned int numberRatesRatiosSets;
    
    // each model use a matrix with a given number of different rates
    // for the JukesCantor model, numberRatesRatios = 0
    // for the GTR model, numberRatesRatios = 5
    unsigned int numberRatesRatios;

    // vector used to make the correspondance between a rate category and its
    // rate matrix.
    vector < int > ratesCategory;

    unsigned int invariant;
    unsigned int numberGammaCategories;
    
    //linear and "S" adjustment
    //rates adjustment : two rate matrices only are used for the
    //lowest and highest category, other matrices are deduced
    //according to their average substitution rate
    double* sR;
    
    // "X" adjustment
    // more frequencies and rate sets are added to model the frequencies
    // at (any) given rate(s)
    vector< double > * xR;    
    vector < double > * ratesRatiosX;
    unsigned int numberRatesRatiosXSets;
    
    unsigned char gammaCatAdjust;
    
    
    typedef enum{
        UNIFORM,
        VAR,
        NB_TYPE
    } PriorType;
     
    PriorType priorType;
    //stored prior
    double lnPrior;
    
    //store the covariance matrix of the gaussian process associated with each
    //frequency. Symmetric matrix stored in compact format (size= n(n+1)/2)
    //WARNING, it is assumed that the matrix is definite positive
    //(not just semi-definite)
   // vector< vector< double > > covMatrix;
   // vector< vector< double > > invCovMatrix;
   // vector< double > logdet;
    //store the rates used for linear S and X and the last prior computation
    vector< double > savedRates;
    
    //hyperparameters for each process
    //vector< vector <double> > priorParams;
    //index of the hyperprior with relevant parameter
    //vector< vector <unsigned int> > hyperPriors;
    //vector< vector <PerturbatorParameter*> > hyperPriorsPert;
    //junk space storage for a(C,i)... C is the symbol, i is the rate index
    //array2D< double > activationMatrix;
   
    
};

#endif
