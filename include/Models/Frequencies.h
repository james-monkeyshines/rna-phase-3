#ifndef FREQUENCIES_H
#define FREQUENCIES_H

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

class Frequencies : public Observer< UpdateMessage >, public Subject< UpdateMessage >{
protected:
    Frequencies() {}

    /* a helper function to initialise prior and pertubation parameters
     * for S or X when used. Called in initialiseMCMC                      */
    void registerVariationParameter( ParametersSet& parameters, Perturbator* perturbator );

public:
    Frequencies( ParametersSet& parameters, MatrixModel* model );

    void initPrimitive(ParametersSet& parameters,
            unsigned int numberGammaCategories, unsigned int invariant );
            
    inline unsigned char getGammaCatAdjust(){
        return gammaCatAdjust;
    }
    
    inline unsigned int getNumberFrequenciesSets(){
        return numberFrequenciesSets;
    }
    
    inline unsigned int getNumberFrequenciesXSets(){
        return numberFrequenciesXSets;
    }
    
    inline const vector< double > & operator[](unsigned int rateCat ){
        return frequencies[frequenciesCategory[rateCat]];
    }
    
    inline const vector< double > & operator()(unsigned int frequenciesCat ){
        assert(frequenciesCat<numberFrequenciesSets);
        return frequencies[frequenciesCat];
    }
    
    inline unsigned int frequenciesCat( unsigned int rateCat ){
        return frequenciesCategory[rateCat];
    }
       
    void initialisation( unsigned int numberFrequencies );
    
    virtual void initialisation( const vector<double>& freq);
    
    double* linkS();    
    void assignS( double* sF );   
    vector< double > * linkX();  
    void assignX( vector<double>* xF );
    
    
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
    virtual void getModelParameters( ParametersSet & parameters, const Model* model, unsigned int symbolCategory ) const;
    
    /** ************************************************************************
     * setModelParameters
     ************************************************************************ */
    virtual void setModelParameters( ParametersSet & parameters, Model* model, unsigned int symbolCategory );
    
    /** ************************************************************************
     * printParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the parameters names and values of the model in a
     *             readable format.
     *             The frequencies and matrices cannot be handled at that level
     *             since their corresponding state is unknowm
     *             sF and xF values can be printed however
     ************************************************************************ */
    virtual void printParameters( ostream & outputStream ) const;
    
    /** ************************************************************************
     * initialiseMCMC
     * @input      perturbation/prior parameters
     * @input      the perturbator to register
     * @semantics  prepare the frequencies for a MCMC run
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
     * @semantics  please note that basic prior on frequencies
     *             (ie dirichlet prior) are handled by the perturbator given in
     *             parameter to initialiseMCMC. This function only deals with
     *             priors which are related to variation over rate categories
     ************************************************************************ */
    virtual double getLnPrior( const vector<double> & rates );

    
    /** ************************************************************************
     * update
     * @semantics  Observer method called by the Subject
     ************************************************************************ */
    virtual void update(UpdateMessage* subject);

    /** ************************************************************************
     * initialiseML
     * @input      penalty parameters for the frequencies
     * @semantics  prepare the penalty term in the likelihood for a penalized
     *             approach on frequencies variation among gamma categories
     ************************************************************************ */
    virtual void initialiseML( ParametersSet & parameters );

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
    virtual double getLnPenalty( const vector<double> & rates );

    /** ************************************************************************
     * diffLnPenalty
     * @input          none
     * @output         d(Lpenalty)/dparam, partial derivative vector wrt.
     *                 penalty free parameters
     * @semantics      this function returns diff(ln(penalty)) according to the user
     *                 prior specification (for ML runs)
     * @preconditions  ML system previously initialised
     ************************************************************************ */
    virtual void diffLnPenalty( vector<double>& gradVector ) const;

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
    void create( Perturbator* perturbator,
                      const string& name, const PriorField& prior,
                      unsigned int priorField, double& min, double& max, double& mean,
                      ParametersSet& parameters);
                      
    /** ************************************************************************
     * createSumAct
     * @semantics  helper method, read in prior the parameter for the 
     *             hyperparameter index priorField, create a hyperparameter
     *             (and register it to the perturbator) or store a constant value 
     ************************************************************************ */
    void createSumAct( Perturbator* perturbator,
                      const string& name, const PriorField& prior,
                      unsigned int priorField, double& min, double& max, double& mean,
                      ParametersSet& parameters);
        
    /** ************************************************************************
     * computeCovMatrix
     * @semantics  set up the covariance matrices for the numberNucleotides
     *             gaussian process. Should be updated each time
     *             hyperpriors are changing
     ************************************************************************ */
    void computeCovMatrix();

    /** ************************************************************************
     * computeLnPrior
     * @return     the new prior
     * @semantics  used to update the prior when hyperparameters or frequencies
     *             are changing
     ************************************************************************ */
    double computeLnPrior();
         
   
    // vector used to store the sets of frequencies 
    vector < double > * frequencies;

    // each gamma category can be assigned its own set of states frequencies
    // Usually numberFrequenciesSets == 1 (the usual case)
    unsigned int numberFrequenciesSets;
    
    // each model allows a different number of variable frequencies
    unsigned int numberFrequencies;
    
    // vector used to make the correspondance between a rate category
    // (invariant and gamma category) and the states frequencies.
    vector < int > frequenciesCategory;

    unsigned int invariant;
    unsigned int numberGammaCategories;
    
    //linear and "S" adjustment
    //frequencies adjustment : two frequencies sets only are used for the
    //lowest and highest category, other frequencies are deduced
    //according to their average substitution rate
    double* sF;

    // "X" adjustment
    // more frequencies and rate sets are added to model the frequencies
    // at (any) given rate(s)
    vector< double > * xF;
    vector < double > * frequenciesX;
    unsigned int numberFrequenciesXSets;
    
    unsigned char gammaCatAdjust;
    
    typedef enum{
        UNIFORM,
        GP,
        NB_TYPE
    } VariationType;
     
    VariationType variationType;
    
    
    //latent variable, prior parameter which is the sum of the activation
    //at each point
    vector< double > sumAct;
    vector< unsigned int > hyperSumAct;
    vector <PerturbatorParameter*> hyperSumActPert;
    
    //hyperparameters for each process
    vector< vector <double> > priorParams;
    //index of the hyperprior with relevant parameter
    vector< vector <unsigned int> > hyperPriors;
    vector< vector <PerturbatorParameter*> > hyperPriorsPert;
    //junk space storage for a(C,i)... C is the symbol, i is the rate index
    array2D< double > activationMatrix;
    
    //stored prior
    double lnPrior;
    //store the inverse of the covariance matrix of the gaussian process
    //associated with each frequency. (assumed to be equal for the moment)
    vector< array2D< double > > invCovMatrix;
    //store log ( det(CovMatrix) )
    vector< double > logdet;
    double logdetsum;
    //store the rates used for linear S and X and the last prior computation
    vector< double > savedRates;
};

#endif
