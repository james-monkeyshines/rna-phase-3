#ifndef TWOSTATEMODEL_H
#define TWOSTATEMODEL_H

#include "Models/MatrixModel.h"
class Perturbator;
class ParametersSet;
class SequenceTable;


class TwoStateModel : public MatrixModel {
  
protected:
    TwoStateModel( const string & registrationName );

    TwoStateModel( ParametersSet & parameter );
    
    void initEquivalencyTable();
    
public:
    virtual ~TwoStateModel();

    
    /** ************************************************************************
     * getSymbolNumber
     * @input   base, a string representing a symbol of the model
     * @return  The basenumber of the symbol
     *           e.g. "A"=0 , "C"=1 , "G"=2 and "U"=3 for the DNA models
     ************************************************************************ */
    int getSymbolNumber( const string & base, unsigned int symbolCategory = 0 ) const;


    /** ************************************************************************
     * getSymbol
     * @input   symbolNumber, the index of a symbol of the model
     * @return  The string which represents the state
     ************************************************************************ */
    string getSymbol( unsigned int symbolNumber, unsigned int symbolCategory = 0 ) const;

    /** ************************************************************************
     * getState
     * @input   stateNumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int stateNumber, unsigned int symbolCategory = 0 ) const;
    
    /** ************************************************************************
     * getNumberSymbols
     * @return  The number of distinct symbols recognized by the model
     ************************************************************************ */
    inline unsigned int getNumberSymbols( unsigned int = 0 ) const {
        return 82;
    }

    /** ************************************************************************
     * getNumberStates
     * @return  The number of states a model has
     ************************************************************************ */
    inline unsigned int getNumberStates( unsigned int = 0 ) const {
        return 2;
    }
 
    
    /** ************************************************************************
     * getExchangeability
     * @input  residue1, the index of a state of the model
     * @input  residue2, the index of a state of the model
     * @input  gammaCategory, a discrete gamma category
     * @return  The residue1<->residue2 exchangeability parameter for the given
     *          gamma category
     * @semantics  A two state model does not define matrixIndex
     *             (just one rate) and the base method has to be redefined
     ************************************************************************ */
    virtual double getExchangeability( unsigned int residue1,
        unsigned int residue2, unsigned int gammaCategory ) const;
                
        
    /** ************************************************************************
     * retrieveEmpiricalFrequencies
     * @input      sequenceTable, the SequenceTable instance containing the
     *             species' sequences used to perform the inference.
     * @semantics  initialises the frequencies vector from the sequence
     ************************************************************************ */
    vector < double > retrieveEmpiricalFrequencies( SequenceTable * sequenceTable, int modelId );


    /** ************************************************************************
     * getEmpiricalPropInvariant
     * @input      modelId
     * @return     a rate ratio vector of length 5 initialised according to the
     *             sequences in the sequenceTable
     ************************************************************************ */
    double getEmpiricalPropInvariant(SequenceTable* sequenceTable, int modelId );

        
    
    /** ************************************************************************
     * probability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    double probability( unsigned int oldState, unsigned int newState, double time,
    unsigned int category, unsigned int symbolCategory = 0) const;

    /** ************************************************************************
     * diffProbability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    double diffProbability( unsigned int oldState, unsigned int newState, double time,
    unsigned int category, unsigned int symbolCategory = 0 ) const;

    /** ************************************************************************
     * secondDiffProbability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    double secondDiffProbability( unsigned int oldState, unsigned int newState, double time,
    unsigned int category, unsigned int symbolCategory = 0 ) const;


    /** ************************************************************************
     * updateEigenMatrix
     * @input  update the eigen matrices when a parameter has changed
     ************************************************************************ */
    virtual void updateEigenMatrix();
    
    
    /** ************************************************************************
     * printParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the parameters names and values of the model in a
     *             readable format.
     *             -frequencies of the state
     *             -gamma shape parameter (if relevant)
     *             -proportion of invariant sites (if relevant)
     *             -substitution rate for each category
     ************************************************************************ */
    virtual void printParameters( ostream & outputStream ) const;
   
protected:
    virtual void setEigenMatrix() = 0;

};

#endif //TWOSTATEMODEL_H
