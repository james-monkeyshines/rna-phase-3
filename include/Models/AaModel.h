#ifndef AAMODEL_H
#define AAMODEL_H

#include "Models/CodonParent.h"
class Perturbator;
class ParametersSet;
class SequenceTable;


class AaModel : public CodonParent {
  
protected:
    AaModel();
    
    AaModel( const string & registrationName );

    AaModel( ParametersSet & parameter );

    void initEquivalencyTable();
        
public:
    virtual ~AaModel();

    /** ************************************************************************
     * getNumberSymbols
     * @return  The number of distinct symbols recognized by the model
     ************************************************************************ */
    virtual inline unsigned int getNumberSymbols( unsigned int = 0 ) const {
        return (geneticCode.empty()) ? 22 : 22+CodonParent::getNumberSymbols();
    }
    
    /** ************************************************************************
     * getSymbolNumber
     * @input   base, a string representing a symbol of the model
     * @return  The basenumber of the symbol
     *           e.g. "Ala (A)"=0 , "Arg (R)"=1 ,... "Tyr (V)"=19
     ************************************************************************ */
    virtual int getSymbolNumber( const string & base, unsigned int = 0 ) const;

    /** ************************************************************************
     * getSymbol
     * @input   symbolNumber, the index of a symbol of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getSymbol( unsigned int symbolNumber, unsigned int = 0 ) const;

    /** ************************************************************************
     * getState
     * @input   stateNumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int stateNumber, unsigned int = 0 ) const;
    
    /** ************************************************************************
     * getNumberStates
     * @return  The number of states a model has
     ************************************************************************ */
    virtual inline unsigned int getNumberStates( unsigned int = 0 ) const {
        return 20;
    }
 
    /** ************************************************************************
     * retrieveEmpiricalFrequencies
     * @input      sequenceTable, the SequenceTable instance containing the
     *             species' sequences used to perform the inference.
     * @semantics  initialises the frequencies vector from the sequence
     ************************************************************************ */
    virtual vector<double> retrieveEmpiricalFrequencies( SequenceTable * sequenceTable, int modelId );


    /** ************************************************************************
     * retrieveEmpiricalRates
     * @input      sequenceTable, the SequenceTable instance containing the
     *             species' sequences used to perform the inference.
     * @input      a pointer to a double to store the proportion of invariant
     *             site (not compulsory)
     * @return     a rate ratio vector of length 5 initialised according to the
     *             sequences in the sequenceTable
     ************************************************************************ */
    virtual array2D<double> retrieveEmpiricalRates( SequenceTable * sequenceTable,
    int modelId, double * propInvariant = NULL );
    
    /** ************************************************************************
     * probability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    virtual double probability( unsigned int oldState, unsigned int newState, double time,
    unsigned int category, unsigned int symbolCategory = 0) const;

    /** ************************************************************************
     * diffProbability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    virtual double diffProbability( unsigned int oldState, unsigned int newState, double time,
    unsigned int category, unsigned int symbolCategory = 0 ) const;

    /** ************************************************************************
     * secondDiffProbability
     * @input   oldState, a state of the model
     * @input   newState, a state of the model
     * @input   time, a delay value
     * @input   category, a category
     * @return  The probability of ....
     ************************************************************************ */
    virtual double secondDiffProbability( unsigned int oldState, unsigned int newState, double time,
    unsigned int category, unsigned int symbolCategory = 0 ) const;

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
};

#endif //AAMODEL_H
