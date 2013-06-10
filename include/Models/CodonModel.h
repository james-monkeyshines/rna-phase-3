#ifndef CODONMODEL_H
#define CODONMODEL_H

#include <vector>
#include <utility>
#include "Util/array2D.h"
#include "Models/CodonParent.h"


class Perturbator;
class ParametersSet;
class SequenceTable;

using namespace std;


class CodonModel : public CodonParent{

protected :
    /** ************************************************************************
     * CodonModel
     * @semantics  bypass the constructor of MatrixModel and CodonModel
     *
     ************************************************************************ */
    CodonModel();

    CodonModel( const string & registrationName );

    CodonModel( ParametersSet & parameter);


public :
    virtual ~CodonModel();

    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName() const;

    /** ************************************************************************
     * getNumberStates
     * @return  The number of states a model has
     ************************************************************************ */
    virtual unsigned int getNumberStates( unsigned int = 0) const{
        return geneticCode.size();
    }
    
    /** ************************************************************************
     * getState
     * @input   stateNumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int stateNumber, unsigned int = 0 ) const;

    /** ************************************************************************
     * retrieveEmpiricalFrequencies
     * @input      sequenceTable, the SequenceTable instance containing the
     *             species' sequences used to perform the inference.
     * @semantics  initialises the frequencies vector from the sequence
     ************************************************************************ */
    virtual vector< double > retrieveEmpiricalFrequencies( SequenceTable * sequenceTable, int modelId );

    /** ************************************************************************
     * retrieveEmpiricalRates
     * @input      sequenceTable, the SequenceTable instance containing the
     *             species' sequences used to perform the inference.
     * @input      a pointer to a double to store the proportion of invariant
     *             site (not compulsory)
     * @return     a rate ratio vector initialised according to the sequences
     *             in the sequenceTable
     ************************************************************************ */
    virtual array2D< double > retrieveEmpiricalRates( SequenceTable * sequenceTable,
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


protected:
    unsigned int matrixSize;

    virtual void setEigenMatrix();

    void initEquivalencyTable();

    string geneticCodeFileName;
};

#endif //CODONMODEL_H
