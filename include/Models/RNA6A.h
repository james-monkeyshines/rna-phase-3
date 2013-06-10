#ifndef RNA6A_H
#define RNA6A_H

#include <math.h>

#include "RnaModel.h"

#include "Sequence/SequenceTable.h"
#include "Util/array2D.h"


// The general time reversible model for nucleotide substitution with
// using the discrete gamma model with some sites invariant


class RNA6A : public RnaModel {
protected:

    /** ***********************************************************************
     * RNA6A
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA General Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA6A( const string & registrationName );

    /** ************************************************************************
     * RNA6A
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA7A substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA6A( ParametersSet & parameters );

public:

    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const;

    /** ************************************************************************
     * initialisation
     * @input      sequenceTable, the SequenceTable instance containing the
     *             species' sequences used to perform the inference
     *             (not compulsory)
     * @semantics  Initialise the model and its parameters according to
     *             the sequence table if provided.
     ************************************************************************ */
    virtual void initialisation( SequenceTable * sequenceTable = NULL, int modelId = 0 );

    /** ************************************************************************
     * ~RNA6A
     * @semantics  Destructor of the class RNA6A, free dynamic variables
     ************************************************************************ */
    virtual ~RNA6A();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

    /** ************************************************************************
     * getState
     * @input   stateNumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int stateNumber, unsigned int symbolCategory = 0 ) const;


protected:
    void initEquivalencyTable();
    virtual void initMatrixIndex();
            
private:
    static RNA6A prototype;

};

#endif //RNA6A_H




