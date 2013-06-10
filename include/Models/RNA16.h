#ifndef RNA16_H
#define RNA16_H

#include <math.h>
#include "Models/RnaModel.h"

class RNA16 : public RnaModel {
protected:

    /** ***********************************************************************
     * RNA16
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the GTR16 substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16( const string & registrationName );

    /** ************************************************************************
     * RNA16
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16 substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16( ParametersSet & parameters );

public:

    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA16( parameters );
    }

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
     * ~RNA16
     * @semantics  Destructor of the class Rna6a, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16();

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
    virtual void initMatrixIndex();

private:
    static RNA16 prototype;

};

#endif //RNA16_H

