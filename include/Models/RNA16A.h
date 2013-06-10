#ifndef RNA16A_H
#define RNA16A_H

#include <math.h>

#include "RNA16.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"


class RNA16A : public RNA16 {
protected:

    /** ***********************************************************************
     * RNA16A
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA16A Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16A( const string & registrationName );

    /** ************************************************************************
     * RNA16A
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16A substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16A( ParametersSet & parameters );

public:

    /* ***********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA16A( parameters );
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
     * ~RNA16A
     * @semantics  Destructor of the class Rna6a, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16A();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

protected:
    virtual void initMatrixIndex();
            
private:
    static RNA16A prototype;

};

#endif //RNA16A_H




