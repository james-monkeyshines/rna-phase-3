#ifndef RNA16J_H
#define RNA16J_H

#include <math.h>

#include "RNA16.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"


class RNA16J : public RNA16 {
protected:

    /** ***********************************************************************
     * RNA16J
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA16J Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16J( const string & registrationName );

    /** ************************************************************************
     * RNA16J
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16J substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16J( ParametersSet & parameters );

public:

    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA16J( parameters );
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
     * ~RNA16J
     * @semantics  Destructor of the class RNA16J, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16J();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

protected:
    virtual void initMatrixIndex();
            
private:
    static RNA16J prototype;

};

#endif //RNA16J_H




