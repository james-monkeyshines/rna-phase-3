#ifndef RNA6E_H
#define RNA6E_H

#include <math.h>

#include "RNA6A.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"



class RNA6E : public RNA6A {
protected:
    /** ***********************************************************************
     * RNA6E
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA6E substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA6E( const string & registrationName );

    /** ************************************************************************
     * RNA6E
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA6E substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA6E( ParametersSet & parameters );

public:
    /* ***********************************************************************
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
     * ~RNA6E
     * @semantics  Destructor of the class RNA6E, free dynamic variables
     ************************************************************************ */
    virtual ~RNA6E();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;


protected:

    virtual void initMatrixIndex();

private:
    static RNA6E prototype;

};

#endif //RNA6E_H




