#ifndef RNA6B_H
#define RNA6B_H

#include <math.h>

#include "RNA6A.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"


// The general time reversible model for nucleotide substitution with
// using the discrete gamma model with some sites invariant


class RNA6B : public RNA6A {
protected:
    /** ***********************************************************************
     * RNA6B
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA 6B substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA6B( const string & registrationName );

    /** ************************************************************************
     * RNA6B
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA6B substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA6B( ParametersSet & parameters );

public:

    /* ***********************************************************************
     * clone
     * @input      a ParametersSet object to initialise the new model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  function used by the ModelFactory to clone the prototype and
     *             create a real model.
     ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA6B( parameters );
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
     * ~RNA6B
     * @semantics  Destructor of the class Rna7a, free dynamic variables
     ************************************************************************ */
    virtual ~RNA6B();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;


protected:

    virtual void initMatrixIndex();

private:
    static RNA6B prototype;

};

#endif //RNA6B_H




