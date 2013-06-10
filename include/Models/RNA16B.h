#ifndef RNA16B_H
#define RNA16B_H

#include <math.h>

#include "RNA16.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"

/**
 * The F81-like paired-model with 15 free parameters (16 frequencies)
 */
class RNA16B : public RNA16 {
protected:

    /** ***********************************************************************
     * RNA16B
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA16HKY85 Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16B( const string & registrationName );

    /** ************************************************************************
     * RNA16B
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16HKY85 substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16B( ParametersSet & parameters );

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
     * ~RNA16B
     * @semantics  Destructor of the class RNA16B, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16B();
        
    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;


    /** ************************************************************************
     * printParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the parameters names and values of the model in a
     *             readable format.
     *             -frequencies of the state
     *             -gamma shape parameter (if relevant)
     *             -proportion of invariant sites (if relevant)
     ************************************************************************ */
    virtual void printParameters( ostream & outputStream ) const;
    
protected:
    virtual void initMatrixIndex();

private:
    static RNA16B prototype;
};

#endif //RNA16B_H




