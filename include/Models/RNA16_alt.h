#ifndef RNA16_alt_H
#define RNA16_alt_H

#include "Models/RNA16.h"

class RNA16_alt : public RNA16 {
protected:

    /** ***********************************************************************
     * RNA16_alt
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the GTR16 substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16_alt( const string & registrationName );

    /** ************************************************************************
     * RNA16_alt
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16_alt substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16_alt( ParametersSet & parameters );

public:

    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA16_alt( parameters );
    }

    /** ************************************************************************
     * ~RNA16_alt
     * @semantics  Destructor of the class Rna6a, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16_alt();

    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

protected:
    virtual void initMatrixIndex();

private:
    static RNA16_alt prototype;

};

#endif //RNA16_alt_H

