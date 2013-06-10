#ifndef RNA7D_H
#define RNA7D_H

#include "RNA7A.h"

class RNA7D : public RNA7A {
protected:

    /** ***********************************************************************
     * RNA7D
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA7D substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA7D( const string & registrationName );

    /** ************************************************************************
     * RNA7D
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA7D substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA7D( ParametersSet & parameters );

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
     * ~RNA7D
     * @semantics  Destructor of the class RNA7D, free dynamic variables
     ************************************************************************ */
    virtual ~RNA7D();

    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

protected:
    virtual void initMatrixIndex();
	virtual bool setParameterNumbers();
	virtual vector < double > condenseRates( array2D < double > & rates16x16 );

private:
    static RNA7D prototype;

};

#endif //RNA7D_H

