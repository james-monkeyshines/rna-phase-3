#ifndef RNA7C_H
#define RNA7C_H

#include "RNA7A.h"

class RNA7C : public RNA7A {
protected:

    /** ***********************************************************************
     * RNA7D
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA7C substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA7C( const string & registrationName );

    /** ************************************************************************
     * RNA7D
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA7C substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA7C( ParametersSet & parameters );

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
     * @semantics  Destructor of the class RNA7C, free dynamic variables
     ************************************************************************ */
    virtual ~RNA7C();

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
    static RNA7C prototype;

};

#endif //RNA7C_H

