#ifndef RNA7E_H
#define RNA7E_H

#include "RNA7D.h"

class RNA7E : public RNA7D {
protected:

    /** ***********************************************************************
     * RNA7E
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA7E substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA7E( const string & registrationName );

    /** ************************************************************************
     * RNA7E
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA7E substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA7E( ParametersSet & parameters );

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
     * ~RNA7E
     * @semantics  Destructor of the class RNA7D, free dynamic variables
     ************************************************************************ */
    virtual ~RNA7E();

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
    static RNA7E prototype;

};

#endif //RNA7E_H

