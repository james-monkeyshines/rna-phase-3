#ifndef RNA16_bps_alt_H
#define RNA16_bps_alt_H

#include "RNA16_bps.h"

class RNA16_bps_alt : public RNA16_bps {
protected:

    /** ***********************************************************************
     * RNA16_bps_alt
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA16_bps_alt Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16_bps_alt( const string & registrationName );

    /** ************************************************************************
     * RNA16_bps_alt
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16_bps_alt substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16_bps_alt( ParametersSet & parameters );

public:

    /* ***********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA16_bps_alt( parameters );
    }

    /** ************************************************************************
     * ~RNA16_bps_alt
     * @semantics  Destructor of the class Rna6a, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16_bps_alt();

    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

protected:
    virtual void initMatrixIndex();

private:
    static RNA16_bps_alt prototype;

};

#endif //RNA16_bps_alt_H

