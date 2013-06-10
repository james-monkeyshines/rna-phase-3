#ifndef RNA7B_H
#define RNA7B_H

#include "RNA7A.h"

class RNA7B : public RNA7A {
protected:

    /** ***********************************************************************
     * RNA7B
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA 7B substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA7B( const string & registrationName );

    /** ************************************************************************
     * RNA7B
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA7B substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA7B( ParametersSet & parameters );

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
     * ~RNA7B
     * @semantics  Destructor of the class Rna7a, free dynamic variables
     ************************************************************************ */
    virtual ~RNA7B();

    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

    /** ************************************************************************
     * getFrequency
     * @intput     residue, the index of a state of the model
     * @intput     rateCategory, the discrete gamma category
     * @intput     symbolCategory, must be 0
     * @return     The frequencies of the given state of the model
     * @semantics  replace the common method from matrix model where
     *             each state has its own frequency
     ************************************************************************ */
    virtual double getFrequency( unsigned int residue,
        unsigned int rateCategory, unsigned int = 0 ) const;
 
    /** ************************************************************************
     * getFrequencyState
     * @input       an index of frequency
     * @return      a string
     * @semantics   return a name to give to the frequency
     ************************************************************************ */
    virtual string getFrequencyState( unsigned int freqState, unsigned int = 0) const;

protected:
	virtual bool setParameterNumbers();

private:
    static RNA7B prototype;

};

#endif //RNA7B_H

