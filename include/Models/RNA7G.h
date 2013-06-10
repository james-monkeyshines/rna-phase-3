#ifndef RNA7G_H
#define RNA7G_H

#include "RNA7E.h"

class RNA7G : public RNA7E {
private:
    /** ***********************************************************************
     * RNA7G
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA 7G substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA7G( const string & registrationName );

protected:
    /** ************************************************************************
     * RNA7G
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA7E substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA7G( ParametersSet & parameters );

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
     * ~RNA7G
     * @semantics  Destructor of the class RNA7G, free dynamic variables
     ************************************************************************ */
    virtual ~RNA7G();


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
    virtual string getFrequencyState( unsigned int freqState, unsigned int = 0  ) const;

protected:
	virtual bool setParameterNumbers();

private:
    static RNA7G prototype;

};

#endif //RNA7G_H

