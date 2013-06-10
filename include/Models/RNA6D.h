#ifndef RNA6D_H
#define RNA6D_H

#include <math.h>

#include "RNA6E.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"



class RNA6D : public RNA6E {
private:
    /** ***********************************************************************
     * RNA6D
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA6D substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA6D( const string & registrationName );

protected:
    /** ************************************************************************
     * RNA6D
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA6D substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA6D( ParametersSet & parameters );

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
     * ~RNA6D
     * @semantics  Destructor of the class RNA6D, free dynamic variables
     ************************************************************************ */
    virtual ~RNA6D();


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
        unsigned int rateCategory, unsigned int symbolCategory = 0 ) const;

    /** ************************************************************************
     * getFrequencyState
     * @input       an index of frequency
     * @return      a string
     * @semantics   return a name to give to the frequency
     ************************************************************************ */
    virtual string getFrequencyState( unsigned int freqState, unsigned int = 0 ) const;

private:
    static RNA6D prototype;

};

#endif //RNA6D_H




