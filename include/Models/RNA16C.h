#ifndef RNA16C_H
#define RNA16C_H

#include <math.h>

#include "RNA16A.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"


class RNA16C : public RNA16A {
protected:

    /** ***********************************************************************
     * RNA16C
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA16C Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16C( const string & registrationName );

    /** ************************************************************************
     * RNA16C
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16C substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16C( ParametersSet & parameters );

public:

    /* ***********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA16C( parameters );
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
     * ~RNA16C
     * @semantics  Destructor of the class RNA16C, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16C();


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
    virtual string getFrequencyState( unsigned int freqState, unsigned int = 0 ) const;

            
private:
    static RNA16C prototype;

};

#endif //RNA16C_H




