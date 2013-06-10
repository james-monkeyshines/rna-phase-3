#ifndef RNA6C_H
#define RNA6C_H

#include <math.h>

#include "RNA6B.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"


// The general time reversible model for nucleotide substitution with
// using the discrete gamma model with some sites invariant


class RNA6C : public RNA6B {
private:
    /** .**********************************************************************
     * RNA6C
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA6C substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA6C( const string & registrationName );

protected:
    /** .***********************************************************************
     * RNA6C
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA6C substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA6C( ParametersSet & parameters );

public:

    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA6C( parameters );
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
     * ~RNA6C
     * @semantics  Destructor of the class RNA6C, free dynamic variables
     ************************************************************************ */
    virtual ~RNA6C();


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
    static RNA6C prototype;

};

#endif //RNA6C_H




