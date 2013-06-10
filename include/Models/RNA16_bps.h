#ifndef RNA16_bps_H
#define RNA16_bps_H

#include "RNA16.h"

class RNA16_bps : public RNA16 {
protected:

    /** ***********************************************************************
     * RNA16_bps
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA16_bps Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16_bps( const string & registrationName );

    /** ************************************************************************
     * RNA16_bps
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16_bps substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16_bps( ParametersSet & parameters );

public:

    /* ***********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA16_bps( parameters );
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
     * ~RNA16_bps
     * @semantics  Destructor of the class Rna6a, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16_bps();

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

private:
    static RNA16_bps prototype;

};

#endif //RNA16_bps_H

