#ifndef RNA16D_H
#define RNA16D_H

#include <math.h>

#include "RNA16REVEQ.h"

#include "Util/matrixmath.h"
#include "Util/array2D.h"


class RNA16D : public RNA16REVEQ {
protected:

    /** ***********************************************************************
     * RNA16D
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA16D Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA16D( const string & registrationName );

    /** ************************************************************************
     * RNA16D
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA16D substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA16D( ParametersSet & parameters );

public:

    /** ************************************************************************
     * getExchangeability
     * @input  residue1, the index of a state of the model
     * @input  residue2, the index of a state of the model
     * @input  gammaCategory, a discrete gamma category
     * @return  The residue1<->residue2 exchangeability parameter for the given
     *          gamma category redefined for this model.
     ************************************************************************ */
    virtual double getExchangeability( unsigned int residue1,
        unsigned int residue2, unsigned int gammaCategory ) const;
        
        
   /* ***********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new RNA16D( parameters );
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
     * ~RNA16D
     * @semantics  Destructor of the class RNA16D, free dynamic variables
     ************************************************************************ */
    virtual ~RNA16D();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

    /** ************************************************************************
     * getFrequency
     * @intput  residue, the index of a state of the model
     * @intput  rateCategory, the discrete gamma category
     * @intput  symbolCategory, must be 0
     * @return  The frequencies of the given state of the model
     ************************************************************************ */
    virtual double getFrequency( unsigned int residue,
        unsigned int rateCategory, unsigned int = 0 ) const;
    
protected:
    virtual void initMatrixIndex();

    double getInvKappa( unsigned int rateCategory ) const;
            
private:
    static RNA16D prototype;

};

#endif //RNA16D_H



