#ifndef RNA7A_H
#define RNA7A_H

#include <math.h>
#include "Models/RnaModel.h"

class RNA7A : public RnaModel {
protected:
    RNA7A();

    /** ***********************************************************************
     * RNA7A
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the RNA General Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    RNA7A( const string & registrationName );

    /** ************************************************************************
     * RNA7A
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the RNA7A substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    RNA7A( ParametersSet & parameters );

public:
   /* ***********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model. The behaviour of this function is
    *             slightly different compared to all other models:
    *             7-state models can have extra parameters to allow
    *             for mismatch frequency and rate variations among
    *             gamma category. These are constraints on the more general
    *             linear variations implement in MatrixModel
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
     * ~RNA7A
     * @semantics  Destructor of the class RNA7A, free dynamic variables
     ************************************************************************ */
    virtual ~RNA7A();

    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

    /** ************************************************************************
     * getState
     * @input   stateNumber, the index of a state of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getState( unsigned int stateNumber, unsigned int symbolCategory = 0 ) const;

protected:
    virtual void initMatrixIndex();
	virtual bool setParameterNumbers();
    vector < double > mergeFrequencies( vector < double > & freqs16, bool bp_symmetry = false );
	virtual vector < double > condenseRates( array2D < double > & rates16x16 );
	virtual void lnLAdjustment( vector < double > & countSite );
	void setEquivalencyTable();
	
private:
    static RNA7A prototype;
    
};

#endif //RNA7A_H

