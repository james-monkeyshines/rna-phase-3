#ifndef TWOSTATE_H
#define TWOSTATE_H

/** ****************************************************************************
 * The General two-state Time Reversible Model 
 *
 *       0          1    
 * A|    -    |   r0f1  |
 * C|   r0f0  |    -    |
 *
 * r0 = 1
 *
 * One can use a discrete gamma distribution of substitution rate over sites
 * and include a proportion of invariant sites.
 *
 * If we do not use a gamma distribution or invariant sites, we use a factor mr 
 * with that matrix in order to have :
 * sum(matrix * [f0,f1]' * mr)=average rate of substitutions per unit of
 * evolutionary time.
 * with invariant sites :
 * sum(matrix * [f0,f1]' * mr) * (1.0 - proportionInvariantSites)= average rate of substitutions
 * with a gamma dist :
 * sum(matrix * [f0,f1]' * mr) = average rate of substitutions for the
 * category (computed with the average substitution rate and alpha and the
 * proportion of invariant sites if an invariant category is used)
 *                       
 * this model is used for RY-coding, (AU/GU/GC)-(UA/UG/UC) coding and
 * binary coding
 *
 * f0 = frequencies[0]
 * f1 = frequencies[1]
 *
 **************************************************************************** */


#include <math.h>
#include <vector>

#include "Models/TwoStateModel.h"
#include "Util/matrixmath.h"
#include "Util/array2D.h"

class ParametersSet;

class Perturbator;

class SequenceTable;

class TWOSTATE : public TwoStateModel {
private:
    /** .**********************************************************************
     * TWOSTATE
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the General Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    TWOSTATE( const string & registrationName );

protected:
    /** .***********************************************************************
     * TWOSTATE
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the TWOSTATE substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    TWOSTATE( ParametersSet & parameters );

public:
    /* .**********************************************************************
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
     * ~TWOSTATE
     * @semantics  Destructor of the class REV, free dynamic variables
     ************************************************************************ */
    virtual ~TWOSTATE();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

protected:
    virtual void setEigenMatrix();
 
private:
    static TWOSTATE prototype;
};
#endif //TWOSTATE_H




