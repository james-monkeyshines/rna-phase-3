#ifndef TN93_H
#define TN93_H

/** ****************************************************************************
 * The TN93 Model for DNA substitution :
 * substitution rate matrix :
 *
 *       A          C          G          T
 * A|    -     |  r0fC  |  r5fG  |  r0fT  |
 * C|  r0fC  |    -     |  r0fG  |  r4fT  |
 * G|  r5fA  |  r0fC  |    -     |  r0fT  |
 * T|  r0fA  |  r4fC  |  r0fG |    -     |
 *
 * r0 = rateRatios[0] : transversion
 * r4 = rateRatios[1] : 2nd transition
 * r5 = 1 (1st transversion, reference)
 *
 * One can use a discrete gamma distribution of substitution rate over sites
 * and include a proportion of invariant sites.
 *
 * If we do not use a gamma distribution or invariant sites, we use a factor mr 
 * with that matrix in order to have :
 * sum(matrix * [fA,fC,fG,fT]' * mr)=average rate of substitutions per unit of
 * evolutionary time.
 * with invariant sites :
 * sum(matrix * [fA,fC,fG,fT]' * mr) * (1.0 - proportionInvariantSites)= average rate of substitutions
 * with a gamma dist :
 * sum(matrix * [fA,fC,fG,fT]' * mr) = average rate of substitutions for the
 * category (computed with the average substitution rate and alpha and the
 * proportion of invariant sites if an invariant category is used)
 *                       
 * "A" = 0, "C" = 1, "G" = 2 and "U" = 3
 *
 * fA = frequencies[0]
 * fC = frequencies[1]
 * fG = frequencies[2]
 * fU = frequencies[3]
 *
 **************************************************************************** */


#include <math.h>
#include <vector>

#include "Models/DnaModel.h"
#include "Util/matrixmath.h"
#include "Util/array2D.h"
class ParametersSet;

class Perturbator;

class SequenceTable;

class TN93 : public DnaModel {
private:
    /** .**********************************************************************
     * TN93
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the Timura Nei 93 substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    TN93( const string & registrationName );

protected:
    /** .***********************************************************************
     * TN93
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the TN93 substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    TN93( ParametersSet & parameters );

public:
    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new TN93( parameters );
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
     * ~TN93
     * @semantics  Destructor of the class TN93, free dynamic variables
     ************************************************************************ */
    virtual ~TN93();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

protected:
    virtual void setEigenMatrix();
    void initMatrixIndex();
 
private:
    static TN93 prototype;
};
#endif //TN93_H

