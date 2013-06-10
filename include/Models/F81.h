#ifndef F81_H
#define F81_H

/** ****************************************************************************
 * The F81 Model for DNA substitution :
 * substitution rate matrix :
 *
 *      A      C      G      T
 * A|  -   |  fC  |  fG  |  fT  |
 * C|  fA  |  -   |  fG  |  fT  |
 * G|  fA  |  fC  |  -   |  fT  |
 * T|  fA  |  fC  |  fG  |   -  |
 *
 * r5 = 1
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

class F81 : public DnaModel {
private:
    /** .**********************************************************************
     * F81
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the HKY 85 substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    F81( const string & registrationName );

protected:
    /** .***********************************************************************
     * F81
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the F81 substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    F81( ParametersSet & parameters );

public:
    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new F81( parameters );
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
     * ~F81
     * @semantics  Destructor of the class F81, free dynamic variables
     ************************************************************************ */
    virtual ~F81();


    /** ************************************************************************
     * getExchangeability
     * @input  residue1, the index of a state of the model
     * @input  residue2, the index of a state of the model
     * @input  gammaCategory, a discrete gamma category
     * @return  The residue1<->residue2 exchangeability parameter for the given
     *          gamma category
     * @semantics  A JC69 model does not define matrixIndex
     *             (just one rate) and the base method has to be redefined
     *             this method is probably not used
     ************************************************************************ */
    virtual double getExchangeability( unsigned int residue1,
				unsigned int residue2, unsigned int gammaCategory ) const {
    	return 1.0;
    }

        
    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;


protected:
    virtual void setEigenMatrix();
 
private:
    static F81 prototype;
};
#endif //F81_H

