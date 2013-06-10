#ifndef K80_H
#define K80_H

/** ****************************************************************************
 * The K80 Model for DNA substitution :
 * substitution rate matrix :
 *
 *       A          C          G          T
 * A|    -    |  .25r0  |  .25r5  |  r0.25  |
 * C|  .25r0  |    -    |  .25r0  |  .25r5  |
 * G|  .25r5  |  .25r0  |    -    |  .25r0  |
 * T|  .25r0  |  .25r5  |  .25r0  |    -    |
 *
 * r0 = rateRatios[0] : transversion
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
 * fA = fC = fG = fU = .25
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

class K80 : public DnaModel {
private:
    /** .**********************************************************************
     * K80
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the Kimura 2 parameters substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    K80( const string & registrationName );

protected:
    /** .***********************************************************************
     * K80
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the K2P substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    K80( ParametersSet & parameters );

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
     * ~K80
     * @semantics  Destructor of the class K80, free dynamic variables
     ************************************************************************ */
    virtual ~K80();

    
    /** ************************************************************************
     * getFrequency
     * @intput  residue, the index of a state of the model
     * @intput  rateCategory, the discrete gamma category
     * @intput  symbolCategory, must be 0
     * @return  The frequencies of the given state of the model
     ************************************************************************ */
    inline virtual double getFrequency( unsigned int, unsigned int,
                                                     unsigned int = 0 ) const {
        return .25;
    }


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

    /** ************************************************************************
     * updateEigenMatrix
     * @input  none
     * @semantics  compute the multiplier in front of each rate matrix to have
     *             the right averageRate for each of them. Call setEigenMatrix
     *             once it is done.
     ************************************************************************ */
    virtual void updateEigenMatrix();


    /** ************************************************************************
     * printParameters
     * @input      outputStream, the stream used to output the parameters
     * @semantics  Output the parameters names and values of the model in a
     *             readable format.
     *             -frequencies of the state
     *             -gamma shape parameter (if relevant)
     *             -proportion of invariant sites (if relevant)
     *             -substitution rate for each category
     ************************************************************************ */
    virtual void printParameters( ostream & outputStream ) const;


protected:
    virtual void setEigenMatrix();
    void initMatrixIndex();
 
private:
    static K80 prototype;
};
#endif //K80_H

