#ifndef YNH98_H
#define YNH98_H

/** ****************************************************************************
 * A simple codon substitution model (Yang, Nielsen and Hasegawa (1998)),
 * a simplified version of the model of Goldman and Yang (1994)
 * substitution rate matrix :
 *
 * Qij =
 *   0,                if the two codons differ at more than one position
 *   Pi_j,             for synonymous transversion
 *   Kappa*Pi_j,       for synonymous transition
 *   Omega*Pi_j,       for nonsynonymous transversion
 *   Omega*Kappa*Pi_j, for nonsynonymous transition
 **************************************************************************** */




 /*
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
 */



#include <math.h>
#include <vector>

#include "Models/CodonModel.h"
#include "Util/matrixmath.h"
#include "Util/array2D.h"
class ParametersSet;

class Perturbator;

class SequenceTable;

class YNH98 : public CodonModel {
private:
    /** .**********************************************************************
     * YNH98
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the simple Codon Model substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    YNH98( const string & registrationName );

protected:
    /** .***********************************************************************
     * YNH98
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the YNH98 substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    YNH98( ParametersSet & parameters );

public:
    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new YNH98( parameters );
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
     * ~YNH98
     * @semantics  Destructor of the class YNH98, free dynamic variables
     ************************************************************************ */
    virtual ~YNH98();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;

    /** ************************************************************************
     * getExchangeability
     * @input  residue1, the index of a state of the model
     * @input  residue2, the index of a state of the model
     * @input  gammaCategory, a discrete gamma category
     * @return  The residue1<->residue2 exchangeability parameter for the given
     *          gamma category
     ************************************************************************ */
    virtual double getExchangeability( unsigned int residue1,
        unsigned int residue2, unsigned int gammaCategory ) const;


private:
    static YNH98 prototype;

    /***************************************************************************
    * defModelMatrix(i,j) = -1 if i=j (diagonal)
    * defModelMatrix(i,j) = 0  if the subsitution rate between i and j is ' 0 '
    * defModelMatrix(i,j) = 1  if the subs rate between i and j is ' PI_j '
    * defModelMatrix(i,j) = 2  if ... is ' omega * PI_j '
    * defModelMatrix(i,j) = 3  if ... is ' kappa * PI_j '
    * defModelMatrix(i,j) = 4  if ... is ' omega * kappa * PI_j '
    *
    * PI_j is the frequence of the codon j
    * omega is the ratio of nonsynonymous/synonymous subsitution rates
    * kappa is the ratio of transition/transversion
    ***************************************************************************/
    array2D <int> defModelMatrix;

    /***************************************************************************
     * initDefModelMatrix()
     * @semantics Initialize the matrix which defines the model
    ***************************************************************************/
    void initDefModelMatrix();

    /***************************************************************************
     * isTransitionOrTransversion
     * @input  c1, a character to compare (U, C, G or T)
     * @input  c2, a character to compare (U, C, G or T)
     * @return  1 if the substitution from c1 to c2 is a transition,
     *         -1 if the substitution from c1 to c2 is a transversion,
     *          0 if c1=c2
    ***************************************************************************/
    int isTransitionOrTransversion(char c1, char c2);
};
#endif //YNH98_H

