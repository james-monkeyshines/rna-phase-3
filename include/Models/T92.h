#ifndef T92_H
#define T92_H

/** ****************************************************************************
 * The Galtier Gouy model 
 *
 *            A              C            G              T
 * A|         -       |   r0fCG/2  |   r5fCG/2  |   r0(1-fCG)/2  |
 * C|   r0(1-fCG)/2   |      -     |   r0fCG/2  |   r5(1-fCG)/2  |
 * G|   r5(1-fCG)/2   |   r0fCG/2  |      -     |   r0(1-fCG)/2  |
 * T|   r0(1-fCG)/2   |   r5fCG/2  |   r0fCG/2  |        -       |
 *
 * r0 = transversion ratio
 * r5 = 1 (transition used as a reference
 *
 * One can use a discrete gamma distribution of substitution rate over sites
 * and include a proportion of invariant sites.
 *
 * fA = frequencies[0]/2
 * fC = frequencies[1]/2
 * fG = frequencies[1]/2
 * fU = frequencies[0]/2
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

class T92 : public DnaModel {
private:
    /** .**********************************************************************
     * T92
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the General Time Reversible substitution
     *             model class, this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    T92( const string & registrationName );

protected:
    /** .***********************************************************************
     * T92
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the GTR substitution model object :
     *             a model must be constructed via the model factory
     ************************************************************************ */
    T92( ParametersSet & parameters );

    void initMatrixIndex();

public:
        
    /* .**********************************************************************
    * clone
    * @input      a ParametersSet object to initialise the new model
    * @pre        the set of parameters is consistent (see documentation)
    * @semantics  function used by the ModelFactory to clone the prototype and
    *             create a real model.
    ************************************************************************ */
    virtual Model * clone( ParametersSet & parameters ) const {
        return new T92( parameters );
    }

    /** ************************************************************************
     * ~T92
     * @semantics  Destructor of the class T92, free dynamic variables
     ************************************************************************ */
    virtual ~T92();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;
    
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
  
protected:
    virtual void setEigenMatrix();
 
private:
    static T92 prototype;
};
#endif //T92_H




