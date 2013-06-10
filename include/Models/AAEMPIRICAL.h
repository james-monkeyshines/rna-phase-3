#ifndef AAEMPIRICAL_H
#define AAEMPIRICAL_H

/** ****************************************************************************
 * The Empirical Models for amino acid substitution :
 *
 * Dayhoff (1978), Jones (1992), mtREV24, ...
 *
 * The rates are not estimated. They are defined in files.
 * Unless it is a +F model, frequencies should be in the file as well.
 **************************************************************************** */

#include <math.h>
#include <vector>

#include "Models/AaModel.h"
#include "Util/matrixmath.h"
#include "Util/array2D.h"
class ParametersSet;

class Perturbator;

class SequenceTable;

class AAEMPIRICAL : public AaModel {
private:
    /** .**********************************************************************
     * AAEMPIRICAL
     * @input      registrationName, the name used to register to the model
     *             factory
     * @semantics  private constructor of the substitution model class,
     *             this construtor is used only once at launch
     *             time to register the prototype to the ModelFactory
     ************************************************************************ */
    AAEMPIRICAL( const string & registrationName );

protected:
    /** .***********************************************************************
     * AAEMPIRICAL
     * @input      a ParametersSet object to initialise the model
     * @pre        the set of parameters is consistent (see documentation)
     * @semantics  protected constructor of the AAEMPIRICAL substitution model
     *             object: a model must be constructed via the model factory
     ************************************************************************ */
    AAEMPIRICAL( ParametersSet & parameters );

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
     * ~AAEMPIRICAL
     * @semantics  Destructor of the class AAEMPIRICAL, free dynamic variables
     ************************************************************************ */
    virtual ~AAEMPIRICAL();


    /** ************************************************************************
     * getName
     * @return  The name of the model
     ************************************************************************ */
    virtual string getName( void ) const;
    
    
    /** ************************************************************************
     * getFrequency
     * @input  residue, the index of a state of the model
     * @input  residue2, the index of a state of the model
     * @input  gammaCategory, a discrete gamma category
     * @return  The residue1<->residue2 exchangeability parameter for the given
     *          gamma category
     ************************************************************************ */
    virtual double getFrequency( unsigned int residue,
        unsigned int gammaCategory, unsigned int symbolCategory = 0 ) const;

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


protected:
    virtual void setEigenMatrix();
    void readEmpiricalRates();
    void readEmpiricalFrequencies();
    bool plusF;
    string empiricalRatesFileName;
    
private:
    static AAEMPIRICAL prototype;
    ifstream empiricalValuesFile;
    array2D<double> empiricalRates;
    vector<double> empiricalFrequencies;
};
#endif //AAEMPIRICAL_H

