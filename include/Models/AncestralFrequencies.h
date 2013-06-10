#ifndef ANCESTRALFREQUENCIES_H
#define ANCESTRALFREQUENCIES_H

#include "Models/Frequencies.h"

class AncestralFrequencies: public Frequencies{
public:
    AncestralFrequencies(ParametersSet& parameters,
        unsigned int numberGammaCategories, unsigned int invariant );

//    void initialiseMCMC( ParametersSet & parameters,
//            Perturbator* perturbator );

    /** ************************************************************************
     * getModelParameters
     * @semantics   ancestral frequency parameters depends on the number of
     *              state, not the number of frequencies
     *              cf Heterogeneous::initialisation for the issue...
     ************************************************************************ */
    virtual void getModelParameters( ParametersSet & parameters, const Model* model, unsigned int symbolCategory ) const;

};
#endif
