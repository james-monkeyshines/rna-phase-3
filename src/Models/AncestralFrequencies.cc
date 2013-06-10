#include "Models/AncestralFrequencies.h"

#include "Models/Perturbator.h"
#include "Models/PerturbatorGaussParameter.h"
#include "Models/Model.h"
#include "Util/ParametersSet.h"

AncestralFrequencies::AncestralFrequencies(ParametersSet& parameters,
             unsigned int numberGammaCategories, unsigned int invariant ){
    initPrimitive( parameters, numberGammaCategories, invariant);
}


void AncestralFrequencies::getModelParameters( ParametersSet & parameters, const Model* model, unsigned int symbolCategory ) const{
    ParametersSet* currentCategory;
    char tmpString[40];
    char dblString[40];

    if (numberFrequenciesSets>=1){
        if (numberFrequenciesSets>=2){
            sprintf(tmpString,"%d",numberFrequenciesSets);
            parameters("FREQUENCIES")["Number of frequencies sets"] = tmpString;
            if (gammaCatAdjust & LS){
                sprintf(dblString, "%.8f", *sF);
                parameters("FREQUENCIES")["Variation parameter s"] = dblString;
            }
            unsigned int increment = (gammaCatAdjust&L) ? numberFrequenciesSets - 1 : 1;
            for ( unsigned int i = 0; i < numberFrequenciesSets; i += increment ){
                sprintf(tmpString,"FREQUENCIESCATEGORY%d",i+1);
                currentCategory = &parameters("FREQUENCIES")(tmpString);
                for (unsigned int j = 0; j < numberFrequencies; ++j){
                    sprintf(dblString, "%.8f",frequencies[i][j]);
                    (*currentCategory)["F("+model->getState(j,symbolCategory)+")"] = dblString;
                }
            }
            if (gammaCatAdjust & LX){
                double ratio = 1.0;
                for ( int i = numberFrequenciesXSets-1; i >= 0; --i ){
                    sprintf(tmpString,"FREQUENCIESCATEGORYX%d",i+1);
                    currentCategory = &parameters("FREQUENCIES")(tmpString);
                    ratio = ratio * (*xF)[i];
                    sprintf(dblString, "%.8f", ratio);
                    (*currentCategory)["Variation parameter x"] = dblString;
                    for (unsigned int j = 0; j < numberFrequencies; ++j){
                        sprintf(dblString, "%.8f",frequenciesX[i][j]);
                        (*currentCategory)["F("+model->getState(j,symbolCategory)+")"] = dblString;
                    }
                }
            }
            //store the matching between rate category and set of frequencies
            //in user defined matching case
            if (!gammaCatAdjust&L){
                for ( unsigned int i = 0; i < numberFrequenciesSets; ++i ){
                    sprintf(tmpString,"RATESCATEGORY%d",i+1);
                    currentCategory = &parameters("RATESCATEGORIES")(tmpString);
                    sprintf(dblString,"%d",frequenciesCategory[i]+1);
                   (*currentCategory)["Frequencies set for the category"] = dblString;
                }
            }
        }
        else{
            for (unsigned int j = 0; j < numberFrequencies; ++j){
                sprintf(dblString, "%.8f",frequencies[0][j]);
                parameters("FREQUENCIES")["F("+model->getState(j,symbolCategory)+")"] = dblString;
            }
        }
    }

}
