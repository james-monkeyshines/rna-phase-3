#include <assert.h>

#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Models/Perturbator.h"
#include "Models/TN93.h"

#include "Sequence/SequenceTable.h"

#include "Util/statlib.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"


#define SQR(n) (n*n)

using namespace std;

//register the model to the model factory with the name TN93
TN93 TN93::prototype( "TN93" );

TN93::TN93( const string & registrationName ) : DnaModel( registrationName ){
//the private constructor is called by the prototype only
}

TN93::TN93( ParametersSet & parameters ) : DnaModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void TN93::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 2;
    numberFrequencies = 4;

    //initialise the matrixIndex for correspondance between a matrix element
    //and an index in the rate ratios vector
    initMatrixIndex();
        
    //call the basic initialisation method
    DnaModel::initialisation( sequenceTable, modelId );
    
    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();
    
    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if (sequenceTable){
        //estimation of the model frequencies from the empirical frequencies
        vector<double> freq4 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        frequencies->initialisation( freq4 );

        //estimation of all the rate ratios (and proportion of invariant sites
        //if relevant) from the empirical sequences
        vector < double > rates5;
        if ( invariantCategory ) {
            rates5 = retrieveEmpiricalRates( sequenceTable, modelId,
            & proportionInvariantSites );
        }
        else {
            rates5 = retrieveEmpiricalRates( sequenceTable, modelId );
        }
        //the 2 transversion rates are the only free parameter
        vector < double > rates2;
        rates2.push_back( (rates5[0]+rates5[1]+rates5[2]+rates5[4])/4.0 );
        rates2.push_back( rates5[3] );
        ratesRatios->initialisation( rates2 );        
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();    
}
  
TN93::~TN93() {
}

void TN93::initMatrixIndex() {
    matrixIndex.resize( 4, 4 );
    for ( int i = 0 ; i < 4 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
    //reference
    matrixIndex(0, 2) = matrixIndex(2, 0) = -2;
    //transversion
    matrixIndex(0, 1) = matrixIndex(1, 0) = 0;
    matrixIndex(0, 3) = matrixIndex(3, 0) = 0;
    matrixIndex(1, 2) = matrixIndex(2, 1) = 0;
    matrixIndex(2, 3) = matrixIndex(3, 2) = 0;
    //second transition
    matrixIndex(1, 3) = matrixIndex(3, 1) = 1;
}

string TN93::getName( void ) const {
    return ( "TN93" + MatrixModel::getName() );
}


void TN93::setEigenMatrix() {
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
    category++ ) {
        // fill the rateMatrix (cf Hky.h)
        int freqCat = frequencies->frequenciesCat(category+invariantCategory);
        int rateCat = ratesRatios->ratesRatiosCat(category);
        rateMatrix[category] ( 0, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [1] *
        (*ratesRatios)(rateCat) [0]; // A->C
        rateMatrix[category] ( 0, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [2];
        rateMatrix[category] ( 0, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [3] *
        (*ratesRatios)(rateCat) [0]; // A->U
        rateMatrix[category] ( 0, 0 ) =
        -( ( rateMatrix[category] ) ( 0, 1 ) +
        ( rateMatrix[category] ) ( 0, 2 ) +
        ( rateMatrix[category] ) ( 0, 3 ) );
        rateMatrix[category] ( 1, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [0] *
        (*ratesRatios)(rateCat) [0]; // C->A
        rateMatrix[category] ( 1, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [2] *
        (*ratesRatios)(rateCat) [0]; // C->G
        rateMatrix[category] ( 1, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [3] *
        (*ratesRatios)(rateCat) [1]; // C->U
        rateMatrix[category] ( 1, 1 ) =
        -( ( rateMatrix[category] ) ( 1, 0 ) +
        ( rateMatrix[category] ) ( 1, 2 ) +
        ( rateMatrix[category] ) ( 1, 3 ) );
        rateMatrix[category] ( 2, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [0]; // G->A
        rateMatrix[category] ( 2, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [1] *
        (*ratesRatios)(rateCat) [0]; // G->C
        rateMatrix[category] ( 2, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [3] *
        (*ratesRatios)(rateCat) [0]; // G->U
        rateMatrix[category] ( 2, 2 ) =
        -( ( rateMatrix[category] ) ( 2, 0 ) +
        ( rateMatrix[category] ) ( 2, 1 ) +
        ( rateMatrix[category] ) ( 2, 3 ) );
        rateMatrix[category] ( 3, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [0] *
        (*ratesRatios)(rateCat) [0]; // U->A
        rateMatrix[category] ( 3, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [1] *
        (*ratesRatios)(rateCat) [1]; // U->C
        rateMatrix[category] ( 3, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [2] *
        (*ratesRatios)(rateCat) [0]; // U->G
        rateMatrix[category] ( 3, 3 ) =
        -( ( rateMatrix[category] ) ( 3, 0 ) +
        ( rateMatrix[category] ) ( 3, 1 ) +
        ( rateMatrix[category] ) ( 3, 2 ) );

        for ( int i = 0; i < 4; ++i ) {
            eigenValues[category] [i] = 0.0;
            for ( int j = 0; j < 4; ++j ) {
                eigenMatrix[category] ( i, j ) = 0.0;
                ieigenMatrix[category] ( i, j ) = 0.0;
            }
        }
                
        double purineFreq = (*frequencies)(freqCat)[0] + (*frequencies)(freqCat)[2];
        double pyrimidineFreq = (*frequencies)(freqCat)[1] + (*frequencies)(freqCat)[3];

        eigenValues[category] [0] = 0.0;
        eigenValues[category] [1] = - (*ratesRatios)(rateCat) [0] * substitutionRate[category];
        eigenValues[category] [2] = - substitutionRate[category] *
        ( (*ratesRatios)(rateCat) [0] * pyrimidineFreq + purineFreq) ;
        eigenValues[category] [3] =  - substitutionRate[category] *
        ( (*ratesRatios)(rateCat) [0] * purineFreq +
          (*ratesRatios)(rateCat) [1] * pyrimidineFreq );


        // Update the eigenmatrix by hand
        eigenMatrix[category] ( 0, 0 ) = (*frequencies)(freqCat)[0];
        eigenMatrix[category] ( 0, 1 ) = (*frequencies)(freqCat)[1];
        eigenMatrix[category] ( 0, 2 ) = (*frequencies)(freqCat)[2];
        eigenMatrix[category] ( 0, 3 ) = (*frequencies)(freqCat)[3];

        eigenMatrix[category] ( 1, 0 ) = -( (*frequencies)(freqCat)[0] * pyrimidineFreq ) / purineFreq;
        eigenMatrix[category] ( 1, 1 ) = (*frequencies)(freqCat)[1];
        eigenMatrix[category] ( 1, 2 ) = -( (*frequencies)(freqCat)[2] * pyrimidineFreq ) / purineFreq;
        eigenMatrix[category] ( 1, 3 ) = (*frequencies)(freqCat)[3];

        eigenMatrix[category] ( 2, 0 ) = -1.0;
        eigenMatrix[category] ( 2, 1 ) = 0.0;
        eigenMatrix[category] ( 2, 2 ) = 1.0;
        eigenMatrix[category] ( 2, 3 ) = 0.0;

        eigenMatrix[category] ( 3, 0 ) = 0.0;
        eigenMatrix[category] ( 3, 1 ) = -1.0;
        eigenMatrix[category] ( 3, 2 ) = 0.0;
        eigenMatrix[category] ( 3, 3 ) = 1.0;

        // Inverse of the eigenmatrix
        ieigenMatrix[category] ( 0, 0 ) = 1.0;
        ieigenMatrix[category] ( 0, 1 ) = -1.0;
        ieigenMatrix[category] ( 0, 2 ) = -(*frequencies)(freqCat)[2] / purineFreq;
        ieigenMatrix[category] ( 0, 3 ) = 0.0;

        ieigenMatrix[category] ( 1, 0 ) = 1.0;
        ieigenMatrix[category] ( 1, 1 ) = purineFreq / pyrimidineFreq;
        ieigenMatrix[category] ( 1, 2 ) = 0.0;
        ieigenMatrix[category] ( 1, 3 ) = -(*frequencies)(freqCat)[3] / pyrimidineFreq;

        ieigenMatrix[category] ( 2, 0 ) = 1.0;
        ieigenMatrix[category] ( 2, 1 ) = -1.0;
        ieigenMatrix[category] ( 2, 2 ) = (*frequencies)(freqCat)[0] / purineFreq;
        ieigenMatrix[category] ( 2, 3 ) = 0.0;

        ieigenMatrix[category] ( 3, 0 ) = 1.0;
        ieigenMatrix[category] ( 3, 1 ) = purineFreq / pyrimidineFreq;
        ieigenMatrix[category] ( 3, 2 ) = 0.0;
        ieigenMatrix[category] ( 3, 3 ) = (*frequencies)(freqCat)[1] / pyrimidineFreq;
    } //for each category (except invariant)
}

