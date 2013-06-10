#include "Models/F81.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name F81
F81 F81::prototype( "F81" );

F81::F81( const string & registrationName ) : DnaModel( registrationName ){
//the private constructor is called by the prototype only
}

F81::F81( ParametersSet & parameters ) : DnaModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void F81::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 0;
    numberFrequencies = 4;
    
    //initialise the matrixIndex for correspondance between a matrix element
    //and an index in the rate ratios vector
    //initMatrixIndex();
    
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
        ratesRatios->initialisation( 0 );
        if ( invariantCategory ) {
            retrieveEmpiricalRates( sequenceTable, modelId,
                         & proportionInvariantSites );
        }
        else {
            //not useful
        }
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();    
}
  
F81::~F81() {
}

string F81::getName( void ) const {
    return ( "F81" + MatrixModel::getName() );
}

void F81::setEigenMatrix() {
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
    category++ ) {
        // fill the rateMatrix (cf F81.h)
        int freqCat = frequencies->frequenciesCat(category+invariantCategory);
        rateMatrix[category] ( 0, 1 ) = substitutionRate[category] * (*frequencies)(freqCat)[1]; // A->C
        rateMatrix[category] ( 0, 2 ) = substitutionRate[category] * (*frequencies)(freqCat)[2]; // A->G
        rateMatrix[category] ( 0, 3 ) = substitutionRate[category] * (*frequencies)(freqCat)[3]; // A->U
        rateMatrix[category] ( 0, 0 ) = -( ( rateMatrix[category] ) ( 0, 1 ) + ( rateMatrix[category] ) ( 0, 2 ) + ( rateMatrix[category] ) ( 0, 3 ) );

        rateMatrix[category] ( 1, 0 ) = substitutionRate[category] * (*frequencies)(freqCat)[0]; // C->A
        rateMatrix[category] ( 1, 2 ) = substitutionRate[category] * (*frequencies)(freqCat)[2]; // C->G
        rateMatrix[category] ( 1, 3 ) = substitutionRate[category] * (*frequencies)(freqCat)[3]; // C->U
        rateMatrix[category] ( 1, 1 ) = -( ( rateMatrix[category] ) ( 1, 0 ) + ( rateMatrix[category] ) ( 1, 2 ) + ( rateMatrix[category] ) ( 1, 3 ) );

        rateMatrix[category] ( 2, 0 ) = substitutionRate[category] * (*frequencies)(freqCat)[0]; // G->A
        rateMatrix[category] ( 2, 1 ) = substitutionRate[category] * (*frequencies)(freqCat)[1]; // G->C
        rateMatrix[category] ( 2, 3 ) = substitutionRate[category] * (*frequencies)(freqCat)[3]; // G->U
        rateMatrix[category] ( 2, 2 ) = -( ( rateMatrix[category] ) ( 2, 0 ) + ( rateMatrix[category] ) ( 2, 1 ) + ( rateMatrix[category] ) ( 2, 3 ) );

        rateMatrix[category] ( 3, 0 ) = substitutionRate[category] * (*frequencies)(freqCat)[0]; // U->A
        rateMatrix[category] ( 3, 1 ) = substitutionRate[category] * (*frequencies)(freqCat)[1]; // U->C
        rateMatrix[category] ( 3, 2 ) = substitutionRate[category] * (*frequencies)(freqCat)[2]; // U->G
        rateMatrix[category] ( 3, 3 ) = -( ( rateMatrix[category] ) ( 3, 0 ) + ( rateMatrix[category] ) ( 3, 1 ) + ( rateMatrix[category] ) ( 3, 2 ) );

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
        eigenValues[category] [1] = -substitutionRate[category];
        eigenValues[category] [2] = -substitutionRate[category];
        eigenValues[category] [3] = -substitutionRate[category];


        // Update the eigenmatrix by hand
        eigenMatrix[category] ( 0, 0 ) = (*frequencies)(freqCat)[0];
        eigenMatrix[category] ( 0, 1 ) = (*frequencies)(freqCat)[1];
        eigenMatrix[category] ( 0, 2 ) = (*frequencies)(freqCat)[2];
        eigenMatrix[category] ( 0, 3 ) = (*frequencies)(freqCat)[3];

        eigenMatrix[category] ( 1, 0 ) = 0.0;
        eigenMatrix[category] ( 1, 1 ) = -1.0;
        eigenMatrix[category] ( 1, 2 ) = 0.0;
        eigenMatrix[category] ( 1, 3 ) = 1.0;

        eigenMatrix[category] ( 2, 0 ) = -1.0;
        eigenMatrix[category] ( 2, 1 ) = 0.0;
        eigenMatrix[category] ( 2, 2 ) = 1.0;
        eigenMatrix[category] ( 2, 3 ) = 0.0;

        eigenMatrix[category] ( 3, 0 ) = -( (*frequencies)(freqCat)[0] * pyrimidineFreq ) / purineFreq;
        eigenMatrix[category] ( 3, 1 ) = (*frequencies)(freqCat)[1];
        eigenMatrix[category] ( 3, 2 ) = -( (*frequencies)(freqCat)[2] * pyrimidineFreq ) / purineFreq;
        eigenMatrix[category] ( 3, 3 ) = (*frequencies)(freqCat)[3];

        // Inverse of the eigenmatrix
        ieigenMatrix[category] ( 0, 0 ) = 1.0;
        ieigenMatrix[category] ( 0, 1 ) = 0.0;
        ieigenMatrix[category] ( 0, 2 ) = -(*frequencies)(freqCat)[2] / purineFreq;
        ieigenMatrix[category] ( 0, 3 ) = -1.0;

        ieigenMatrix[category] ( 1, 0 ) = 1.0;
        ieigenMatrix[category] ( 1, 1 ) = -(*frequencies)(freqCat)[3] / pyrimidineFreq;
        ieigenMatrix[category] ( 1, 2 ) = 0.0;
        ieigenMatrix[category] ( 1, 3 ) = purineFreq / pyrimidineFreq;

        ieigenMatrix[category] ( 2, 0 ) = 1.0;
        ieigenMatrix[category] ( 2, 1 ) = 0.0;
        ieigenMatrix[category] ( 2, 2 ) = (*frequencies)(freqCat)[0] / purineFreq;
        ieigenMatrix[category] ( 2, 3 ) = -1.0;

        ieigenMatrix[category] ( 3, 0 ) = 1.0;
        ieigenMatrix[category] ( 3, 1 ) = (*frequencies)(freqCat)[1] / pyrimidineFreq;
        ieigenMatrix[category] ( 3, 2 ) = 0.0;
        ieigenMatrix[category] ( 3, 3 ) = purineFreq / pyrimidineFreq;
    } //for each category (except invariant)
}

