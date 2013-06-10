#include "Models/REV.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name REV
REV REV::prototype( "REV" );

REV::REV( const string & registrationName ) : DnaModel( registrationName ){
//the private constructor is called by the prototype only
}


REV::REV( ParametersSet & parameters ) : DnaModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void REV::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 5;
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
        ratesRatios->initialisation( rates5 );
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();    
}  

REV::~REV() {
}

string REV::getName( void ) const {
    return ( "REV" + MatrixModel::getName() );
}

void REV::initMatrixIndex() {
    matrixIndex.resize( 4, 4 );
    for ( int i = 0 ; i < 4 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
    //reference
    matrixIndex(0, 2) = matrixIndex(2, 0) = -2;
    //transversion
    matrixIndex(0, 1) = matrixIndex(1, 0) = 0;
    matrixIndex(0, 3) = matrixIndex(3, 0) = 1;
    matrixIndex(1, 2) = matrixIndex(2, 1) = 2;
    matrixIndex(2, 3) = matrixIndex(3, 2) = 4;
    //second transition
    matrixIndex(1, 3) = matrixIndex(3, 1) = 3;
}

void REV::setEigenMatrix() {  
    vector< double > junk( 4 );
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
    category++ ) {
        // fill the rateMatrix (cf REV.h)
        int freqCat = frequencies->frequenciesCat(category+invariantCategory);
        int rateCat = ratesRatios->ratesRatiosCat(category);
        rateMatrix[category] ( 0, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1] *
        (*ratesRatios)(rateCat)[0]; // A->C
        rateMatrix[category] ( 0, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[2];
        rateMatrix[category] ( 0, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[3] *
        (*ratesRatios)(rateCat)[1]; // A->U
        rateMatrix[category] ( 0, 0 ) =
        -( ( rateMatrix[category] ) ( 0, 1 ) +
        ( rateMatrix[category] ) ( 0, 2 ) +
        ( rateMatrix[category] ) ( 0, 3 ) );
        rateMatrix[category] ( 1, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0] *
        (*ratesRatios)(rateCat)[0]; // C->A
        rateMatrix[category] ( 1, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[2] *
        (*ratesRatios)(rateCat)[2]; // C->G
        rateMatrix[category] ( 1, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[3] *
        (*ratesRatios)(rateCat)[3]; // C->U
        rateMatrix[category] ( 1, 1 ) =
        -( ( rateMatrix[category] ) ( 1, 0 ) +
        ( rateMatrix[category] ) ( 1, 2 ) +
        ( rateMatrix[category] ) ( 1, 3 ) );
        rateMatrix[category] ( 2, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0]; // G->A
        rateMatrix[category] ( 2, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1] *
        (*ratesRatios)(rateCat)[2]; // G->C
        rateMatrix[category] ( 2, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[3] *
        (*ratesRatios)(rateCat)[4]; // G->U
        rateMatrix[category] ( 2, 2 ) =
        -( ( rateMatrix[category] ) ( 2, 0 ) +
        ( rateMatrix[category] ) ( 2, 1 ) +
        ( rateMatrix[category] ) ( 2, 3 ) );
        rateMatrix[category] ( 3, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0] *
        (*ratesRatios)(rateCat)[1]; // U->A
        rateMatrix[category] ( 3, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1] *
        (*ratesRatios)(rateCat)[3]; // U->C
        rateMatrix[category] ( 3, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[2] *
        (*ratesRatios)(rateCat)[4]; // U->G
        rateMatrix[category] ( 3, 3 ) =
        -( ( rateMatrix[category] ) ( 3, 0 ) +
        ( rateMatrix[category] ) ( 3, 1 ) +
        ( rateMatrix[category] ) ( 3, 2 ) );

        for ( int i = 0; i < 4; ++i ) {
            eigenValues[category] [i] = 0.0;
            junk[i] = 0.0;
            for ( int j = 0; j < 4; ++j ) {
                eigenMatrix[category] ( i, j ) = 0.0;
                ieigenMatrix[category] ( i, j ) = 0.0;
            }
        }
        // Find the eigenValues and eigenVectors of the rate matrix
        int err = LeftEigenSystem( rateMatrix[category], eigenMatrix[category],
                            eigenValues[category], junk );
        if (err){
            cerr << "Error during matrix exponentiation" << endl;
            exit(EXIT_FAILURE);
        }
        // Find the inverse of the matrix whose rows are the left
        // eigen vectors of the rate matrix
        ieigenMatrix[category] = inverse( eigenMatrix[category] );
    } //for each category (except invariant)
}

