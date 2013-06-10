#include <assert.h>

#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Models/Perturbator.h"
#include "Models/THREESTATE.h"

#include "Sequence/SequenceTable.h"

#include "Util/statlib.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"

#define SQR(n) (n*n)

using namespace std;

//register the model to the model factory with the name REV
THREESTATE THREESTATE::prototype( "THREESTATE" );

THREESTATE::THREESTATE( const string & registrationName ) : ThreeStateModel( registrationName ){
//the private constructor is called by the prototype only
}


THREESTATE::THREESTATE( ParametersSet & parameters ) : ThreeStateModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void THREESTATE::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 2;
    numberFrequencies = 3;
    
    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();
    
    //initialise the matrixIndex for correspondance between a matrix element
    //and an index in the rate ratios vector
    initMatrixIndex();
            
    //call the basic initialisation method
    ThreeStateModel::initialisation( sequenceTable, modelId );
        
    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if (sequenceTable){
        //estimation of the model frequencies from the empirical frequencies
        vector<double> freq3 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        frequencies->initialisation( freq3 );

        //estimation of all the rate ratios (and proportion of invariant sites
        //if relevant) from the empirical sequences
        vector < double > rates2;
        if ( invariantCategory ) {
            rates2 = retrieveEmpiricalRates( sequenceTable, modelId,
            & proportionInvariantSites );
        }
        else {
            rates2 = retrieveEmpiricalRates( sequenceTable, modelId );
        }
        ratesRatios->initialisation( rates2 );
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}  

THREESTATE::~THREESTATE() {
}

string THREESTATE::getName( void ) const {
    return ( "THREESTATE" + MatrixModel::getName() );
}

void THREESTATE::initMatrixIndex() {
    matrixIndex.resize( 3, 3 );
    for ( int i = 0 ; i < 3 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
    //reference
    matrixIndex(0, 2) = matrixIndex(2, 0) = -2;
    //transversion
    matrixIndex(0, 1) = matrixIndex(1, 0) = 0;
    matrixIndex(1, 2) = matrixIndex(2, 1) = 1;
}


void THREESTATE::setEigenMatrix() {
    vector< double > junk( 3 );
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
          ++category ) {
        // fill the rateMatrix (cf REV.h)
        int freqCat = frequencies->frequenciesCat(category+invariantCategory);
        int rateCat = ratesRatios->ratesRatiosCat(category);
        rateMatrix[category] ( 0, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1] *
                (*ratesRatios)(rateCat)[0]; // A->C/T
        rateMatrix[category] ( 0, 2 ) =
            substitutionRate[category] * (*frequencies)(freqCat)[2]; //A->G
        rateMatrix[category] ( 0, 0 ) = - rateMatrix[category] ( 0, 1 )
                                        - rateMatrix[category] ( 0, 2 );
        rateMatrix[category] ( 1, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0] *
                (*ratesRatios)(rateCat)[0]; // C/T->A
        rateMatrix[category] ( 1, 2 ) =
            substitutionRate[category] * (*frequencies)(freqCat)[2] *
                (*ratesRatios)(rateCat)[1]; // C/T->G
        rateMatrix[category] ( 1, 1 ) = - rateMatrix[category] ( 1, 0 )
                                        - rateMatrix[category] ( 1, 2 );
        
        rateMatrix[category] ( 2, 0 ) =
            substitutionRate[category] * (*frequencies)(freqCat)[0]; //G->A
        rateMatrix[category] ( 2, 1 ) =
            substitutionRate[category] * (*frequencies)(freqCat)[1];
                (*ratesRatios)(rateCat)[1]; // G->C/T
        rateMatrix[category] ( 2, 2 ) = - rateMatrix[category] ( 2, 0 )
                                        - rateMatrix[category] ( 2, 1 );

        for ( int i = 0; i < 3; ++i ) {
            eigenValues[category] [i] = 0.0;
            junk[i] = 0.0;
            for ( int j = 0; j < 3; ++j ) {
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
