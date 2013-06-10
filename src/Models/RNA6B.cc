#include "Models/RNA6B.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name RNA6B
RNA6B RNA6B::prototype( "RNA6B" );

RNA6B::RNA6B( const string & registrationName ) :
RNA6A( registrationName ) {
//the private constructor is called by the prototype only
//RNA6B is a restriction of RNA6A
}

RNA6B::RNA6B( ParametersSet & parameters )
: RNA6A( parameters ) {
//the normal constructor is called by the prototype's clone method
//RNA6B is a restriction of RNA6A
}


void RNA6B::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //The initialisation procedure for the 6B model is slightly different
    //than for the 6A model, therefore we redefine it
    //initialise the number of parameters
    numberRatesRatios = 2;
    numberFrequencies = 6;

    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();
    
    //initialise the matrixIndex for correspondance between a matrix element
    //and an index in the rate ratios vector
    initMatrixIndex();

    //call the basic initialisation method
    RnaModel::initialisation( sequenceTable, modelId );
    
    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if ( sequenceTable ) {
        vector<double> freq6;
        array2D<double> rates6;
        vector<double> rates2;
        rates2.resize(2);

        //initialise the vector of frequencies from the empirical frequencies
        //of the 16 possible symbols (normalize to remove MM)
        freq6 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        frequencies->initialisation( freq6 );
        
        //estimation of all the rate ratios (and proportion of invariant sites
        //if relevant from the empirical sequences
        if ( invariantCategory ) {
            rates6 = retrieveEmpiricalRates( sequenceTable, modelId, & proportionInvariantSites );
        }
        else {
            rates6 = retrieveEmpiricalRates( sequenceTable, modelId );
        }

        //fill the rate vector according to the estimation
        // reference ratio
        double ref = ( rates6( 0, 2 ) + rates6( 3, 5 ) ) / 2;
        // 0 - single transition rate ratio
        rates2[0] = ( rates6( 0, 1 ) + rates6( 1, 2 ) + rates6( 3, 4 ) +
                      rates6( 4, 5 ) );
        rates2[0] = rates2[0] / (4*ref);
        // 1 - double transversion rate ratio
        for ( int i = 0; i < 3; ++i ) {
            for ( int j = 3; j < 6; ++j ) {
                rates2[1] += rates6( i, j );
            }
        }
        rates2[1] = rates2[1]/(9*ref);
        ratesRatios->initialisation( rates2 );
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}


// Destructor
RNA6B::~RNA6B() {
}


// The models unique name
string RNA6B::getName( void ) const {
    return ( "RNA6B" + MatrixModel::getName() );
}




void RNA6B::initMatrixIndex() {
    matrixIndex.resize( 6, 6 );
    //the reference ratio is given the negative index of -2
    matrixIndex( 0, 2 ) = matrixIndex(3, 5) = -2;
    // 0 - single transition rate ratio
    matrixIndex(0, 1) = matrixIndex(1, 2) =
        matrixIndex(3, 4) = matrixIndex(4, 5) = 0;
    // 1 - double transversion rate ratio
    for ( int i = 0; i < 3; ++i ) {
        for ( int j = 3; j < 6; ++j ) {
            matrixIndex(i, j) = 1;
        }
    }
    //copy on the other triangle of the matrix
    for ( int i = 1; i < 6; ++i ) {
        for ( int j = 0; j < i; ++j ) {
            matrixIndex(i, j) = matrixIndex(j, i);
        }
    }
    //initialise the diagonal
    for ( int i = 0 ; i < 6 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
}


