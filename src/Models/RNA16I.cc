#include "Models/RNA16I.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name RNA16I
RNA16I RNA16I::prototype( "RNA16I" );

RNA16I::RNA16I( const string & registrationName ) :
RNA16( registrationName ) {
//the private constructor is called by the prototype only
//RNA16I is a restriction of RNA16
}

RNA16I::RNA16I( ParametersSet & parameters )
: RNA16( parameters ) {
//the normal constructor is called by the prototype's clone method
//RNA16I is a restriction of RNA16
}


void RNA16I::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //The initialisation procedure for the RNA16I model is slightly different
    //than for the RNA16 model, therefore we redefine it
    //initialise the number of parameters
    numberRatesRatios = 5;
    numberFrequencies = 16;
    
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
        array2D<double> rates16;
        vector<double> rates5;
        rates5.resize(5);
        
        frequencies->initialisation( retrieveEmpiricalFrequencies( sequenceTable, modelId ) );
        //estimation of all the rate ratios (and proportion of invariant sites
        //if relevant from the empirical sequences
        if ( invariantCategory ) {
            rates16 = retrieveEmpiricalRates( sequenceTable, modelId, & proportionInvariantSites );
        }
        else {
            rates16 = retrieveEmpiricalRates( sequenceTable, modelId );
        }
        //fill the rate vector according to the estimation
        //reference ratio A<->G, it should be divided by 8
        //but since all the rates in rates5 are initialised from
        //8 rates in rates16 it is a useless computation
        double ref = ( rates16( 0, 1 ) + rates16( 2, 8 ) + rates16( 3, 4 ) +
                       rates16( 5, 11 )+ rates16( 6, 7 ) + rates16( 6, 9 ) +
                       rates16( 7, 10 )+ rates16( 9, 10 ) );
        // second transition rate ratio C<->U
        rates5[3] = ( rates16( 0, 8 ) + rates16( 1, 2 ) + rates16( 3, 11 ) +
                      rates16( 4, 5 )+ rates16( 12, 13 ) + rates16( 12, 14 ) +
                      rates16( 13, 15 )+ rates16( 14, 15 ) ) / ref;
        // first transversion rate ratio A<->C
        rates5[0] = ( rates16( 0, 13 ) + rates16( 2, 9 ) + rates16( 3, 14 ) +
                      rates16( 5, 7 )+ rates16( 6, 8 ) + rates16( 6, 11 ) +
                      rates16( 8, 12 )+ rates16( 11, 12 ) ) / ref;
        // second transversion rate ratio A<->U
        rates5[1] = ( rates16( 0, 6 ) + rates16( 0, 15 ) + rates16( 1, 9 ) +
                      rates16( 3, 6 )+ rates16( 3, 15 ) + rates16( 4, 7 ) +
                      rates16( 8, 14 )+ rates16( 11, 13 ) ) / ref;
        // third transversion rate ratio C<->G
        rates5[2] = ( rates16( 1, 13 ) + rates16( 2, 10 ) + rates16( 2, 12 ) +
                      rates16( 4, 14 )+ rates16( 5, 10 ) + rates16( 5, 12 ) +
                      rates16( 7, 8 )+ rates16( 9, 11 ) ) / ref;
        // fourth transversion rate ratio G<->U
        rates5[4] = ( rates16( 0, 7 ) + rates16( 1, 10 ) + rates16( 1, 15 ) +
                      rates16( 2, 14 )+ rates16( 3, 9 ) + rates16( 4, 10 ) +
                      rates16( 4, 15 )+ rates16( 5, 13 ) ) / ref;
        ratesRatios->initialisation( rates5 );
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}


// Destructor
RNA16I::~RNA16I() {
}


// The models unique name
string RNA16I::getName( void ) const {
    return ( "RNA16I" + MatrixModel::getName() );
}



void RNA16I::initMatrixIndex() {
    matrixIndex.resize( 16, 16 );  
    //initialise the diagonal (-3)
    for ( int i = 0 ; i < 16 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
    //fill the matrix for double substitution (-1 for rate = 0)
    for ( int i = 0; i < 15 ; ++i ) {
        for ( int j = i+1; j < 16 ; ++j ) {
            matrixIndex(i, j) = -1;
        }
    }
    
    //the reference ratio is given the negative index of -2
    matrixIndex( 0, 1 ) = matrixIndex(2, 8) =
        matrixIndex(3, 4) = matrixIndex(5, 11) =
        matrixIndex(6, 7) = matrixIndex(6, 9) =
        matrixIndex(7, 10) = matrixIndex(9, 10) = -2;
    // second transition rate ratio C<->U
    matrixIndex( 0, 8 ) = matrixIndex( 1, 2 ) =
        matrixIndex( 3, 11 ) = matrixIndex( 4, 5 ) =
        matrixIndex( 12, 13 ) = matrixIndex( 12, 14 ) =
        matrixIndex( 13, 15 ) = matrixIndex( 14, 15 ) = 3;
    matrixIndex( 0, 13 ) = matrixIndex( 2, 9 ) =
        matrixIndex( 3, 14 ) = matrixIndex( 5, 7 ) =
        matrixIndex( 6, 8 ) = matrixIndex( 6, 11 ) =
        matrixIndex( 8, 12 ) = matrixIndex( 11, 12 ) = 0;
   matrixIndex( 0, 6 ) = matrixIndex( 0, 15 ) =
        matrixIndex( 1, 9 ) = matrixIndex( 3, 6 ) =
        matrixIndex( 3, 15 ) = matrixIndex( 4, 7 ) =
        matrixIndex( 8, 14 ) = matrixIndex( 11, 13 ) = 1;
   matrixIndex( 1, 13 ) = matrixIndex( 2, 10 ) =
        matrixIndex( 2, 12 ) = matrixIndex( 4, 14 ) =
        matrixIndex( 5, 10 ) = matrixIndex( 5, 12 ) =
        matrixIndex( 7, 8 ) = matrixIndex( 9, 11 ) = 2;
    matrixIndex( 0, 7 ) = matrixIndex( 1, 10 ) =
        matrixIndex( 1, 15 ) = matrixIndex( 2, 14 ) =
        matrixIndex( 3, 9 ) = matrixIndex( 4, 10 ) =
        matrixIndex( 4, 15 ) = matrixIndex( 5, 13 ) = 4;
    
    //symetrize the matrix
    for ( int i = 0; i < 15 ; ++i ) {
        for ( int j = i+1; j < 16 ; ++j ) {
            matrixIndex(j, i) = matrixIndex(i, j);
        }
    }
}
