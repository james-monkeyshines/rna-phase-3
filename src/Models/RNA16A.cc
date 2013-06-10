#include "Models/RNA16A.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cstdio>

using namespace std;

//register the model to the model factory with the name RNA16A
RNA16A RNA16A::prototype( "RNA16A" );

RNA16A::RNA16A( const string & registrationName ) :
RNA16( registrationName ) {
//the private constructor is called by the prototype only
//RNA16A is a restriction of RNA16
}

RNA16A::RNA16A( ParametersSet & parameters )
: RNA16( parameters ) {
//the normal constructor is called by the prototype's clone method
//RNA16A is a restriction of RNA16
}


void RNA16A::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //The initialisation procedure for the 16A model is slightly different
    //than for the 16 model, therefore we redefine it
    //initialise the number of parameters
    numberRatesRatios = 4;
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
        vector<double> rates4;
        rates4.resize(4);
        
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
        //reference ratio
        double ref = ( rates16( 0, 2 ) + rates16( 3, 5 ) ) / 2;
        // 0 - single transition rate ratio
        rates4[0] = ( rates16( 0, 1 ) + rates16( 1, 2 ) + rates16( 3, 4 ) +
                                          rates16( 4, 5 ) );
        rates4[0] = rates4[0] / (4*ref);
        // 1 - double transversion rate ratio
        for ( int i = 0; i < 3; ++i ) {
            for ( int j = 3; j < 6; ++j ) {
                rates4[1] += rates16( i, j );
            }
        }
        rates4[1] = rates4[1]/(9*ref);
        // 2 - pair/mismatch rate ratio
        for ( int i = 0; i < 6; ++i ) {
            for ( int j = 6; j < 16; ++j ) {
                rates4[2] += rates16( i, j );
            }
        }
        rates4[2] = rates4[2]/(60*ref);
        // 3 - mismatch/mismatch rate ratio
        // to be reviewed because of the null rate
        for ( int i = 6; i < 16; ++i ) {
            for ( int j = 6; j < 16; ++j ) {
                rates4[3] += rates16( i, j );
            }
        }
        rates4[3] = rates4[3]/(60*ref);
        ratesRatios->initialisation( rates4 );
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}


// Destructor
RNA16A::~RNA16A() {
}


// The models unique name
string RNA16A::getName( void ) const {
    return ( "RNA16A" + MatrixModel::getName() );
}



void RNA16A::initMatrixIndex() {
    matrixIndex.resize( 16, 16 );  
    //the reference ratio is given the negative index of -2
    matrixIndex( 0, 2 ) = matrixIndex(3, 5) = -2;
    // 0 - ingle transition rate ratio
    matrixIndex(0, 1) = matrixIndex(1, 2) =
        matrixIndex(3, 4) = matrixIndex(4, 5) = 0;
    // 1 - double transversion rate ratio
    for ( int i = 0; i < 3; ++i ) {
        for ( int j = 3; j < 6; ++j ) {
            matrixIndex(i, j) = 1;
        }
    }

    string pair1;
    string pair2;
    // 3 - pair/mismatch rate ratio
    for ( int i = 0; i < 6; ++i ) {
        pair1 = getState(i);
        for ( int j = 6; j < 16; ++j ) {
            pair2 = getState(j);
            assert(pair1.size() == 2);
            assert(pair2.size() == 2);
            //if one nucleotide is the same rate = pair/mismatch rate ratio
            if ( (pair1[0]==pair2[0]) || (pair1[1]==pair2[1]) ){
                matrixIndex(i, j) = 2;
            }
            // else null rate ratio -> index = -1
            else{
                matrixIndex(i, j) = -1;
            }
        }
    }
    // 4 - mismatch/mismatch rate ratio
    for ( int i = 6; i < 15; ++i ) {
        pair1 = getState(i);
        for ( int j = i + 1; j < 16; ++j ) {
            pair2 = getState(j);
            assert(pair1.size() == 2);
            assert(pair2.size() == 2);
            //if one nucleotide is the same rate = pair/mismatch rate ratio
            if ( (pair1[0]==pair2[0]) || (pair1[1]==pair2[1]) ){
                matrixIndex(i, j) = 3;
            }
            // else null rate ratio, it is given the index of -1
            else{
                matrixIndex(i, j) = -1;
            }
        }
    }
    //copy on the other triangle of the matrix
    for ( int i = 1; i < 16; ++i ) {
        for ( int j = 0; j < i; ++j ) {
            matrixIndex(i, j) = matrixIndex(j, i);
        }
    }
    //initialise the diagonal
    for ( int i = 0 ; i < 16 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
}
