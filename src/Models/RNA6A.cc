#include "Models/RNA6A.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name RNA6A
RNA6A RNA6A::prototype( "RNA6A" );

RNA6A::RNA6A( const string & registrationName ) :
RnaModel( registrationName ) {
//the private constructor is called by the prototype only
}

RNA6A::RNA6A( ParametersSet & parameters )
: RnaModel( parameters, 6 ) {
//the normal constructor is called by the prototype's clone method
//RNA6A is a RnaModel, its matrix size is 6*6
}


void RNA6A::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 14;
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
        vector<double> rates14;
        rates14.resize(14);

        //initialise the vector of frequencies from the empirical frequencies
        //of the 6 possible symbols
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
        for ( int i = 0; i < 5 ; ++i ) {
            for ( int j = i + 1; j < 6 ; ++j ) {
                int index = matrixIndex(i, j);
                //if the rate is not the reference rate (0, 2)
                if ( index >= 0 ) {
                    rates14[index] = rates6( i, j );
                }
            }
        }
        ratesRatios->initialisation( rates14 );
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}


// Destructor
RNA6A::~RNA6A() {
}


// The models unique name
string RNA6A::getName( void ) const {
    return ( "RNA6A" + MatrixModel::getName() );
}




void RNA6A::initMatrixIndex() {
    matrixIndex.resize( 6, 6 );
    int index = -1;
    for ( int i = 0 ; i < 5 ; ++i ) {
        for ( int j = i + 1 ; j < 6 ; ++j ) {
            if ( ( i == 0 ) && ( j == 2 ) ) {
                //the reference ratio is given the negative index of -2
                matrixIndex(i, j) = matrixIndex(j, i) = -2;
            }
            else {
                ++index;
                matrixIndex(i, j) = matrixIndex(j, i) = index;
            }
        }
    }
    assert( index == 13 );
    //initialise the diagonal
    for ( int i = 0 ; i < 6 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
}


void RNA6A::initEquivalencyTable() {
    int i, j;

    // model states : AU, GU, GC, UA, UG, CG .
    equivalencyTable.resize( 64, 6 );

    // Set the patterns for valid states
    for ( i = 0; i < 6; ++i ) {
        for ( j = 0; j < 6; ++j ) {
            equivalencyTable( i, j ) = 0.0;
        }
        equivalencyTable( i, i ) = 1.0;
    }

    // Match all other pairs (unknown)
    for ( i = 6; i < 64; ++i ) {
        for ( j = 0; j < 6; ++j ) {
            equivalencyTable( i, j ) = 1.0;
        }
    }
}



string RNA6A::getState( unsigned int stateNumber, unsigned int ) const {
    switch ( stateNumber ) {
        // Watson-Crick base pairs
        case 0 :
            return ( string( "AU" ) );
        case 1 :
            return ( string( "GU" ) );
        case 2 :
            return ( string( "GC" ) );
        case 3 :
            return ( string( "UA" ) );
        case 4 :
            return ( string( "UG" ) );
        case 5 :
            return ( string( "CG" ) );
        default :
            return("");
    }
}

Model * RNA6A::clone( ParametersSet & parameters ) const{
    return new RNA6A( parameters );
}
