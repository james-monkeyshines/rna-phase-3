#include "Models/RNA6C.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name RNA6B
RNA6C RNA6C::prototype( "RNA6C" );

RNA6C::RNA6C( const string & registrationName ) :
RNA6B( registrationName ) {
//the private constructor is called by the prototype only
//RNA6C is a restriction of RNA6B
}

RNA6C::RNA6C( ParametersSet & parameters )
: RNA6B( parameters ) {
//the normal constructor is called by the prototype's clone method
//RNA6C is a restriction of RNA6B
}


void RNA6C::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //The initialisation procedure for the 6C model is slightly different
    //than for the 6B model, therefore we redefine it
    //initialise the number of parameters
    numberRatesRatios = 2;
    numberFrequencies = 3;

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
        vector<double> freq3;
        freq3.resize(3);

        //initialise the vector of frequencies from the empirical frequencies
        //of the 16 possible symbols
        freq6 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        freq3[0] = freq6[0] + freq6[3];
        freq3[1] = freq6[1] + freq6[4];
        freq3[2] = freq6[2] + freq6[5];
        frequencies->initialisation( freq3 );
        
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
RNA6C::~RNA6C() {
}


// The models unique name
string RNA6C::getName( void ) const {
    return ( "RNA6C" + MatrixModel::getName() );
}



double RNA6C::getFrequency( unsigned int residue,
    unsigned int rateCategory, unsigned int ) const {
    assert(residue<6);
    return (*frequencies)[rateCategory][residue%3]/2.0;
    return 0.0;
}

string RNA6C::getFrequencyState( unsigned int freqState, unsigned int ) const{
    assert(freqState<3);
    return getState(freqState)+'+'+getState(freqState+3);
}
