#include "Models/RNA16C.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <numeric>

using namespace std;

//register the model to the model factory with the name RNA16A
RNA16C RNA16C::prototype( "RNA16C" );

RNA16C::RNA16C( const string & registrationName ) :
RNA16A( registrationName ) {
//the private constructor is called by the prototype only
//RNA16C is a restriction of RNA16A
}

RNA16C::RNA16C( ParametersSet & parameters )
: RNA16A( parameters ) {
//the normal constructor is called by the prototype's clone method
//RNA16C is a restriction of RNA16A
}


void RNA16C::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //The initialisation procedure for the 16C model is slightly different
    //than for the 16A model, therefore we redefine it
    //initialise the number of parameters
    numberRatesRatios = 4;
    numberFrequencies = 7;
    
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
        vector<double> freq16;
        vector<double> freq7;
        freq7.resize(7);
        array2D<double> rates16;
        vector<double> rates4;
        rates4.resize(4);
        
        freq16 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        for ( unsigned int i = 0; i < 6; ++i ){
            freq7[i]=freq16[i];
        }
        freq7[6] = 0.0;
        freq7[6] = 1.0 - accumulate(freq7.begin(),freq7.end(),0.0);
        frequencies->initialisation( freq7 );
        
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
RNA16C::~RNA16C() {
}


// The models unique name
string RNA16C::getName( void ) const {
    return ( "RNA16C" + MatrixModel::getName() );
}


double RNA16C::getFrequency( unsigned int residue,
    unsigned int rateCategory, unsigned int ) const {
    assert(residue<16);
    if (residue<6){
        return (*frequencies)[rateCategory][residue];
    }
    return (*frequencies)[rateCategory][6]/10.0;
}

string RNA16C::getFrequencyState( unsigned int freqState, unsigned int ) const{
    assert(freqState<7);
    if (freqState==6) return "MM";
    else return getState(freqState);
}
