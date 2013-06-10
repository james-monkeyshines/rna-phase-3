#include "Models/RNA16.h"
using namespace std;

RNA16 RNA16::prototype( "RNA16" );
RNA16::RNA16( const string & registrationName ): RnaModel( registrationName ) {}
RNA16::RNA16( ParametersSet & parameters ): RnaModel( parameters, 16 ) {}

void RNA16::initialisation( SequenceTable * sequenceTable, int modelId ) {
    // Initialise the number of parameters
	numberRatesRatios = 119;
    numberFrequencies = 16;
	
    // Initialise the (16-state) equivalency table,
	// for correspondance between symbols and states.
    initEquivalencyTable();
    
    // Initialise the matrixIndex for correspondance between
    // a matrix element and an index in the rate ratios vector.
    initMatrixIndex();
	
    // Basic initialisation method.
    RnaModel::initialisation( sequenceTable, modelId );
    
    // Initialise rate ratios and frequencies according to the sequences.
    if ( sequenceTable ) {
        vector<double> freqs16;
        array2D<double> rates16x16;
        vector<double> rates(119);
        
        // Initialise the vector of frequencies from the empirical
		// frequencies of the 16 possible symbols.
       	freqs16 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        frequencies->initialisation( freqs16 );
		
        // Estimate the rate ratios (and proportion of invariant sites,
		// if relevant) from the empirical sequences.
        if ( invariantCategory ) {
            rates16x16 = retrieveEmpiricalRates( sequenceTable, modelId, & proportionInvariantSites );
        }
        else {
            rates16x16 = retrieveEmpiricalRates( sequenceTable, modelId );
        }
        for ( int i = 0; i < 15 ; ++i ) {
            for ( int j = i + 1; j < 16 ; ++j ) {
                int index = matrixIndex(i, j);
                if ( index >= 0 ) {
                    rates[index] = rates16x16( i, j );
                }
            }
        }
        ratesRatios->initialisation( rates );
    }
	
    // Initialise substitution rate matrix and eigen system.
    updateAverageRateVector();
    updateEigenMatrix();
}

// Destructor
RNA16::~RNA16() {
}

string RNA16::getName( void ) const {
    return ( "RNA16" + MatrixModel::getName() );
}

void RNA16::initMatrixIndex() {
    matrixIndex.resize( 16, 16 );
    int index = -1;
    for ( int i = 0 ; i < 15 ; ++i ) {
        for ( int j = i + 1 ; j < 16 ; ++j ) {
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
    assert( index == 118 );
    //initialise the diagonal
    for ( int i = 0 ; i < 16 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
}

string RNA16::getState( unsigned int stateNumber, unsigned int ) const {
    switch ( stateNumber ) {
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
        case 6 :
            return ( string( "AA" ) );
        case 7 :
            return ( string( "AG" ) );
        case 8 :
            return ( string( "AC" ) );
        case 9 :
            return ( string( "GA" ) );
        case 10 :
            return ( string( "GG" ) );
        case 11 :
            return ( string( "CA" ) );
        case 12 :
            return ( string( "CC" ) );
        case 13 :
            return ( string( "CU" ) );
        case 14 :
            return ( string( "UC" ) );
        case 15 :
            return ( string( "UU" ) );
        default :
            return("");
    }
}

