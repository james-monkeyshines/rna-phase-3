#include "Models/RNA16_bps.h"
using namespace std;

// This model is very similar to RNA16, with base pair symmetry for canonical pairs only.
RNA16_bps RNA16_bps::prototype( "RNA16_bps" );
RNA16_bps::RNA16_bps( const string & registrationName ): RNA16( registrationName ) {}
RNA16_bps::RNA16_bps( ParametersSet & parameters ): RNA16( parameters ) {}

void RNA16_bps::initialisation( SequenceTable * sequenceTable, int modelId ) {
    // Initialise the number of parameters
	numberRatesRatios = 119;
    numberFrequencies = 13;
    
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
        vector<double> freqs(13);
        array2D<double> rates16x16;
        vector<double> rates(119);
        
		freqs16 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        // Easiest way to avoid messing up these sums is to do it long-hand.
		// Subtract 1e-6 to remove redundant pseudocount done in 'retrieveEmpiricalFrequencies'.
		freqs[0] = freqs16[0] + freqs16[3] - 1e-6;
		freqs[1] = freqs16[1] + freqs16[4] - 1e-6;
		freqs[2] = freqs16[2] + freqs16[5] - 1e-6;
		freqs[3] = freqs16[6];
		freqs[4] = freqs16[7];
		freqs[5] = freqs16[8];
		freqs[6] = freqs16[9];
		freqs[7] = freqs16[10];
		freqs[8] = freqs16[11];
		freqs[9] = freqs16[12];
		freqs[10] = freqs16[13];
		freqs[11] = freqs16[14];
		freqs[12] = freqs16[15];
    	double sum = 0.0;
		for ( int i = 0; i < freqs.size(); ++i ) {
			sum += freqs[i];
		}
		for ( int i = 0; i < freqs.size(); ++i ) {
			freqs[i] /= sum;
		}
		frequencies->initialisation( freqs );
		
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
RNA16_bps::~RNA16_bps() {
}

string RNA16_bps::getName( void ) const {
    return ( "RNA16_bps" + MatrixModel::getName() );
}

double RNA16_bps::getFrequency( unsigned int residue,
    unsigned int rateCategory, unsigned int ) const {
    switch (residue){
        case 0: case 3:
            return (*frequencies)[rateCategory][0]/2.0;
        case 1: case 4:
            return (*frequencies)[rateCategory][1]/2.0;
        case 2: case 5:
            return (*frequencies)[rateCategory][2]/2.0;
        case 6:
            return (*frequencies)[rateCategory][3];
        case 7:
            return (*frequencies)[rateCategory][4];
        case 8:
            return (*frequencies)[rateCategory][5];
        case 9:
            return (*frequencies)[rateCategory][6];
        case 10:
            return (*frequencies)[rateCategory][7];
        case 11:
            return (*frequencies)[rateCategory][8];
        case 12:
            return (*frequencies)[rateCategory][9];
        case 13:
            return (*frequencies)[rateCategory][10];
        case 14:
            return (*frequencies)[rateCategory][11];
        case 15:
            return (*frequencies)[rateCategory][12];
    }
    assert("wrong residue in RNA16_bps:getFrequency()"==0);
    return 0.0;
}

string RNA16_bps::getFrequencyState( unsigned int freqState, unsigned int ) const{
	assert(freqState<10);
	switch (freqState){
	case 0:
		return getState(0)+'+'+getState(3);
	case 1:
		return getState(1)+'+'+getState(4);
	case 2:
		return getState(2)+'+'+getState(5);
	case 3:
		return getState(6);
	case 4:
		return getState(7);
	case 5:
		return getState(8);
	case 6:
		return getState(9);
	case 7:
		return getState(10);
	case 8:
		return getState(11);
	case 9:
		return getState(12);
	case 10:
		return getState(13);
	case 11:
		return getState(14);
	case 12:
		return getState(15);
	}
}


