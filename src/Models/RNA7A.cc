#include "Models/RNA7A.h"
using namespace std;

RNA7A RNA7A::prototype( "RNA7A" );
RNA7A::RNA7A(){}
RNA7A::RNA7A( const string & registrationName ): RnaModel( registrationName ) {}
RNA7A::RNA7A( ParametersSet & parameters ): RnaModel( parameters, 7 ) {}

void RNA7A::initialisation( SequenceTable * sequenceTable, int modelId ) {
	// This initialisation is inherited by all 7-state models.
	// The setParameterNumbers function is redefined in other models; and
	// the condenseRates function may be overwritten, to account for
	// differences between model definitions.
	
    // Initialise the number of parameters
	bool bp_symmetry = setParameterNumbers();
	
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
        vector<double> countSite(64);
        vector<double> freqs;
        array2D<double> rates16x16;
        vector<double> rates;
        
        // Initialise the vector of frequencies from the empirical
		// frequencies of the 16 possible symbols.
        freqs16 = retrieveEmpiricalFrequencies( sequenceTable, modelId, & countSite );
		lnLAdjustment( countSite );
		// Calculate 7-state frequencies, possibly also applying base pair symmetry.
        freqs = mergeFrequencies( freqs16, bp_symmetry );
        frequencies->initialisation( freqs );
		
        // Estimate the rate ratios (and proportion of invariant sites,
		// if relevant) from the empirical sequences.
        if ( invariantCategory ) {
            rates16x16 = retrieveEmpiricalRates( sequenceTable, modelId, & proportionInvariantSites );
        }
        else {
            rates16x16 = retrieveEmpiricalRates( sequenceTable, modelId );
        }
		// Translate a 16x16 matrix of rates into a vector; restrictions of
		// the 7A model might merge rates to get a single transition rate etc.
		rates = condenseRates( rates16x16 );
        ratesRatios->initialisation( rates );
    }
	
	// The equivalency table needed to be 16-state until now, so that the individual
	// mismatch frequencies could be used in the likelihood adjustment. But for
	// subsequent analyses, it's easiest to now treat mismatches as a single state.
	setEquivalencyTable();
	
    // Initialise substitution rate matrix and eigen system.
    updateAverageRateVector();
    updateEigenMatrix();
}

bool RNA7A::setParameterNumbers() {
    numberRatesRatios = 20;
    numberFrequencies = 7;
	bool bp_symmetry = false;
	return bp_symmetry;
}

vector < double > RNA7A::mergeFrequencies( vector < double > & freqs16, bool bp_symmetry ) {
	// When generating the empirical frequencies, a pseudocount of 1e-6 was added
	// to each; so when merging the frequencies, subtract 9 * 1e-6 for the mismatches,
	// and 3 * 1e-6 if applying base-pair symmetry; then rescale so that the numbers total 1.
	vector < double > freqs(7);
	int last_index = 6;
	if (bp_symmetry) {
		freqs.resize(4);
		last_index = 3;
	}
	
    double sum = 0.0;
	for (int i = 0; i < 6 ; ++i ) {
		freqs[i % last_index] += freqs16[i];
		sum += freqs16[i];
	}
	for (int i = 6; i < 16 ; ++i ) {
		freqs[last_index] += freqs16[i];
		sum += freqs16[i];
	}
	freqs[last_index] -= 9e-6;
	sum -= 9e-6;
	if (bp_symmetry) {
		sum -= 3e-6;
	}
	
	for ( int i = 0; i < freqs.size(); ++i ) {
		freqs[i] /= sum;
	}
	
	return freqs;
}

vector < double > RNA7A::condenseRates( array2D < double > & rates16x16 ) {
	vector < double > rates(20);
	
	for ( int i = 0; i < 5 ; ++i ) {
		for ( int j = i + 1; j < 6 ; ++j ) {
			int index = matrixIndex(i, j);
			if ( index >= 0 ) {
				rates[index] = rates16x16( i, j );
			}
		}
	}
	
	for ( int i = 0; i < 6 ; ++i ) {
		double total_mismatch_rate = 0.0;
		for ( int j = 6; j < 16 ; ++j ) {
			total_mismatch_rate += rates16x16( i, j );
		}
		rates[matrixIndex(i, 6)] = total_mismatch_rate /= 10;
	}
	
	return rates;
}

void RNA7A::lnLAdjustment( vector < double > & countSite ) {
	// These next variables are globals...
	lnL_adj_eq_freq = 0.0;
	lnL_adj_emp_freq = 0.0;
	
	double sites = 0.0;
	double mismatches = 0.0;
	for ( int i = 0; i < 16; ++i ) {
		sites += countSite[i];
		if ( i > 5 ) {
			mismatches += countSite[i];
		}
	}
	
	if (mismatches > 0) {
		double mismatch_freq = mismatches / sites;
		lnL_adj_eq_freq += mismatches * log( 0.1 );
		
		for ( int i = 6; i < 16; ++i ) {
			if (countSite[i] > 0) {
				lnL_adj_emp_freq += log( pow(countSite[i]/sites, countSite[i]) / pow(mismatch_freq, countSite[i]) );
			}
		}
	}
}

void RNA7A::setEquivalencyTable() {
    int i, j;

    // Model states : AU, GU, GC, UA, UG, CG, MM
    equivalencyTable.resize( 64, 7 );

    // Clear the table
    for ( j = 0; j < 7; ++j ) {
        for ( i = 0; i < 64; ++i ) {
            equivalencyTable( i, j ) = 0.0;
        }
    }

    // Canonical states
    for ( i = 0; i < 6; ++i ) {
        equivalencyTable( i, i ) = 1.0;
    }

    // Mismatch states
    for ( i = 6; i < 16; ++i ) {
        equivalencyTable( i, 6 ) = 1.0;
    }

    equivalencyTable( 16, 0 ) = 1.0; // AX (AU)
	equivalencyTable( 16, 6 ) = 1.0; // AX (AA, AC or AG)
    equivalencyTable( 17, 6 ) = 1.0; // AR (AG or AA)
    equivalencyTable( 18, 0 ) = 1.0; // AY (AU)
	equivalencyTable( 18, 6 ) = 1.0; // AY (AC)

    equivalencyTable( 19, 5 ) = 1.0; // CX (CG)
	equivalencyTable( 19, 6 ) = 1.0; // CX (CA, CC or CU)
    equivalencyTable( 20, 5 ) = 1.0; // CR (CG)
	equivalencyTable( 20, 6 ) = 1.0; // CR (CA)
    equivalencyTable( 21, 6 ) = 1.0; // CY (CC or CU)

    equivalencyTable( 22, 1 ) = 1.0; // GX (GU)
	equivalencyTable( 22, 2 ) = 1.0; // GX (GC)
    equivalencyTable( 22, 6 ) = 1.0; // GX (GA or GG)
    equivalencyTable( 23, 6 ) = 1.0; // GR (GA or GG)
    equivalencyTable( 24, 1 ) = 1.0; // GY (GU)
	equivalencyTable( 24, 2 ) = 1.0; // GY (GC)

    equivalencyTable( 25, 3 ) = 1.0; // UX (UA)
	equivalencyTable( 25, 4 ) = 1.0; // UX (UG)
    equivalencyTable( 25, 6 ) = 1.0; // UX (UC or UU)
    equivalencyTable( 26, 3 ) = 1.0; // UR (UA)
	equivalencyTable( 26, 4 ) = 1.0; // UR (UG)
    equivalencyTable( 27, 6 ) = 1.0; // UY (UC or UU)

    equivalencyTable( 28, 3 ) = 1.0; // XA (UA)
	equivalencyTable( 28, 6 ) = 1.0; // XA (AA, CA or GA)
    equivalencyTable( 29, 6 ) = 1.0; // RA (GA or AA)
    equivalencyTable( 30, 3 ) = 1.0; // YA (UA)
	equivalencyTable( 30, 6 ) = 1.0; // YA (CA)

    equivalencyTable( 31, 2 ) = 1.0; // XC (GC)
	equivalencyTable( 31, 6 ) = 1.0; // XC (AC, CC or UC)
    equivalencyTable( 32, 2 ) = 1.0; // RC (GC)
	equivalencyTable( 32, 6 ) = 1.0; // RC (AC)
    equivalencyTable( 33, 6 ) = 1.0; // YC (CC or UC)

    equivalencyTable( 34, 4 ) = 1.0; // XG (UG)
	equivalencyTable( 34, 5 ) = 1.0; // XG (CG)
    equivalencyTable( 34, 6 ) = 1.0; // XG (AG or GG)
    equivalencyTable( 35, 6 ) = 1.0; // RG (AG or GG)
    equivalencyTable( 36, 4 ) = 1.0; // YG (UG)
	equivalencyTable( 36, 5 ) = 1.0; // YG (CG)

    equivalencyTable( 37, 0 ) = 1.0; // XU (AU)
	equivalencyTable( 37, 1 ) = 1.0; // XU (GU)
    equivalencyTable( 37, 6 ) = 1.0; // XU (CU or UU)
    equivalencyTable( 38, 0 ) = 1.0; // RU (AU)
	equivalencyTable( 38, 1 ) = 1.0; // RU (GU)
    equivalencyTable( 39, 6 ) = 1.0; // YU (CU or UU)

    for ( i = 0; i < 7; i++ ) {
        equivalencyTable( 40, i ) = 1.0; // XX
    }

    equivalencyTable( 41, 3 ) = 1.0; // XR (UA)
    equivalencyTable( 41, 4 ) = 1.0; // XR (UG)
    equivalencyTable( 41, 5 ) = 1.0; // XR (CG)
    equivalencyTable( 41, 6 ) = 1.0; // XR (AA, CA, GA, AG, GG)

    equivalencyTable( 42, 0 ) = 1.0; // XY (AU)
    equivalencyTable( 42, 1 ) = 1.0; // XY (GU)
    equivalencyTable( 42, 2 ) = 1.0; // XY (GC)
    equivalencyTable( 42, 6 ) = 1.0; // XY (AC, CC, UC, CU, UU)

    equivalencyTable( 43, 3 ) = 1.0; // YX (UA)
    equivalencyTable( 43, 4 ) = 1.0; // YX (UG)
    equivalencyTable( 43, 5 ) = 1.0; // YX (CG)
    equivalencyTable( 43, 6 ) = 1.0; // YX (CA, CC, CU, UC, UU)

    equivalencyTable( 44, 3 ) = 1.0; // YR (UA)
    equivalencyTable( 44, 4 ) = 1.0; // YR (UG)
    equivalencyTable( 44, 5 ) = 1.0; // YR (CG)
    equivalencyTable( 44, 6 ) = 1.0; // YR (CA)

    equivalencyTable( 45, 6 ) = 1.0; // YY (CC, CU, UC, UU)

    equivalencyTable( 46, 0 ) = 1.0; // RX (AU)
	equivalencyTable( 46, 1 ) = 1.0; // RX (GU)
    equivalencyTable( 46, 2 ) = 1.0; // RX (GC)
	equivalencyTable( 46, 6 ) = 1.0; // RX (AA, AC, AG, GA, GG)

    equivalencyTable( 47, 6 ) = 1.0; // RR (AA, AG, GA, GG)

    equivalencyTable( 48, 0 ) = 1.0; // RY (AU)
	equivalencyTable( 48, 1 ) = 1.0; // RY (GU)
    equivalencyTable( 48, 2 ) = 1.0; // RY (GC)
	equivalencyTable( 48, 6 ) = 1.0; // RY (AC)

    equivalencyTable( 49, 3 ) = 1.0; // -A (UA)
	equivalencyTable( 49, 6 ) = 1.0; // -A (AA, CA, GA)
    equivalencyTable( 50, 4 ) = 1.0; // -G (UG)
	equivalencyTable( 50, 5 ) = 1.0; // -G (CG)
    equivalencyTable( 50, 6 ) = 1.0; // -G (AG, GG)
    equivalencyTable( 51, 2 ) = 1.0; // -C (GC)
	equivalencyTable( 51, 6 ) = 1.0; // -C (AC, CC, UC)
    equivalencyTable( 52, 0 ) = 1.0; // -U (AU)
	equivalencyTable( 52, 1 ) = 1.0; // -U (GU)
    equivalencyTable( 52, 6 ) = 1.0; // -U (CU, UU)

    for ( i = 0; i < 7; i++ ) {
        equivalencyTable( 53, i ) = 1.0; // -X
    }

    equivalencyTable( 54, 3 ) = 1.0; // -R (UA)
    equivalencyTable( 54, 4 ) = 1.0; // -R (UG)
    equivalencyTable( 54, 5 ) = 1.0; // -R (CG)
    equivalencyTable( 54, 6 ) = 1.0; // -R (AA, CA, GA, AG, GG)

    equivalencyTable( 55, 0 ) = 1.0; // -Y (AU)
	equivalencyTable( 55, 1 ) = 1.0; // -Y (GU)
    equivalencyTable( 55, 2 ) = 1.0; // -Y (GC)
	equivalencyTable( 55, 6 ) = 1.0; // -Y (AC, CC, UC, CU, UU)

    equivalencyTable( 56, 0 ) = 1.0; // A- (AU)
	equivalencyTable( 56, 6 ) = 1.0; // A- (AA, AC, AG)
    equivalencyTable( 57, 1 ) = 1.0; // G- (GU)
	equivalencyTable( 57, 2 ) = 1.0; // G- (GC)
    equivalencyTable( 57, 6 ) = 1.0; // G- (GA, GG)
    equivalencyTable( 58, 5 ) = 1.0; // C- (CG)
	equivalencyTable( 58, 6 ) = 1.0; // C- (CA, CC, CU)
    equivalencyTable( 59, 3 ) = 1.0; // U- (UA)
	equivalencyTable( 59, 4 ) = 1.0; // U- (UG)
    equivalencyTable( 59, 6 ) = 1.0; // U- (UC, UU)

    equivalencyTable( 61, 0 ) = 1.0; // R- (AU)
	equivalencyTable( 61, 1 ) = 1.0; // R- (GU)
	equivalencyTable( 61, 2 ) = 1.0; // R- (GC)
    equivalencyTable( 61, 6 ) = 1.0; // R- (AA, AC, AG, GA, GG)

    equivalencyTable( 62, 3 ) = 1.0; // Y- (UA)
	equivalencyTable( 62, 4 ) = 1.0; // Y- (UG)
	equivalencyTable( 62, 5 ) = 1.0; // Y- (CG)
    equivalencyTable( 62, 6 ) = 1.0; // Y- (CA, CC, CU, UC, UU)

    for ( i = 0; i < 7; i++ ) {
        equivalencyTable( 60, i ) = 1.0;    // X-
        equivalencyTable( 63, i ) = 1.0;    // --
    }
}

// Destructor
RNA7A::~RNA7A() {}

string RNA7A::getName( void ) const {
    string name = "RNA7A" + MatrixModel::getName();
    return ( name );
}

void RNA7A::initMatrixIndex() {
    matrixIndex.resize( 7, 7 );
    int index = -1;
    for ( int i = 0 ; i < 6 ; ++i ) {
        for ( int j = i + 1 ; j < 7 ; ++j ) {
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
    assert( index == 19 );
    //initialise the diagonal
    for ( int i = 0 ; i < 7 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
}

string RNA7A::getState( unsigned int stateNumber, unsigned int ) const {
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
            return ( string( "MM" ) );
        default :
            return("");
    }
}

Model * RNA7A::clone( ParametersSet & parameters ) const {
    return new RNA7A( parameters );
}
