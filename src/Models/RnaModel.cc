#include "Models/RnaModel.h"

#include <iostream>
#include <iomanip>

#include "Sequence/SequenceTable.h"

#include "PatternDesign/Factory.h"

#include "PatternDesign/Singleton.h"
#include "Util/statlib.h"
#include "Util/matrixmath.h"
#include "Util/array2D.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"


using namespace std;

RnaModel::RnaModel(){
}

RnaModel::RnaModel( const string & registrationName ) :
        MatrixModel( registrationName ){
}

RnaModel::RnaModel( ParametersSet & parameters, unsigned int matrixSize ) :
        MatrixModel( parameters ) {

    this->matrixSize = matrixSize;

    initMatrix( matrixSize );
}

RnaModel::~RnaModel() {
}


void RnaModel::initEquivalencyTable() {
    int i, j;
    
    equivalencyTable.resize( 64, 16 );

    // Clear the table.
    for ( i = 0; i < 64; ++i ) {
        for ( j = 0; j < 16; ++j ) {
            equivalencyTable( i, j ) = 0.0;
        }
    }

    //valid states
    for ( i = 0; i < 16; ++i ) {
        equivalencyTable( i, i ) = 1.0;
    }
    
    equivalencyTable( 16, 0 ) = equivalencyTable( 16, 6 ) =
    equivalencyTable( 16, 7 ) = equivalencyTable( 16, 8 ) = 1.0; // AX
    equivalencyTable( 17, 6 ) = equivalencyTable( 17, 7 ) = 1.0; // AR (AG or AA)
    equivalencyTable( 18, 0 ) = equivalencyTable( 18, 8 ) = 1.0; // AY (AU or AC)

    equivalencyTable( 19, 5 ) = equivalencyTable( 19, 11 ) =
    equivalencyTable( 19, 12 ) = equivalencyTable( 19, 13 ) = 1.0; // CX
    equivalencyTable( 20, 5 ) = equivalencyTable( 20, 11 ) = 1.0; // CR (CG or CA)
    equivalencyTable( 21, 12 ) = equivalencyTable( 21, 13 ) = 1.0; // CY (CC or CU)

    equivalencyTable( 22, 1 ) = equivalencyTable( 22, 2 ) =
    equivalencyTable( 22, 9 ) = equivalencyTable( 22, 10 ) = 1.0; // GX
    equivalencyTable( 23, 9 ) = equivalencyTable( 23, 10 ) = 1.0; // GR (GA or GG)
    equivalencyTable( 24, 1 ) = equivalencyTable( 24, 2 ) = 1.0; // GY (GC or GU)

    equivalencyTable( 25, 3 ) = equivalencyTable( 25, 4 ) =
    equivalencyTable( 25, 14 ) = equivalencyTable( 25, 15 ) = 1.0; // UX
    equivalencyTable( 26, 3 ) = equivalencyTable( 26, 4 ) = 1.0; // UR (UA or UG)
    equivalencyTable( 27, 14 ) = equivalencyTable( 27, 15 ) = 1.0; // UY (UC or UU)

    equivalencyTable( 28, 3 ) = equivalencyTable( 28, 6 ) =
    equivalencyTable( 28, 9 ) = equivalencyTable( 28, 11 ) = 1.0; // XA
    equivalencyTable( 29, 9 ) = equivalencyTable( 29, 10 ) = 1.0; // RA (GA or GG)
    equivalencyTable( 30, 3 ) = equivalencyTable( 30, 11 ) = 1.0; // YA (CA or UA)

    equivalencyTable( 31, 2 ) = equivalencyTable( 31, 8 ) =
    equivalencyTable( 31, 12 ) = equivalencyTable( 31, 14 ) = 1.0; // XC
    equivalencyTable( 32, 2 ) = equivalencyTable( 32, 8 ) = 1.0; // RC (GC or AC)
    equivalencyTable( 33, 12 ) = equivalencyTable( 33, 14 ) = 1.0; // YC (CC or UC)

    equivalencyTable( 34, 4 ) = equivalencyTable( 34, 5 ) =
    equivalencyTable( 34, 7 ) = equivalencyTable( 34, 10 ) = 1.0; // XG
    equivalencyTable( 35, 7 ) = equivalencyTable( 35, 10 ) = 1.0; // RG (AG or GG)
    equivalencyTable( 36, 4 ) = equivalencyTable( 36, 5 ) = 1.0; // YG (CG or UG)

    equivalencyTable( 37, 0 ) = equivalencyTable( 37, 1 ) =
    equivalencyTable( 37, 13 ) = equivalencyTable( 37, 15 ) = 1.0; // XU
    equivalencyTable( 38, 0 ) = equivalencyTable( 38, 1 ) = 1.0; // RU (AU or GU)
    equivalencyTable( 39, 13 ) = equivalencyTable( 39, 15 ) = 1.0; // YU (CU or UU)


    equivalencyTable( 41, 3 ) = equivalencyTable( 41, 4 ) =
    equivalencyTable( 41, 5 ) = equivalencyTable( 41, 6 ) =
    equivalencyTable( 41, 7 ) = equivalencyTable( 41, 9 ) =
    equivalencyTable( 41, 10 ) = equivalencyTable( 41, 11 ) = 1.0; //XR (XA,XG)

    equivalencyTable( 42, 0 ) = equivalencyTable( 42, 1 ) =
    equivalencyTable( 42, 2 ) = equivalencyTable( 42, 8 ) =
    equivalencyTable( 42, 12 ) = equivalencyTable( 42, 13 ) =
    equivalencyTable( 42, 14 ) = equivalencyTable( 42, 15 ) = 1.0; //XY (XC,XU)

    equivalencyTable( 43, 3 ) = equivalencyTable( 43, 4 ) =
    equivalencyTable( 43, 5 ) = equivalencyTable( 43, 11 ) =
    equivalencyTable( 43, 12 ) = equivalencyTable( 43, 13 ) =
    equivalencyTable( 43, 14 ) = equivalencyTable( 43, 15 ) = 1.0; //YX (CX,UX)

    equivalencyTable( 44, 3 ) = equivalencyTable( 44, 4 ) =
    equivalencyTable( 44, 5 ) = equivalencyTable( 44, 11 ) = 1.0; //YR (CA,CG,UA,UG)

    equivalencyTable( 45, 2 ) = equivalencyTable( 45, 5 ) =
    equivalencyTable( 45, 10 ) = equivalencyTable( 45, 12 ) = 1.0; //YY (CC,CG,GC,GG)

    equivalencyTable( 46, 0 ) = equivalencyTable( 46, 1 ) =
    equivalencyTable( 46, 2 ) = equivalencyTable( 46, 6 ) =
    equivalencyTable( 46, 7 ) = equivalencyTable( 46, 8 ) =
    equivalencyTable( 46, 9 ) = equivalencyTable( 46, 10 ) = 1.0; //RX (AX,GX)

    equivalencyTable( 47, 6 ) = equivalencyTable( 47, 7 ) =
    equivalencyTable( 47, 9 ) = equivalencyTable( 47, 10 ) = 1.0; //RR (AA,AG,GA,GG)

    equivalencyTable( 48, 0 ) = equivalencyTable( 48, 1 ) =
    equivalencyTable( 48, 2 ) = equivalencyTable( 48, 8 ) = 1.0; //RY (AC,AU,GC,GU)

    equivalencyTable( 49, 6 ) = equivalencyTable( 49, 9 ) =
    equivalencyTable( 49, 11 ) = equivalencyTable( 49, 13 ) = 1.0; //-A

    equivalencyTable( 50, 7 ) = equivalencyTable( 50, 10 ) =
    equivalencyTable( 50, 5 ) = equivalencyTable( 50, 4 ) = 1.0; //-G

    equivalencyTable( 51, 8 ) = equivalencyTable( 51, 2 ) =
    equivalencyTable( 51, 12 ) = equivalencyTable( 51, 14 ) = 1.0; //-C

    equivalencyTable( 52, 0 ) = equivalencyTable( 52, 1 ) =
    equivalencyTable( 52, 13 ) = equivalencyTable( 52, 15 ) = 1.0; //-U


    equivalencyTable( 54, 3 ) = equivalencyTable( 54, 4 ) =
    equivalencyTable( 54, 5 ) = equivalencyTable( 54, 6 ) =
    equivalencyTable( 54, 7 ) = equivalencyTable( 54, 9 ) =
    equivalencyTable( 54, 10 ) = equivalencyTable( 54, 11 ) = 1.0; //-R (-A,-G)

    equivalencyTable( 55, 0 ) = equivalencyTable( 55, 1 ) =
    equivalencyTable( 55, 2 ) = equivalencyTable( 55, 8 ) =
    equivalencyTable( 55, 12 ) = equivalencyTable( 55, 13 ) =
    equivalencyTable( 55, 14 ) = equivalencyTable( 55, 15 ) = 1.0; //-Y (-C,-U)

    equivalencyTable( 56, 0 ) = equivalencyTable( 56, 6 ) =
    equivalencyTable( 56, 7 ) = equivalencyTable( 56, 8 ) = 1.0; //A-

    equivalencyTable( 57, 1 ) = equivalencyTable( 57, 2 ) =
    equivalencyTable( 57, 9 ) = equivalencyTable( 57, 10 ) = 1.0; //G-

    equivalencyTable( 58, 5 ) = equivalencyTable( 58, 11 ) =
    equivalencyTable( 58, 12 ) = equivalencyTable( 58, 13 ) = 1.0; //C-

    equivalencyTable( 59, 3 ) = equivalencyTable( 59, 4 ) =
    equivalencyTable( 59, 14 ) = equivalencyTable( 59, 15 ) = 1.0; //U-

    equivalencyTable( 61, 0 ) = equivalencyTable( 61, 1 ) =
    equivalencyTable( 61, 2 ) = equivalencyTable( 61, 6 ) =
    equivalencyTable( 61, 7 ) = equivalencyTable( 61, 8 ) =
    equivalencyTable( 61, 9 ) = equivalencyTable( 61, 10 ) = 1.0; //R- (A-,G-)

    equivalencyTable( 62, 3 ) = equivalencyTable( 62, 4 ) =
    equivalencyTable( 62, 5 ) = equivalencyTable( 62, 11 ) =
    equivalencyTable( 62, 12 ) = equivalencyTable( 62, 13 ) =
    equivalencyTable( 62, 14 ) = equivalencyTable( 62, 15 ) = 1.0; //Y- (C-,U-)

    for ( i = 0; i < 16; ++i ) {
        equivalencyTable( 40, i ) = 1.0; //XX
        equivalencyTable( 53, i ) = 1.0; //-X
        equivalencyTable( 60, i ) = 1.0; //X-
        equivalencyTable( 63, i ) = 1.0; //--
    }
}

array2D < double > RnaModel::retrieveEmpiricalRates( SequenceTable * sequenceTable, int modelId, double * propInvariant ) {
    array2D < double > countRates( 64, 64 );
    array2D < double > rates( 16, 16 );
    double reference;

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair< string, unsigned int > > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s1, s2;

    for ( unsigned int i = 0 ; i < 64; ++i ) {
        for ( unsigned int j = 0 ; j < 64; ++j ) {
			countRates( i, j ) = 0.0;
        }
    }
	
    for ( int site = 0; site < (int)(sequences.numberColumns()); ++site ) {
        for ( unsigned int seq1 = 0; seq1 < sequenceTable->getNumberSpecies()-1; ++seq1 ) {
            s1 = getSymbolNumber( sequences( seq1, site ) );
            if (s1==getNumberSymbols()) {
                errorSymbol( sequenceTable, modelId, seq1, site );
            }
			
            for ( unsigned int seq2 = seq1 + 1; seq2 < sequenceTable->getNumberSpecies(); ++seq2 ) {
               s2 = getSymbolNumber( sequences( seq2, site ) );
                if (s2 == getNumberSymbols()){
                    errorSymbol( sequenceTable, modelId, seq2, site );
                }
				if (s1 < s2) {
					countRates( s1, s2 )++;
				} else if (s1 > s2) {
					countRates( s2, s1 )++;
				}
			}
		}
	}
	
	// Exclude gaps from the rate calculations.
	for ( int s1 = 0; s1 < 49; ++s1 ) {
		for ( int s2 = s1+1; s2 < 49; ++s2 ) {
			if (!(s1 < 16 && s2 < 16) && s1 != s2) {
				if (countRates( s1, s2 ) > 0) {
					double sum = 0.0;
					for ( int nuc1 = 0; nuc1 < 16; ++nuc1 ) {
						for ( int nuc2 = 0; nuc2 < 16; ++nuc2 ) {
							if ( ( equivalencyTable(s1, nuc1) == 1 ) &&
								 ( equivalencyTable(s2, nuc2) == 1 ) ) {
								sum++;
							}
						}
					}
					for ( int nuc1 = 0; nuc1 < 16; ++nuc1 ) {
						for ( int nuc2 = 0; nuc2 < 16; ++nuc2 ) {
							if ( ( equivalencyTable(s1, nuc1) == 1 ) &&
								 ( equivalencyTable(s2, nuc2) == 1 ) ) {
								countRates( nuc1, nuc2 ) += countRates( s1, s2 ) * (1/sum);
							}
						}
					}
				}
			}
		}
	}
	
    double sum = 0.0;
	for ( int nuc1 = 0; nuc1 < 16; ++nuc1 ) {
		for ( int nuc2 = nuc1+1; nuc2 < 16; ++nuc2 ) {
			// Use a pseudocount to ensure that we don't start with zero-valued
			// rates, because the optimizer can't then escape from them.
			countRates( nuc1, nuc2 ) += 1;
			sum += countRates( nuc1, nuc2 );
		}
	}
	
    reference = countRates( 0, 2 ) / sum;
	for ( int nuc1 = 0; nuc1 < 16; ++nuc1 ) {
		for ( int nuc2 = nuc1; nuc2 < 16; ++nuc2 ) {
			rates( nuc1, nuc2 ) = countRates( nuc1, nuc2 ) / (sum * reference);
			rates( nuc2, nuc1 ) = rates( nuc1, nuc2 );
        }
    }
	
    if ( propInvariant != NULL ) {
        int numberInvariant = 0;
        for ( vector< pair< string, unsigned int > >::const_iterator iter = invariantBases.begin();
              iter != invariantBases.end(); ++iter ){
            numberInvariant += (*iter).second;
        }
        * propInvariant = (double)numberInvariant/
                          (double)(numberInvariant + sequences.numberColumns());
    }
	
    return rates;
}

vector < double > RnaModel::retrieveEmpiricalFrequencies( SequenceTable * sequenceTable, int modelId, vector < double > * countSiteParam) {
    vector < double > countSite( 64 );
    vector < double > countNuc( 16 );
    vector < double > freqState( 16 );

    for ( int i = 0; i < 64; ++i ) {
        countSite[i] = 0.0;
    }

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair< string, unsigned int > > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s;

    for ( unsigned int sequence = 0; sequence < sequenceTable->getNumberSpecies(); ++sequence ) {
        for ( unsigned int site = 0; site < sequenceTable->getSequencesLength(modelId); ++site ) {
            s = getSymbolNumber( sequences( sequence, site ) );
            if (s == getNumberSymbols()){
               errorSymbol( sequenceTable, modelId, sequence, site );
            }
			countSite[s]++;
		}
	}
	
    for ( unsigned int site = 0; site < invariantBases.size(); ++site ){
        s = getSymbolNumber( invariantBases[site].first );
        if (s == getNumberSymbols()){
            errorSymbol( sequenceTable, modelId, 0, -site-1 );
        }
        countSite[s] += ( invariantBases[site].second * sequenceTable->getNumberSpecies() );
	}
	
    for ( int nucleotide = 0; nucleotide < 16; ++nucleotide ) {
		countNuc[nucleotide] = countSite[nucleotide];
	}
	
	// Exclude double gaps from the frequency calculations (character 63 = '--').
	// To exclude single gaps, end at 49 rather than 63.
	for ( int ambig = 16; ambig < 49; ++ambig ) {
		if (countSite[ambig] > 0) {
    		vector < unsigned int > ambigBases = getAggregateStates( ambig );
			
			double sum = 0.0;
			for ( int i = 0; i < ambigBases.size(); ++i ) {
				sum += countNuc[ambigBases[i]];
			}
			
			for ( int i = 0; i < ambigBases.size(); ++i ) {
				for ( int state = 0; state < 16; ++state ) {
					if ( equivalencyTable(ambigBases[i], state) == 1 ) {
						countSite[state] += countSite[ambig] * countNuc[ambigBases[i]] / sum;
					}
				}
			}
		}
	}
	
    double sum = 0.0;
    for ( int nucleotide = 0; nucleotide < 16; ++nucleotide ) {
        sum +=  countSite[nucleotide];
    }
	
    for ( int nucleotide = 0; nucleotide < 16; ++nucleotide ) {
    	for ( int state = 0; state < 16; ++state ) {
			if ( equivalencyTable(nucleotide, state) == 1 ) {
				// Use a very small pseudocount value, to avoid frequencies of zero.
    			freqState[state] += (countSite[nucleotide] / sum)  + 1e-6;
			}
		}
    }
	
    sum = 0.0;
    for ( int state = 0; state < 16; ++state ) {
        sum +=  freqState[state];
    }
    for ( int state = 0; state < 16; ++state ) {
		// Ensure that frequencies add up to 1, having added a pseudocount.
    	freqState[state] /= sum;
    }
	
    if ( countSiteParam != NULL ) {
		* countSiteParam = countSite;
	}
	
    return freqState;
}

int RnaModel::getSymbolNumber( const string & pairedBase, unsigned int ) const {
    switch ( pairedBase[0] ) {
        case 'a': case 'A':
            switch ( pairedBase[1] ) {
                case 'a': case 'A':
                    return 6;
                case 'g': case 'G':
                    return 7;
                case 'c': case 'C':
                    return 8;
                case 'u': case 'U': case 't': case 'T':
                    return 0;
                case 'x': case 'X': case 'N': case 'n': case '?':
                    return 16;
                case 'r': case 'R':
                    return 17;
                case 'y': case 'Y':
                    return 18;
                case '-':
                    return 56;
                default:
                    return 64;
            }
        case 'c': case 'C':
            switch ( pairedBase[1] ) {
                case 'a':  case 'A':
                    return 11;
                case 'g': case 'G':
                    return 5;
                case 'c': case 'C':
                    return 12;
                case 'u': case 'U': case 't': case 'T':
                    return 13;
                case 'x': case 'X': case 'N': case 'n': case '?':
                    return 19;
                case 'r': case 'R':
                    return 20;
                case 'y': case 'Y':
                    return 21;
                case '-':
                    return 58;
                default:
                    return 64;
            }
        case 'G': case 'g':
            switch ( pairedBase[1] ) {
                case 'a': case 'A':
                    return 9;
                case 'g': case 'G':
                    return 10;
                case 'c': case 'C':
                    return 2;
                case 'u': case 'U': case 't': case 'T':
                    return 1;
                case 'x': case 'X': case 'N': case 'n': case '?':
                    return 22;
                case 'r': case 'R':
                    return 23;
                case 'y': case 'Y':
                    return 24;
                case '-':
                    return 57;
                default:
                    return 64;
            }
        case 'u': case 'U': case 't': case 'T':
            switch ( pairedBase[1] ) {
                case 'a': case 'A':
                    return 3;
                case 'g': case 'G':
                    return 4;
                case 'c': case 'C':
                    return 14;
                case 'u': case 'U': case 't': case 'T':
                    return 15;
                case 'x': case 'X': case 'N': case 'n': case '?':
                    return 25;
                case 'r': case 'R':
                    return 26;
                case 'y': case 'Y':
                    return 27;
                case '-':
                    return 59;
                default:
                    return 64;
            }
        case 'x': case 'X': case 'n': case 'N': case '?':
            switch ( pairedBase[1] ) {
                case 'a': case 'A':
                    return 28;
                case 'g': case 'G':
                    return 34;
                case 'c': case 'C':
                    return 31;
                case 'u': case 'U': case 't': case 'T':
                    return 37;
                case 'x': case 'X': case 'n': case 'N': case '?':
                    return 40;
                case 'r': case 'R':
                    return 41;
                case 'y': case 'Y':
                    return 42;
                case '-':
                    return 60;
                default:
                    return 64;
            }
        case 'r': case 'R':
            switch ( pairedBase[1] ) {
                case 'a': case 'A':
                    return 29;
                case 'g': case 'G':
                    return 35;
                case 'c': case 'C':
                    return 32;
                case 'u': case 'U': case 't': case 'T':
                    return 38;
                case 'x': case 'X': case 'n': case 'N': case '?':
                    return 46;
                case 'r': case 'R':
                    return 47;
                case 'y': case 'Y':
                    return 48;
                case '-':
                    return 61;
                default:
                    return 64;
            }
        case 'y': case 'Y':
            switch ( pairedBase[1] ) {
                case 'a': case 'A':
                    return 30;
                case 'c': case 'C':
                    return 33;
                case 'g': case 'G':
                    return 36;
                case 'u': case 'U': case 't': case 'T':
                    return 39;
                case 'x': case 'X': case 'n': case 'N': case '?':
                    return 43;
                case 'r': case 'R':
                    return 44;
                case 'y': case 'Y':
                    return 45;
                case '-':
                    return 62;
                default:
                    return 64;
            }
        case '-':
            switch ( pairedBase[1] ) {
                case 'a': case 'A':
                    return 49;
                case 'g': case 'G':
                    return 50;
                case 'c': case 'C':
                    return 51;
                case 'u': case 'U': case 't': case 'T':
                    return 52;
                case 'x': case 'X': case 'n': case 'N': case '?':
                    return 53;
                case 'r': case 'R':
                    return 54;
                case 'y': case 'Y':
                    return 55;
                case '-':
                    return 63;
                default:
                    return 64;
            }
        default :
            return 64;
    }
}

string RnaModel::getSymbol( unsigned int symbolNumber, unsigned int ) const {
    switch ( symbolNumber ) {
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
            // Mismatches
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

            // Ambiguous base-pairs
        case 16 :
            return ( string( "AX" ) );
        case 17 :
            return ( string( "AR" ) );
        case 18 :
            return ( string( "AY" ) );
        case 19 :
            return ( string( "CX" ) );
        case 20 :
            return ( string( "CR" ) );
        case 21 :
            return ( string( "CY" ) );
        case 22 :
            return ( string( "GX" ) );
        case 23 :
            return ( string( "GR" ) );
        case 24 :
            return ( string( "GY" ) );
        case 25 :
            return ( string( "UX" ) );
        case 26 :
            return ( string( "UR" ) );
        case 27 :
            return ( string( "UY" ) );
        case 28 :
            return ( string( "XA" ) );
        case 29 :
            return ( string( "RA" ) );
        case 30 :
            return ( string( "YA" ) );
        case 31 :
            return ( string( "XC" ) );
        case 32 :
            return ( string( "RC" ) );
        case 33 :
            return ( string( "YC" ) );
        case 34 :
            return ( string( "XG" ) );
        case 35 :
            return ( string( "RG" ) );
        case 36 :
            return ( string( "YG" ) );
        case 37 :
            return ( string( "XU" ) );
        case 38 :
            return ( string( "RU" ) );
        case 39 :
            return ( string( "YU" ) );
        case 40 :
            return ( string( "XX" ) );
        case 41 :
            return ( string( "XR" ) );
        case 42 :
            return ( string( "XY" ) );
        case 43 :
            return ( string( "YX" ) );
        case 44 :
            return ( string( "YR" ) );
        case 45 :
            return ( string( "YY" ) );
        case 46 :
            return ( string( "RX" ) );
        case 47 :
            return ( string( "RR" ) );
        case 48 :
            return ( string( "RY" ) );

            // "gap" base-pair
        case 49 :
            return ( string( "-A" ) );
        case 50 :
            return ( string( "-G" ) );
        case 51 :
            return ( string( "-C" ) );
        case 52 :
            return ( string( "-U" ) );
        case 53 :
            return ( string( "-X" ) );
        case 54 :
            return ( string( "-R" ) );
        case 55 :
            return ( string( "-Y" ) );
        case 56 :
            return ( string( "A-" ) );
        case 57 :
            return ( string( "G-" ) );
        case 58 :
            return ( string( "C-" ) );
        case 59 :
            return ( string( "U-" ) );
        case 60 :
            return ( string( "X-" ) );
        case 61 :
            return ( string( "R-" ) );
        case 62 :
            return ( string( "Y-" ) );
        case 63 :
            return ( string( "--" ) );
        default :
            return ( string( "ZZ" ) ); // unknown base-pair
    }
}

vector < unsigned int > RnaModel::getAggregateStates( unsigned int symbolNumber ) {
	vector < unsigned int > aggregateStates;
	
    switch ( symbolNumber ) {
        case 0 :
			aggregateStates.push_back(0); break;
        case 1 :
			aggregateStates.push_back(1); break;
        case 2 :
			aggregateStates.push_back(2); break;
        case 3 :
			aggregateStates.push_back(3); break;
        case 4 :
			aggregateStates.push_back(4); break;
        case 5 :
			aggregateStates.push_back(5); break;
        case 6 :
			aggregateStates.push_back(6); break;
        case 7 :
			aggregateStates.push_back(7); break;
        case 8 :
			aggregateStates.push_back(8); break;
        case 9 :
			aggregateStates.push_back(9); break;
        case 10 :
			aggregateStates.push_back(10); break;
        case 11 :
			aggregateStates.push_back(11); break;
        case 12 :
			aggregateStates.push_back(12); break;
        case 13 :
			aggregateStates.push_back(13); break;
        case 14 :
			aggregateStates.push_back(14); break;
        case 15 :
			aggregateStates.push_back(15); break;
        case 16 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(8); break;
        case 17 :
			aggregateStates.push_back(6);
			aggregateStates.push_back(7); break;
        case 18 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(8); break;
        case 19 :
			aggregateStates.push_back(5);
			aggregateStates.push_back(11);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13); break;
        case 20 :
			aggregateStates.push_back(5);
			aggregateStates.push_back(11); break;
        case 21 :
			aggregateStates.push_back(12);
			aggregateStates.push_back(13); break;
        case 22 :
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10); break;
        case 23 :
			aggregateStates.push_back(9);
			aggregateStates.push_back(10); break;
        case 24 :
			aggregateStates.push_back(1);
			aggregateStates.push_back(2); break;
        case 25 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 26 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(4); break;
        case 27 :
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 28 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(6);
			aggregateStates.push_back(9);
			aggregateStates.push_back(11); break;
        case 29 :
			aggregateStates.push_back(6);
			aggregateStates.push_back(9); break;
        case 30 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(11); break;
        case 31 :
			aggregateStates.push_back(2);
			aggregateStates.push_back(8);
			aggregateStates.push_back(12);
			aggregateStates.push_back(14); break;
        case 32 :
			aggregateStates.push_back(2);
			aggregateStates.push_back(8); break;
        case 33 :
			aggregateStates.push_back(12);
			aggregateStates.push_back(14); break;
        case 34 :
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(7);
			aggregateStates.push_back(10); break;
        case 35 :
			aggregateStates.push_back(7);
			aggregateStates.push_back(10); break;
        case 36 :
			aggregateStates.push_back(4);
			aggregateStates.push_back(5); break;
        case 37 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(13);
			aggregateStates.push_back(15); break;
        case 38 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1); break;
        case 39 :
			aggregateStates.push_back(13);
			aggregateStates.push_back(15); break;
        case 40 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(8);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10);
			aggregateStates.push_back(11);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 41 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10);
			aggregateStates.push_back(11); break;
        case 42 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(8);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 43 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(11);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 44 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(11); break;
        case 45 :
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 46 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(8);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10); break;
        case 47 :
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10); break;
        case 48 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(8); break;
        case 49 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(6);
			aggregateStates.push_back(9);
			aggregateStates.push_back(11); break;
        case 50 :
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(7);
			aggregateStates.push_back(10); break;
        case 51 :
			aggregateStates.push_back(2);
			aggregateStates.push_back(8);
			aggregateStates.push_back(12);
			aggregateStates.push_back(14); break;
        case 52 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(13);
			aggregateStates.push_back(15); break;
        case 53 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(8);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10);
			aggregateStates.push_back(11);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 54 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10);
			aggregateStates.push_back(11); break;
        case 55 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(8);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 56 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(8); break;
        case 57 :
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10); break;
        case 58 :
			aggregateStates.push_back(5);
			aggregateStates.push_back(11);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13); break;
        case 59 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 60 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(8);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10);
			aggregateStates.push_back(11);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 61 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(8);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10); break;
        case 62 :
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(11);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
        case 63 :
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(3);
			aggregateStates.push_back(4);
			aggregateStates.push_back(5);
			aggregateStates.push_back(6);
			aggregateStates.push_back(7);
			aggregateStates.push_back(8);
			aggregateStates.push_back(9);
			aggregateStates.push_back(10);
			aggregateStates.push_back(11);
			aggregateStates.push_back(12);
			aggregateStates.push_back(13);
			aggregateStates.push_back(14);
			aggregateStates.push_back(15); break;
	}
	
	return aggregateStates;
}

double RnaModel::probability( unsigned int obase, unsigned int sbase,
                 double time, unsigned int category, unsigned int )  const {
    if ( ( invariantCategory == 1 ) && ( category == 0 ) ) {
        if ( obase == sbase ) {
            return 1.0;
        }
        else {
            return 0.0;
        }
    }

    int cat = category - invariantCategory;
    double probability = 0.0;
    for ( unsigned int k = 0; k < matrixSize; ++k )
        probability += ( ieigenMatrix[cat] ( obase, k ) *
        eigenMatrix[cat] ( k, sbase ) * exp( eigenValues[cat] [k] * time ) );

    // Just make sure
    //assert ( ( probability >= -0.000001 ) && ( probability <= 1.000001 ) );
    probability = MIN( probability, 1.0 );
    probability = MAX( probability, 0.0 );

    return ( probability );
}


double RnaModel::diffProbability( unsigned int obase, unsigned int sbase,
                double time, unsigned int category, unsigned int )  const {
    if ( ( invariantCategory == 1 ) && ( category == 0 ) ) {
        return 0.0;
    }
    int cat = category - invariantCategory;
    double diffProbability = 0.0;
    for ( unsigned int k = 0; k < matrixSize; ++k ) {
        diffProbability += ( eigenValues[cat] [k] *
        ieigenMatrix[cat] ( obase, k ) * eigenMatrix[cat] ( k, sbase ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    assert( !isnan( diffProbability ) );
    return ( diffProbability );
}



double RnaModel::secondDiffProbability( unsigned int obase, unsigned int sbase,
     double time, unsigned int category, unsigned int )  const {
    if ( ( invariantCategory == 1 ) && ( category == 0 ) ){
        return 0.0;
    }
    int cat = category - invariantCategory;
    double diff2Probability = 0.0;
    for ( unsigned int k = 0; k < matrixSize; ++k ) {
        diff2Probability += ( pow( eigenValues[cat] [k], 2.0 ) *
            ieigenMatrix[cat] ( obase, k ) * eigenMatrix[cat] ( k, sbase ) *
            exp( eigenValues[cat] [k] * time ) );
    }
    assert( !isnan( diff2Probability ) );
    return ( diff2Probability );
}


void RnaModel::printParameters( ostream & outputStream ) const {

    //print common parameters
    MatrixModel::printParameters( outputStream );

    //rely on the general method to print matrices
    MatrixModel::printParametersAux( outputStream );

}

void RnaModel::setEigenMatrix() {
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
          ++category ) {

        // Initialise the rate matrix
        for ( unsigned int i = 0; i < matrixSize - 1; ++i ) {
            for ( unsigned int j = i + 1; j < matrixSize; ++j ) {
                rateMatrix[category]( i, j ) =
                        getExchangeability( i, j, category ) *
                        getFrequency( j, category + invariantCategory) *
                        substitutionRate[category];
                rateMatrix[category] ( j, i ) =
                        getExchangeability( j, i, category ) *
                        getFrequency( i, category + invariantCategory) *
                        substitutionRate[category];
#ifdef DEBUG
                assert( !isnan(rateMatrix[category]( i, j )) );
                assert( !isnan(rateMatrix[category]( j, i )) );
#endif
            }
        }
        for ( unsigned int i = 0; i < matrixSize; ++i ) {
            double sum = 0.0;
            for ( unsigned int j = 0; j < matrixSize; ++j ) {
                if ( i != j ){
                    sum += rateMatrix[category] ( i, j );
                }
            }
            rateMatrix[category] ( i, i ) = -sum;
        }


#ifdef DEBUG4
        for ( unsigned int i = 0; i < matrixSize; ++i ) {
            for ( unsigned int j = 0; j < matrixSize; ++j ) {
                cout << setw(6) << setprecision(3) << rateMatrix[category] ( i, j ) << " ";
            }
            cout << endl;
        }
#endif

        vector < double > junk( matrixSize );
        for ( unsigned int i = 0; i < matrixSize; ++i ) {
            eigenValues[category] [i] = 0.0;
            junk[i] = 0.0;
            for ( unsigned int j = 0; j < matrixSize; ++j ){
                eigenMatrix[category] ( i, j ) = 0.0;
                ieigenMatrix[category] ( i, j ) = 0.0;
            }
        }
        // Find the eigenvalues and eigenvectors of the rate matrix
        int err = LeftEigenSystem( rateMatrix[category], eigenMatrix[category],
        eigenValues[category], junk );
        if (err){
            cerr << "Error during matrix exponentiation" << endl;
            exit(EXIT_FAILURE);
        }
        // Find the inverse of the matrix whose rows are the left eigen
        // vectors of the rate matrix
        ieigenMatrix[category] = inverse( eigenMatrix[category] );
    }
}
