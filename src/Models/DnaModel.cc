#include "Models/DnaModel.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "Sequence/SequenceTable.h"

#include "Models/Perturbator.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/statlib.h"
#include "Util/array2D.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"


DnaModel::DnaModel( const string & registrationName ) :
        MatrixModel( registrationName ){
}  

DnaModel::DnaModel( ParametersSet & parameters ) : MatrixModel( parameters ){
    initMatrix( 4 );
}

DnaModel::~DnaModel() {
}

string DnaModel::getState( unsigned int stateNumber, unsigned int ) const {
    switch ( stateNumber ) {
        // Watson-Crick base pairs
        case 0 :
            return ( string( "A" ) );
        case 1 :
            return ( string( "C" ) );
        case 2 :
            return ( string( "G" ) );
        case 3 :
            return ( string( "T" ) );
        default :
            return("");
    }
}
          
int DnaModel::getSymbolNumber( const string & base, unsigned int ) const {
    if (base.length() != 1)
        return 16;
    switch ( base[0] ) {
        case 'a': case 'A':
            return 0;
        case 'c': case 'C':
            return 1;
        case 'g': case 'G':
            return 2;
        case 'u': case 't': case 'U': case 'T':
            return 3;
        case 'r': case 'R':
            return 4; // Unknown purine       (A or G)
        case 'y': case 'Y':
            return 5; // Unknown pyramidine   (C or U)
        case 'm': case 'M':
            return 6; //                      (A or C)
        case 'k': case 'K':
            return 7; //                      (U or G)
        case 'w': case 'W':
            return 8; //                      (U or A)
        case 's': case 'S':
            return 9; //                      (C or G)
        case 'b': case 'B':
            return 10; //                     (U, C or G)
        case 'd': case 'D':
            return 11; //                     (U, A or G)
        case 'h': case 'H':
            return 12; //                     (U, A or C)
        case 'v': case 'V':
            return 13; //                     (A, C or G)
        case 'n': case 'x': case '?': case 'N': case 'X':
            return 14; // Unknown nucleotide  (A, C, G or T)
        case '-' :
            return 15; // gaps are treated as unknown nucleotides
        default :
            return 16; // Out of range (unrecognized symbol)
    }
}



vector < double > DnaModel::retrieveEmpiricalRates( SequenceTable * sequenceTable,
int modelId, double * propInvariant ) {

    array2D < double > countRates( 16, 16 );
    vector < double > rates5( 5 );
    double reference;
    
    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair<string, unsigned int> > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s1, s2;

    for ( unsigned int i = 0 ; i < 16; ++i ) {
        for ( unsigned int j = i+1 ; j < 16; ++j ) {
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
	
	// Exclude gaps from the rate calculations (character 15 = '-').
	for ( int s1 = 0; s1 < 15; ++s1 ) {
		for ( int s2 = s1+1; s2 < 15; ++s2 ) {
			if (!(s1 < 4 && s2 < 4) && s1 != s2) {
				if (countRates( s1, s2 ) > 0) {
					double sum = 0.0;
					for ( int nuc1 = 0; nuc1 < 4; ++nuc1 ) {
						for ( int nuc2 = 0; nuc2 < 4; ++nuc2 ) {
							if ( ( equivalencyTable(s1, nuc1) == 1 ) &&
								 ( equivalencyTable(s2, nuc2) == 1 ) ) {
								sum++;
							}
						}
					}
					for ( int nuc1 = 0; nuc1 < 4; ++nuc1 ) {
						for ( int nuc2 = 0; nuc2 < 4; ++nuc2 ) {
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
	for ( int nuc1 = 0; nuc1 < 4; ++nuc1 ) {
		for ( int nuc2 = nuc1+1; nuc2 < 4; ++nuc2 ) {
			// Use a pseudocount to ensure that we don't start with zero-valued
			// rates, because the optimizer can't then escape from them.
			countRates( nuc1, nuc2 ) += 1;
			sum += countRates( nuc1, nuc2 );
		}
	}
	
	reference = countRates( 0, 2 ) / sum; // A-G is the reference
	rates5[0] = countRates( 0, 1 ) / (sum * reference);
	rates5[1] = countRates( 0, 3 ) / (sum * reference);
	rates5[2] = countRates( 1, 2 ) / (sum * reference);
	rates5[3] = countRates( 1, 3 ) / (sum * reference);
	rates5[4] = countRates( 2, 3 ) / (sum * reference);
	
    if ( propInvariant != NULL ) {
        int numberInvariant = 0;
        for ( vector< pair< string, unsigned int > >::const_iterator iter = invariantBases.begin();
              iter != invariantBases.end(); ++iter ){
            numberInvariant += (*iter).second;
        }
        * propInvariant = (double)numberInvariant/
                          (double)(numberInvariant + sequences.numberColumns());
    }
	
    return rates5;
}

vector < double > DnaModel::retrieveEmpiricalFrequencies
( SequenceTable * sequenceTable, int modelId ) {

    vector < double > countSite( 16 );
    vector < double > countNuc( 4 );
    vector < double > freqState( 4 );

    for ( int i = 0; i < 16; ++i ) {
        countSite[i] = 0.0;
    }

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair< string, unsigned int > > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s;
    
    for ( unsigned int sequence = 0; sequence < sequenceTable->getNumberSpecies(); ++sequence ){
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
	
    for ( int i = 0; i < 4; ++i ) {
        countNuc[i] = countSite[i];
    }
	
	// Exclude gaps from the frequency calculations (character 15 = '-').
	for ( int ambig = 4; ambig < 15; ++ambig ) {
		if (countSite[ambig] > 0) {
			double sum = 0.0;
			for ( int nucleotide = 0; nucleotide < 4; ++nucleotide ) {
				if ( equivalencyTable(ambig, nucleotide) == 1 ) {
					sum += countNuc[nucleotide];
				}
			}
			for ( int nucleotide = 0; nucleotide < 4; ++nucleotide ) {
				if ( equivalencyTable(ambig, nucleotide) == 1 ) {
					countSite[nucleotide] += countSite[ambig] * countNuc[nucleotide] / sum;
				}
			}
		}
	}
    
    double sum = 0.0;
    for ( int nucleotide = 0; nucleotide < 4; ++nucleotide ) {
		sum += countSite[nucleotide];
	}
    for ( int nucleotide = 0; nucleotide < 4; ++nucleotide ) {
		// Use a very small pseudocount value, to avoid frequencies of zero.
        freqState[nucleotide] = (countSite[nucleotide] / sum)  + 1e-6;
    }
	
    sum = 0.0;
    for ( int state = 0; state < 4; ++state ) {
        sum +=  freqState[state];
    }
    for ( int state = 0; state < 4; ++state ) {
		// Ensure that frequencies add up to 1, having added a pseudocount.
    	freqState[state] /= sum;
    }
	
    return freqState;
}

string DnaModel::getSymbol( unsigned int basenumber, unsigned int ) const {
    switch ( basenumber ) {
        case 0 :
            return ( string( "A" ) );
        case 1 :
            return ( string( "C" ) );
        case 2 :
            return ( string( "G" ) );
        case 3 :
            return ( string( "T" ) );
        case 4 :
            return ( string( "R" ) ); // Unknown Purine     (A or G)
        case 5 :
            return ( string( "Y" ) ); // Unknown Pyramidine (C or T)
        case 6 :
            return ( string( "M" ) ); //                    (A or C)
        case 7 :
            return ( string( "K" ) ); //                    (U or G)
        case 8 :
            return ( string( "W" ) ); //                    (U or A)
        case 9 :
            return ( string( "S" ) ); //                    (C or G)
        case 10 :
            return ( string( "B" ) ); //                 (C, G or T)
        case 11 :
            return ( string( "D" ) ); //                 (A, G or T)
        case 12 :
            return ( string( "H" ) ); //                 (A, C or T)
        case 13 :
            return ( string( "V" ) ); //                 (A, C or G)
        case 14 :
            return ( string( "N" ) ); // Unknown    (A , C , G or T)
        case 15 :
            return ( string( "-" ) ); // Gaps
        default :
            return ( string( "Z" ) ); // Unrecognized symbol
    }
}

double DnaModel::probability( unsigned int oldState, unsigned int newState, double time,
unsigned int category, unsigned int ) const {
    int k, cat;
    double probability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) ) {
        if ( oldState == newState )
            return 1.0;
        else
            return 0.0;
    }
    cat = category - invariantCategory;
    for ( k = 0; k < 4; k++ ) {
        probability += ( ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    //just to be sure
    //assert ( ( probability >= -0.000001 ) && ( probability <= 1.000001 ) );
    probability = MIN( probability, 1.0 );
    probability = MAX( probability, 0.0 );
    return ( probability );
}


double DnaModel::diffProbability( unsigned int oldState, unsigned int newState, double time,
unsigned int category, unsigned int ) const {
    int k, cat;
    double diffProbability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) )
        return 0.0;
    cat = category - invariantCategory;
    for ( k = 0; k < 4; k++ ) {
        diffProbability += ( eigenValues[cat] [k] *
        ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    return ( diffProbability );
}


double DnaModel::secondDiffProbability( unsigned int oldState, unsigned int newState, double time,
unsigned int category, unsigned int ) const {
    int k, cat;
    double diff2Probability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) )
        return 0.0;
    cat = category - invariantCategory;
    for ( k = 0; k < 4; k++ ) {
        diff2Probability +=
        ( pow( eigenValues[cat] [k], 2.0 ) *
        ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    return ( diff2Probability );
}


void DnaModel::initEquivalencyTable(){
    int i, j;
    equivalencyTable.resize( 16, 4 );
    // Set the patterns for valid states
    for ( i = 0; i < 4; i++ ) {
        for ( j = 0; j < 4; j++ ){
            equivalencyTable( i, j ) = 0.0;
        }
        equivalencyTable( i, i ) = 1.0;
    }
    // Match a Purine
    equivalencyTable( 4, 0 ) = 1.0; // A
    equivalencyTable( 4, 1 ) = 0.0;
    equivalencyTable( 4, 2 ) = 1.0; // G
    equivalencyTable( 4, 3 ) = 0.0;
    // Match a Pyrimidine
    equivalencyTable( 5, 0 ) = 0.0;
    equivalencyTable( 5, 1 ) = 1.0; // C
    equivalencyTable( 5, 2 ) = 0.0;
    equivalencyTable( 5, 3 ) = 1.0; // T/U
    // M match A/C
    equivalencyTable( 6, 0 ) = 1.0; // A
    equivalencyTable( 6, 1 ) = 1.0; // C
    equivalencyTable( 6, 2 ) = 0.0;
    equivalencyTable( 6, 3 ) = 0.0;
    // K match G/T
    equivalencyTable( 7, 0 ) = 0.0;
    equivalencyTable( 7, 1 ) = 0.0;
    equivalencyTable( 7, 2 ) = 1.0; // G
    equivalencyTable( 7, 3 ) = 1.0; // U
    // W match A/T
    equivalencyTable( 8, 0 ) = 1.0; // A
    equivalencyTable( 8, 1 ) = 0.0;
    equivalencyTable( 8, 2 ) = 0.0;
    equivalencyTable( 8, 3 ) = 1.0; // T/U
    // S match C/G
    equivalencyTable( 9, 0 ) = 0.0; 
    equivalencyTable( 9, 1 ) = 1.0; // C
    equivalencyTable( 9, 2 ) = 1.0; // G
    equivalencyTable( 9, 3 ) = 0.0;
    // B match C/G/U
    equivalencyTable( 10, 0 ) = 0.0;
    equivalencyTable( 10, 1 ) = 1.0; // C
    equivalencyTable( 10, 2 ) = 1.0; // G
    equivalencyTable( 10, 3 ) = 1.0; // T/U
    // D match A/G/U
    equivalencyTable( 11, 0 ) = 1.0; // A
    equivalencyTable( 11, 1 ) = 0.0;
    equivalencyTable( 11, 2 ) = 1.0; // G
    equivalencyTable( 11, 3 ) = 1.0; // T/U
    // H match A/C/U
    equivalencyTable( 12, 0 ) = 1.0; // A
    equivalencyTable( 12, 1 ) = 1.0; // C
    equivalencyTable( 12, 2 ) = 0.0;
    equivalencyTable( 12, 3 ) = 1.0; // T/U
    // V match A/C/G
    equivalencyTable( 13, 0 ) = 1.0; // A
    equivalencyTable( 13, 1 ) = 1.0; // C
    equivalencyTable( 13, 2 ) = 1.0; // G
    equivalencyTable( 13, 3 ) = 0.0;
    // Match all nucleotides (unknown)
    equivalencyTable( 14, 0 ) = 1.0; // A
    equivalencyTable( 14, 1 ) = 1.0; // C
    equivalencyTable( 14, 2 ) = 1.0; // G
    equivalencyTable( 14, 3 ) = 1.0; // T/U
    // Gap (No information)
    equivalencyTable( 15, 0 ) = 1.0;
    equivalencyTable( 15, 1 ) = 1.0;
    equivalencyTable( 15, 2 ) = 1.0;
    equivalencyTable( 15, 3 ) = 1.0;
}





// print the parameters in readable format
void DnaModel::printParameters( ostream & outputStream ) const {
  
    //WARNING, redefined in many children GG, ..
  
    //print common parameters
    MatrixModel::printParameters( outputStream );
    
    //rely on the general method to print matrices
    MatrixModel::printParametersAux( outputStream );
}
