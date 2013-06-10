#include "Models/ThreeStateModel.h"

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


ThreeStateModel::ThreeStateModel( const string & registrationName ) :
        MatrixModel( registrationName ){
}  

ThreeStateModel::ThreeStateModel( ParametersSet & parameters ) : MatrixModel( parameters ){
    if ( discreteGamma ){
        eigenMatrix = new array2D < double > [numberGammaCategories];
        ieigenMatrix = new array2D < double > [numberGammaCategories];
        eigenValues = new vector < double > [numberGammaCategories];
        rateMatrix = new array2D < double > [numberGammaCategories];
        for (unsigned int i = 0; i < numberGammaCategories; ++i){
            eigenMatrix[i].resize(3,3);
            ieigenMatrix[i].resize(3,3);
            eigenValues[i].resize(3);
            rateMatrix[i].resize(3,3);
        }
    }
    else{
        eigenMatrix = new array2D < double >;
        ieigenMatrix = new array2D < double >;
        eigenValues = new vector < double >;
        rateMatrix = new array2D < double >;
        eigenMatrix->resize(3,3);
        ieigenMatrix->resize(3,3);
        eigenValues->resize(3);
        rateMatrix->resize(3,3);
    }
}

ThreeStateModel::~ThreeStateModel() {
    if ( discreteGamma ){
        if ( eigenMatrix != NULL ){
            delete[] eigenMatrix;
            eigenMatrix = NULL;
        }
        if ( ieigenMatrix != NULL ){
            delete[] ieigenMatrix;
            ieigenMatrix = NULL;
        }
        if ( eigenValues != NULL ){
            delete[] eigenValues;
            eigenValues = NULL;
        }
        if ( rateMatrix != NULL ){
            delete[] rateMatrix;
            rateMatrix = NULL;
        }
    }
    else{
        if ( eigenMatrix != NULL ){
            delete eigenMatrix;
            eigenMatrix = NULL;
        }
        if ( ieigenMatrix != NULL ){
            delete ieigenMatrix;
            ieigenMatrix = NULL;
        }
        if ( eigenValues != NULL ){
            delete eigenValues;
            eigenValues = NULL;
        }
        if ( rateMatrix != NULL ){
            delete rateMatrix;
            rateMatrix = NULL;
        }
    }
}

string ThreeStateModel::getState( unsigned int stateNumber, unsigned int ) const {
    switch ( stateNumber ) {
        // Watson-Crick base pairs
        case 0 :
            return ( string( "A" ) );
        case 1 :
            return ( string( "Y" ) );
        case 2 :
            return ( string( "G" ) );
        default :
            return("");
    }
}
          
int ThreeStateModel::getSymbolNumber( const string & base, unsigned int ) const {
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
            return 5; // Unknown pyrimidine   (C or U)
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
            return 12; //                     (A or C)
        case 'v': case 'V':
            return 13; //                     (A or C)
        case 'n': case 'x': case '?': case 'N': case 'X':
            return 14; // Unknown nucleotide  (A, C, G or T)
        case '-' :
            return 15; // gaps are treated as unknown nucleotides
        default :
            return 16; // Out of range (unrecognized symbol)
    }
}


vector < double > ThreeStateModel::retrieveEmpiricalFrequencies
( SequenceTable * sequenceTable, int modelId ) {

    vector < double > countSite( 16 );
    vector < double > countNuc( 4 );
    vector < double > freqState( 3 );

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
	
    for ( int nucleotide = 0; nucleotide < 4; ++nucleotide ) {
		countNuc[nucleotide] = countSite[nucleotide];
	}
	
	// Exclude gaps from the frequency calculations (character 15 = '-').
	for ( int ambig = 4; ambig < 15; ++ambig ) {
		if (countSite[ambig] > 0) {
			// This can get a little tricky if an ambiguity character or a gap
			// matches a 'Y' site, since we may need to make the count
			// proportional to the relevant C or T count, not the combined count.
    		vector < unsigned int > ambigBases = getAggregateStates( ambig );
			
            double sum = 0.0;
			for ( int i = 0; i < ambigBases.size(); ++i ) {
				sum += countNuc[ambigBases[i]];
			}
			
			for ( int i = 0; i < ambigBases.size(); ++i ) {
				for ( int state = 0; state < 3; ++state ) {
					if ( equivalencyTable(ambigBases[i], state) == 1 ) {
						countSite[state] += countSite[ambig] * countNuc[ambigBases[i]] / sum;
					}
				}
			}
		}
	}
    
    double sum = 0.0;
    for ( int nucleotide = 0; nucleotide < 4; ++nucleotide ) {
		sum += countSite[nucleotide];
	}
    for ( int nucleotide = 0; nucleotide < 4; ++nucleotide ) {
		for ( int state = 0; state < 3; ++state ) {
			if ( equivalencyTable(nucleotide, state) == 1 ) {
				// Use a very small pseudocount value, to avoid frequencies of zero.
    			freqState[state] += (countSite[nucleotide] / sum)  + 1e-6;
			}
		}
	}
	
    sum = 0.0;
    for ( int state = 0; state < 3; ++state ) {
        sum +=  freqState[state];
    }
    for ( int state = 0; state < 3; ++state ) {
		// Ensure that frequencies add up to 1, having added a pseudocount.
    	freqState[state] /= sum;
    }
	
    return freqState;
}


vector < double > ThreeStateModel::retrieveEmpiricalRates( SequenceTable * sequenceTable,
int modelId, double * propInvariant ) {

    double reference;
    vector < double > rate2( 2 );
    
    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair<string, unsigned int> > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s1, s2;

    reference = 1e-6; // A<->G
    rate2[0] = 1e-6; // A<->C/T
    rate2[1] = 1e-6; // G<->C/T
        
    for ( int site = -((int)invariantBases.size());
            site < (int)(sequences.numberColumns()); ++site ) {
        for ( unsigned int seq = 0;
              seq < sequenceTable->getNumberSpecies()-1; ++seq ) {
            if ( site >= 0 ){
                s1 = getSymbolNumber( sequences( seq, site ) );
            }
            else{
                s1 = getSymbolNumber(
                        invariantBases[-1-site].first );
            }
            if (s1==getNumberSymbols()){
                errorSymbol( sequenceTable, modelId, seq, site );
            }
            for ( unsigned int seq2 = seq + 1;
                  seq2 < sequenceTable->getNumberSpecies(); ++seq2 ) {
                if ( site >= 0 ){
                    s2 = getSymbolNumber( sequences( seq2, site ) );
                }
                else{
                    s2 = s1;
                }
                if (s2 == getNumberSymbols()){
                    errorSymbol( sequenceTable, modelId, seq2, site );
                }
                double sum = 0.0;
                for ( unsigned int state1 = 0; state1 < 4; ++state1 ){
                    for ( unsigned int state2 = 0; state2 < 4; ++state2 ){
                        if (state2!=state1){
                            if ( ( equivalencyTable(s1,state1) == 1 ) &&
                                 ( equivalencyTable(s2,state2) == 1 ) ){
                                sum += 1.0;
                            }
                        }
                    }
                }
                if ( site < 0 ){
                    sum /= (double)invariantBases[-1-site].second;
                }
                for ( unsigned int state1 = 0; state1 < 4; ++state1 ){
                    for ( unsigned int state2 = 0; state2 < 4; ++state2 ){
                        if (state2!=state1){
                            if ( ( equivalencyTable(s1,state1) == 1 ) &&
                                 ( equivalencyTable(s2,state2) == 1 ) ){
                                unsigned int st1, st2;
                                if (state1<state2){
                                    st1 = state1;
                                    st2 = state2;
                                }
                                else{
                                    st1 = state2;
                                    st2 = state1;
                                }
                                if ( ((st1==0)&&(st2==1)) ||
                                     ((st1==0)&&(st2==3)) )
                                    rate2[0] += 1.0/sum;
                                if ( (st1==0)&&(st2==2) )
                                    reference += 1.0/sum;
                                if ( ((st1==1)&&(st2==2)) ||
                                     ((st1==2)&&(st2==3)) )
                                    rate2[1] += 1.0/sum;
                            }
                        }
                    }
                }
            }
        }
    }
    // The rate A<->G is the reference rate, 0 and 1 are slow transversion rates
    rate2[0] /= reference; // A<->C/U
    rate2[1] /= reference; // G<->C/U
    if ( propInvariant != NULL ) {
        int numberInvariant = 0;
        for ( vector< pair< string, unsigned int > >::const_iterator iter = invariantBases.begin();
              iter != invariantBases.end(); ++iter ){
            numberInvariant += (*iter).second;
        }
        * propInvariant = (double)numberInvariant/
                          (double)(numberInvariant + sequences.numberColumns());
    }
    return rate2;
}


string ThreeStateModel::getSymbol( unsigned int basenumber, unsigned int ) const {
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


double ThreeStateModel::probability( unsigned int oldState, unsigned int newState,
                 double time, unsigned int category, unsigned int ) const {
    int k, cat;
    double probability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) ) {
        if ( oldState == newState )
            return 1.0;
        else
            return 0.0;
    }
    cat = category - invariantCategory;
    for ( k = 0; k < 3; k++ ) {
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


double ThreeStateModel::diffProbability( unsigned int oldState, unsigned int newState,
                 double time,unsigned int category, unsigned int ) const {
    int k, cat;
    double diffProbability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) )
        return 0.0;
    cat = category - invariantCategory;
    for ( k = 0; k < 3; k++ ) {
        diffProbability += ( eigenValues[cat] [k] *
        ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    return ( diffProbability );
}


double ThreeStateModel::secondDiffProbability( unsigned int oldState, unsigned int newState,
                double time, unsigned int category, unsigned int ) const {
    int k, cat;
    double diff2Probability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) )
        return 0.0;
    cat = category - invariantCategory;
    for ( k = 0; k < 3; k++ ) {
        diff2Probability +=
        ( pow( eigenValues[cat] [k], 2.0 ) *
        ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    return ( diff2Probability );
}


void ThreeStateModel::initEquivalencyTable(){
    equivalencyTable.resize( 16, 3 );

    // Set the patterns for usual state
    // A
    equivalencyTable( 0, 0 ) = 1.0;
    equivalencyTable( 0, 1 ) = 0.0;
    equivalencyTable( 0, 2 ) = 0.0;
    // C
    equivalencyTable( 1, 0 ) = 0.0;
    equivalencyTable( 1, 1 ) = 1.0;
    equivalencyTable( 1, 2 ) = 0.0;
    // G
    equivalencyTable( 2, 0 ) = 0.0;
    equivalencyTable( 2, 1 ) = 0.0;
    equivalencyTable( 2, 2 ) = 1.0;
    // T/U
    equivalencyTable( 3, 0 ) = 0.0;
    equivalencyTable( 3, 1 ) = 1.0;
    equivalencyTable( 3, 2 ) = 0.0;
    // Match a Purine, A/G
    equivalencyTable( 4, 0 ) = 1.0;
    equivalencyTable( 4, 1 ) = 0.0;
    equivalencyTable( 4, 2 ) = 1.0;
    // Match a Pyramidine, C/U
    equivalencyTable( 5, 0 ) = 0.0;
    equivalencyTable( 5, 1 ) = 1.0;
    equivalencyTable( 5, 2 ) = 0.0;
    // Match a M, A/C
    equivalencyTable( 6, 0 ) = 1.0;
    equivalencyTable( 6, 1 ) = 1.0;
    equivalencyTable( 6, 2 ) = 0.0;
    // Match a K, G/T
    equivalencyTable( 7, 0 ) = 0.0;
    equivalencyTable( 7, 1 ) = 1.0;
    equivalencyTable( 7, 2 ) = 1.0;
    // Match a W, A/T
    equivalencyTable( 8, 0 ) = 1.0;
    equivalencyTable( 8, 1 ) = 1.0;
    equivalencyTable( 8, 2 ) = 0.0;
    // Match a S, C/G
    equivalencyTable( 9, 0 ) = 0.0;
    equivalencyTable( 9, 1 ) = 1.0;
    equivalencyTable( 9, 2 ) = 1.0;
    // Match a B, C/G/U
    equivalencyTable( 10, 0 ) = 0.0;
    equivalencyTable( 10, 1 ) = 1.0;
    equivalencyTable( 10, 2 ) = 1.0;
    // Match a D, A/G/U
    equivalencyTable( 11, 0 ) = 1.0;
    equivalencyTable( 11, 1 ) = 1.0;
    equivalencyTable( 11, 2 ) = 1.0;
    // Match a H, A/C/U
    equivalencyTable( 12, 0 ) = 1.0;
    equivalencyTable( 12, 1 ) = 1.0;
    equivalencyTable( 12, 2 ) = 0.0;
    // Match a V, A/C/G
    equivalencyTable( 13, 0 ) = 1.0;
    equivalencyTable( 13, 1 ) = 1.0;
    equivalencyTable( 13, 2 ) = 1.0;
    // Match all nucleotides (unknown) and gaps
    equivalencyTable( 14, 0 ) = 1.0;
    equivalencyTable( 14, 1 ) = 1.0;
    equivalencyTable( 14, 2 ) = 1.0;
    equivalencyTable( 15, 0 ) = 1.0;
    equivalencyTable( 15, 1 ) = 1.0;
    equivalencyTable( 15, 2 ) = 1.0;
}

vector < unsigned int > ThreeStateModel::getAggregateStates( unsigned int symbolNumber ) {
	vector < unsigned int > aggregateStates;
	
    switch ( symbolNumber ) {
        case 0 : // A
			aggregateStates.push_back(0); break;
        case 1 : // C
			aggregateStates.push_back(1); break;
        case 2 : // G
			aggregateStates.push_back(2); break;
        case 3 : // T
			aggregateStates.push_back(3); break;
        case 4 : // R
			aggregateStates.push_back(0);
			aggregateStates.push_back(2); break;
        case 5 : // Y
			aggregateStates.push_back(1);
			aggregateStates.push_back(3); break;
        case 6 : // M
			aggregateStates.push_back(0);
			aggregateStates.push_back(1); break;
        case 7 : // K
			aggregateStates.push_back(2);
			aggregateStates.push_back(3); break;
        case 8 : // W
			aggregateStates.push_back(0);
			aggregateStates.push_back(3); break;
        case 9 : // S
			aggregateStates.push_back(1);
			aggregateStates.push_back(2); break;
        case 10 : // B
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(3); break;
        case 11 : // D
			aggregateStates.push_back(0);
			aggregateStates.push_back(2);
			aggregateStates.push_back(3); break;
        case 12 : // H
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(3); break;
        case 13 : // V
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2); break;
        case 14 : // N
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(3); break;
        case 15 : // -
			aggregateStates.push_back(0);
			aggregateStates.push_back(1);
			aggregateStates.push_back(2);
			aggregateStates.push_back(3); break;
	}
	return aggregateStates;
}

// print the parameters in readable format
void ThreeStateModel::printParameters( ostream & outputStream ) const {
  
    //print common parameters
    MatrixModel::printParameters( outputStream );

    //rely on the general method to print matrices
    MatrixModel::printParametersAux( outputStream );
}
