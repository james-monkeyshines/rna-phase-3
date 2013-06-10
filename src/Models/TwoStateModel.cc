#include "Models/TwoStateModel.h"

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


TwoStateModel::TwoStateModel( const string & registrationName ) :
        MatrixModel( registrationName ){
}  

TwoStateModel::TwoStateModel( ParametersSet & parameters ) : MatrixModel( parameters ){
    if ( discreteGamma ){
        eigenMatrix = new array2D < double > [numberGammaCategories];
        ieigenMatrix = new array2D < double > [numberGammaCategories];
        eigenValues = new vector < double > [numberGammaCategories];
        rateMatrix = new array2D < double > [numberGammaCategories];
        for (unsigned int i = 0; i < numberGammaCategories; ++i){
            eigenMatrix[i].resize(2,2);
            ieigenMatrix[i].resize(2,2);
            eigenValues[i].resize(2);
            rateMatrix[i].resize(2,2);
        }
    }
    else{
        eigenMatrix = new array2D < double >;
        ieigenMatrix = new array2D < double >;
        eigenValues = new vector < double >;
        rateMatrix = new array2D < double >;
        eigenMatrix->resize(2,2);
        ieigenMatrix->resize(2,2);
        eigenValues->resize(2);
        rateMatrix->resize(2,2);
    }
}

TwoStateModel::~TwoStateModel() {
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

string TwoStateModel::getState( unsigned int stateNumber, unsigned int ) const {
    switch ( stateNumber ) {
        // Watson-Crick base pairs
        case 0 :
            return ( string( "0" ) );
        case 1 :
            return ( string( "1" ) );
        default :
            return("");
    }
}
          
int TwoStateModel::getSymbolNumber( const string & base, unsigned int ) const {
    if (base.length() == 1){
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
                return 12; //                         (A or C)
            case 'v': case 'V':
                return 13; //                     (A or C)
            case 'n': case 'x': case '?': case 'N': case 'X':
                return 14; // Unknown nucleotide  (A, C, G or T)
            case '-' :
                return 15; // gaps are treated as unknown nucleotides
            case '0' :
                return 16; // binary notation 0
            case '1' :
                return 17; // binary notation 0
            default :
                return 82; // Out of range (unrecognized symbol)
        }
    }
    //else, base length == 2
    if ( base.length() != 2 ){
        return 82;
    }
   switch ( base[0] ) {
        case 'a': case 'A':
            switch ( base[1] ) {
                case 'a': case 'A':
                    return 24;
                case 'g': case 'G':
                    return 25;
                case 'c': case 'C':
                    return 26;
                case 'u': case 'U': case 't': case 'T':
                    return 18;
                case 'x': case 'X': case 'N': case 'n': case '?':
                    return 34;
                case 'r': case 'R':
                    return 35;
                case 'y': case 'Y':
                    return 36;
                case '-':
                    return 74;
                default:
                    return 82;
            }
        case 'c': case 'C':
            switch ( base[1] ) {
                case 'a':  case 'A':
                    return 29;
                case 'g': case 'G':
                    return 23;
                case 'c': case 'C':
                    return 30;
                case 'u': case 'U': case 't': case 'T': 
                    return 31;
                case 'x': case 'X': case 'N': case 'n': case '?':
                    return 37;
                case 'r': case 'R':
                    return 38;
                case 'y': case 'Y':
                    return 39;
                case '-':
                    return 76;
                default:
                    return 82;
            }
        case 'G': case 'g':
            switch ( base[1] ) {
                case 'a': case 'A':
                    return 27;
                case 'g': case 'G':
                    return 28;
                case 'c': case 'C':
                    return 20;
                case 'u': case 'U': case 't': case 'T':
                    return 19;
                case 'x': case 'X': case 'N': case 'n': case '?':
                    return 40;
                case 'r': case 'R':
                    return 41;
                case 'y': case 'Y':
                    return 42;
                case '-':
                    return 75;
                default:
                    return 82;
            }
        case 'u': case 'U': case 't': case 'T':
            switch ( base[1] ) {
                case 'a': case 'A':
                    return 21;
                case 'g': case 'G':
                    return 22;
                case 'c': case 'C':
                    return 32;
                case 'u': case 'U': case 't': case 'T':
                    return 33;
                case 'x': case 'X': case 'N': case 'n': case '?':
                    return 43;
                case 'r': case 'R':
                    return 44;
                case 'y': case 'Y':
                    return 45;
                case '-':
                    return 77;
                default:
                    return 82;
            }
        case 'x': case 'X': case 'n': case 'N': case '?':
            switch ( base[1] ) {
                case 'a': case 'A':
                    return 46;
                case 'g': case 'G':
                    return 52;
                case 'c': case 'C':
                    return 49;
                case 'u': case 'U': case 't': case 'T': 
                    return 55;
                case 'x': case 'X': case 'n': case 'N': case '?':
                    return 58;
                case 'r': case 'R':
                    return 59;
                case 'y': case 'Y':
                    return 60;
                case '-':
                    return 78;
                default:
                    return 82;
            }
        case 'r': case 'R':
            switch ( base[1] ) {
                case 'a': case 'A':
                    return 47;
                case 'g': case 'G':
                    return 53;
                case 'c': case 'C':
                    return 50;
                case 'u': case 'U': case 't': case 'T':
                    return 56;
                case 'x': case 'X': case 'n': case 'N': case '?':
                    return 64;
                case 'r': case 'R':
                    return 65;
                case 'y': case 'Y':
                    return 66;
                case '-':
                    return 79;
                default:
                    return 82;
            }
        case 'y': case 'Y':
            switch ( base[1] ) {
                case 'a': case 'A':
                    return 48;
                case 'c': case 'C':
                    return 51;
                case 'g': case 'G':
                    return 54;
                case 'u': case 'U': case 't': case 'T':
                    return 57;
                case 'x': case 'X': case 'n': case 'N': case '?':
                    return 61;
                case 'r': case 'R':
                    return 62;
                case 'y': case 'Y':
                    return 63;
                case '-':
                    return 80;
                default:
                    return 82;
            }
        case '-':
            switch ( base[1] ) {
                case 'a': case 'A':
                    return 67;
                case 'g': case 'G':
                    return 68;
                case 'c': case 'C':
                    return 69;
                case 'u': case 'U': case 't': case 'T':
                    return 70;
                case 'x': case 'X': case 'n': case 'N': case '?':
                    return 71;
                case 'r': case 'R':
                    return 72;
                case 'y': case 'Y':
                    return 73;
                case '-':
                    return 81;
                default:
                    return 82;
            }
        default :
            return 82;
    }
}

double TwoStateModel::getEmpiricalPropInvariant(SequenceTable* sequenceTable, int modelId ) {

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair<string, unsigned int> > & invariantBases = sequenceTable->getInvariantBases( modelId );
    int numberInvariant = 0;
    for ( vector< pair< string, unsigned int > >::const_iterator iter = invariantBases.begin();
          iter != invariantBases.end(); ++iter ){
        numberInvariant += (*iter).second;
    }
    return ( (double)numberInvariant/
             (double)(numberInvariant + sequences.numberColumns()) );
}

vector < double > TwoStateModel::retrieveEmpiricalFrequencies
( SequenceTable * sequenceTable, int modelId ) {

    vector < double > freq2( 2 );

    for ( int i = 0; i < 2; ++i ) {
        freq2[i] = 1e-6;
    }

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair< string, unsigned int > > & invariantBases = sequenceTable->getInvariantBases( modelId );
    const array2D< double > & equivalencyTable = getEquivalencyTable();
    unsigned int s;
    
    // use the equivalency table to match symbol/state    
    for ( unsigned int sequence = 0; sequence < sequenceTable->getNumberSpecies(); ++sequence ){
        for ( unsigned int site = 0; site < sequenceTable->getSequencesLength(modelId); ++site ) {
            s = getSymbolNumber( sequences( sequence, site ) );
            if (s == getNumberSymbols()){
                errorSymbol( sequenceTable, modelId, sequence, site );
            }
            double sum = 0.0;
            for ( int nucleotide = 0; nucleotide < 2; ++nucleotide ) {
                if ( equivalencyTable(s,nucleotide) == 1 ) {
                    sum += 1.0;
                }
            }
            for ( int nucleotide = 0; nucleotide < 2; ++nucleotide ) {
                if ( equivalencyTable(s,nucleotide) == 1 ) {
                    freq2[nucleotide] += 1.0/sum;
                }
            }
        }
    }
    for ( unsigned int site = 0; site < invariantBases.size(); ++site ){
        s = getSymbolNumber( invariantBases[site].first );
        if (s == getNumberSymbols()){
            errorSymbol( sequenceTable, modelId, 0, -site-1 );
        }
        double sum = 0.0;
        for ( int nucleotide = 0; nucleotide < 2; ++nucleotide ) {
            if ( equivalencyTable(s,nucleotide) == 1 ) {
                sum += 1.0;
            }
        }
        for ( int nucleotide = 0; nucleotide < 2; ++nucleotide ) {   
            if ( equivalencyTable(s,nucleotide) == 1 ) {
                freq2[nucleotide] += (invariantBases[site].second *
                                 sequenceTable->getNumberSpecies() / sum );
            }
        }
    }
    
    double sum = freq2[0] + freq2[1];
    for ( int nucleotide = 0; nucleotide < 2; ++nucleotide ) {
        freq2[nucleotide] /= sum;
    }
    return freq2;
}

string TwoStateModel::getSymbol( unsigned int basenumber, unsigned int ) const {
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
        case 16 :
            return ( string( "0" ) );
        case 17 :
            return ( string( "1" ) );
        // Watson-Crick base pairs
        case 18 :
            return ( string( "AU" ) );
        case 19 :
            return ( string( "GU" ) );
        case 20 :
            return ( string( "GC" ) );
        case 21 :
            return ( string( "UA" ) );
        case 22 :
            return ( string( "UG" ) );
        case 23 :
            return ( string( "CG" ) );
            // Mismatches
        case 24 :
            return ( string( "AA" ) );
        case 25 :
            return ( string( "AG" ) );
        case 26 :
            return ( string( "AC" ) );
        case 27 :
            return ( string( "GA" ) );
        case 28 :
            return ( string( "GG" ) );
        case 29 :
            return ( string( "CA" ) );
        case 30 :
            return ( string( "CC" ) );
        case 31 :
            return ( string( "CU" ) );
        case 32 :
            return ( string( "UC" ) );
        case 33 :
            return ( string( "UU" ) );
            // Ambiguous base-pairs
        case 34 :
            return ( string( "AX" ) );
        case 35 :
            return ( string( "AR" ) );
        case 36 :
            return ( string( "AY" ) );
        case 37 :
            return ( string( "CX" ) );
        case 38 :
            return ( string( "CR" ) );
        case 39 :
            return ( string( "CY" ) );
        case 40 :
            return ( string( "GX" ) );
        case 41 :
            return ( string( "GR" ) );
        case 42 :
            return ( string( "GY" ) );
        case 43 :
            return ( string( "UX" ) );
        case 44 :
            return ( string( "UR" ) );
        case 45 :
            return ( string( "UY" ) );
        case 46 :
            return ( string( "XA" ) );
        case 47 :
            return ( string( "RA" ) );
        case 48 :
            return ( string( "YA" ) );
        case 49 :
            return ( string( "XC" ) );
        case 50 :
            return ( string( "RC" ) );
        case 51 :
            return ( string( "YC" ) );
        case 52 :
            return ( string( "XG" ) );
        case 53 :
            return ( string( "RG" ) );
        case 54 :
            return ( string( "YG" ) );
        case 55 :
            return ( string( "XU" ) );
        case 56 :
            return ( string( "RU" ) );
        case 57 :
            return ( string( "YU" ) );
        case 58 :
            return ( string( "XX" ) );
        case 59 :
            return ( string( "XR" ) );
        case 60 :
            return ( string( "XY" ) );
        case 61 :
            return ( string( "YX" ) );
        case 62 :
            return ( string( "YR" ) );
        case 63 :
            return ( string( "YY" ) );
        case 64 :
            return ( string( "RX" ) );
        case 65 :
            return ( string( "RR" ) );
        case 66 :
            return ( string( "RY" ) );
            // "gap" base-pair
        case 67 :
            return ( string( "-A" ) );
        case 68 :
            return ( string( "-G" ) );
        case 69 :
            return ( string( "-C" ) );
        case 70 :
            return ( string( "-U" ) );
        case 71 :
            return ( string( "-X" ) );
        case 72 :
            return ( string( "-R" ) );
        case 73 :
            return ( string( "-Y" ) );
        case 74 :
            return ( string( "A-" ) );
        case 75 :
            return ( string( "G-" ) );
        case 76 :
            return ( string( "C-" ) );
        case 77 :
            return ( string( "U-" ) );
        case 78 :
            return ( string( "X-" ) );
        case 79 :
            return ( string( "R-" ) );
        case 80 :
            return ( string( "Y-" ) );
        case 81 :
            return ( string( "--" ) );
        default :
            return ( string( "Z" ) ); // unknown base-pair
    }
}

double TwoStateModel::getExchangeability( unsigned int,
                                         unsigned int, unsigned int ) const {
    return 1.0;
}
    

void TwoStateModel::updateEigenMatrix(){
    if ( discreteGamma ) {
        for ( unsigned int gammaCategory = 0 ;
                  gammaCategory < numberGammaCategories; ++gammaCategory ) {
            int freqCat = gammaCategory + invariantCategory;
            substitutionRate[gammaCategory] = averageRate[gammaCategory]/
                ( 2 * getFrequency( 0, freqCat ) * getFrequency( 1, freqCat ) );
        }
    }
    else {
        int freqCat = invariantCategory;
        substitutionRate[0] = averageRate[0] /
                ( 2 * getFrequency( 0, freqCat ) * getFrequency( 1, freqCat ) );
    }
    setEigenMatrix();
}

double TwoStateModel::probability( unsigned int oldState, unsigned int newState, double time,
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
    for ( k = 0; k < 2; k++ ) {
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


double TwoStateModel::diffProbability( unsigned int oldState, unsigned int newState, double time,
unsigned int category, unsigned int ) const {
    int k, cat;
    double diffProbability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) )
        return 0.0;
    cat = category - invariantCategory;
    for ( k = 0; k < 2; k++ ) {
        diffProbability += ( eigenValues[cat] [k] *
        ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    return ( diffProbability );
}


double TwoStateModel::secondDiffProbability( unsigned int oldState, unsigned int newState, double time,
unsigned int category, unsigned int ) const {
    int k, cat;
    double diff2Probability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) )
        return 0.0;
    cat = category - invariantCategory;
    for ( k = 0; k < 2; k++ ) {
        diff2Probability +=
        ( pow( eigenValues[cat] [k], 2.0 ) *
        ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    return ( diff2Probability );
}


void TwoStateModel::initEquivalencyTable(){
    equivalencyTable.resize( 82, 2 );
    // Set the patterns for valid states
    equivalencyTable( 0, 0 ) = 1.0; // A
    equivalencyTable( 0, 1 ) = 0.0;
    equivalencyTable( 1, 0 ) = 0.0; // C
    equivalencyTable( 1, 1 ) = 1.0;
    equivalencyTable( 2, 0 ) = 1.0; // G
    equivalencyTable( 2, 1 ) = 0.0;
    equivalencyTable( 3, 0 ) = 0.0; // T/U
    equivalencyTable( 3, 1 ) = 1.0;
    // Match a Purine
    equivalencyTable( 4, 0 ) = 1.0; // A/G
    equivalencyTable( 4, 1 ) = 0.0;
    // Match a Pyramidine
    equivalencyTable( 5, 0 ) = 0.0;
    equivalencyTable( 5, 1 ) = 1.0; // C/U
    // Match all other possible doublet/triplet/unknown/gap combination
    for (unsigned int i = 6; i < 16; ++i){
        equivalencyTable( i, 0 ) = 1.0;
        equivalencyTable( i, 1 ) = 1.0;
    }
    //binary values
    equivalencyTable( 16, 0 ) = 1.0; // 0
    equivalencyTable( 16, 1 ) = 0.0;
    equivalencyTable( 17, 0 ) = 0.0; // 1
    equivalencyTable( 17, 1 ) = 1.0;
    // AU, GU, GC
    equivalencyTable( 18, 0 ) = 1.0;
    equivalencyTable( 18, 1 ) = 0.0;
    equivalencyTable( 19, 0 ) = 1.0;
    equivalencyTable( 19, 1 ) = 0.0;
    equivalencyTable( 20, 0 ) = 1.0;
    equivalencyTable( 20, 1 ) = 0.0;
    // UA, UG, CG
    equivalencyTable( 21, 0 ) = 0.0;
    equivalencyTable( 21, 1 ) = 1.0;
    equivalencyTable( 22, 0 ) = 0.0;
    equivalencyTable( 22, 1 ) = 1.0;
    equivalencyTable( 23, 0 ) = 0.0;
    equivalencyTable( 23, 1 ) = 1.0;
    //all other base-pairs are either mismatches, ambiguity or with gaps
    //treat them as unknown
    for (unsigned int i = 24; i < 82; ++i){
        equivalencyTable( i, 0 ) = 1.0;
        equivalencyTable( i, 1 ) = 1.0;
    }
    //except two of them ( RU = AU or GU, and UR = UA or UG)
    equivalencyTable( 56, 0 ) = 1.0;
    equivalencyTable( 56, 1 ) = 0.0;
    equivalencyTable( 44, 0 ) = 0.0;
    equivalencyTable( 44, 1 ) = 1.0;    
}


// print the parameters in readable format
void TwoStateModel::printParameters( ostream & outputStream ) const {
  
    //print common parameters
    MatrixModel::printParameters( outputStream );

    outputStream.setf(ios::fixed);
    outputStream << setprecision(5);
    for ( unsigned int freqCat = 0; freqCat < frequencies->getNumberFrequenciesSets(); ++freqCat ) {
        if (frequencies->getNumberFrequenciesSets()>=2){
            outputStream << "set of frequencies " << freqCat + 1 << ":" << endl;
        }
        else{
            outputStream << "frequencies:" << endl;
        }
        for ( unsigned int i = 0; i < numberFrequencies; ++i ) {
            outputStream << "f[" << getFrequencyState( i ) << "] = " <<
                getFrequency(i, freqCat) << "  ";
        }
        outputStream << endl;
    }
    outputStream << endl;
    outputStream << "no rate parameter with two-state models" << endl;   
    outputStream << endl;
    outputStream.unsetf(ios::fixed);
}
