#include "Models/AaModel.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <algorithm>

#include "Sequence/SequenceTable.h"

#include "Models/Perturbator.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/statlib.h"
#include "Util/array2D.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"


using namespace std;

AaModel::AaModel():CodonParent(){
}

AaModel::AaModel( const string & registrationName ) :
        CodonParent( registrationName ){
}

AaModel::AaModel( ParametersSet & parameters ) : CodonParent( parameters ){
    if (parameters.findParameter("Genetic code file")){
        string geneticCodeFileName = parameters.stringParameter("Genetic code file");
        readGeneticCode(geneticCodeFileName);
    }
    initMatrix( 20 );
}

AaModel::~AaModel() {
}

string AaModel::getState( unsigned int stateNumber, unsigned int ) const {
    assert(stateNumber<20);
    //create a string with 1 char
    return string(&aminoAcid[stateNumber],1);
}

int AaModel::getSymbolNumber( const string & base, unsigned int ) const {

    if (base.length() == 1){
        vector<char>::const_iterator iter = find(aminoAcid.begin(),aminoAcid.end(),base[0]);
        if (iter!=aminoAcid.end()){
            return distance(aminoAcid.begin(),iter);
        }
        else{
            switch ( base[0] ) {
                case 'x': case 'X': case '?':
                    return 20; // Unknown amino acid
                case '-' :
                    return 21; // gaps are treated as unknown amino acids
                default:
                    return getNumberSymbols();
            }
        }
    }
        
    if (base.length() == 3){
        if (geneticCode.empty()){
            cerr << "Genetic code not loaded..." << endl;
            exit(EXIT_FAILURE);
        }
        unsigned int symbolNumber = CodonParent::getSymbolNumber( base );
        return (symbolNumber < CodonParent::getNumberSymbols() ) ? symbolNumber+22 : this->getNumberSymbols();
    }
    
    // Out of range (unrecognized symbol)
    return this->getNumberSymbols();
}

array2D < double > AaModel::retrieveEmpiricalRates( SequenceTable * sequenceTable,
int modelId, double * propInvariant ) {

    array2D < double > ratesRat( 20, 20 );

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair< string, unsigned int > > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s1, s2;

    for ( unsigned int i = 0 ; i < 20; ++i ) {
        for ( unsigned int j = 0 ; j < 20; ++j ) {
            ratesRat( i, j ) = 1e-6;
        }
    }

    for ( int site = -((int)invariantBases.size());
            site < (int)(sequences.numberColumns()); ++site ) {
        for ( unsigned int seq = 0;
              seq < sequenceTable->getNumberSpecies(); ++seq ) {
            if ( site >= 0 ){
                s1 = getSymbolNumber( sequences( seq, site ) );
            }
            else{
                s1 = getSymbolNumber(
                        invariantBases[-1-site].first );
            }
            if (s1 == getNumberSymbols()){
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
                for ( unsigned int state1 = 0; state1 < 20; ++state1 ){
                    for ( unsigned int state2 = 0; state2 < 20; ++state2 ){
                        if (state2!=state1){
                            if ( ( equivalencyTable(s1,state1) == 1.0 ) &&
                                 ( equivalencyTable(s2,state2) == 1.0 ) ){
                                sum += 1.0;
                            }
                        }
                    }
                }
                if ( site < 0 ){
                    sum /= (double)invariantBases[-1-site].second;
                }
                for ( unsigned int state1 = 0; state1 < 20; ++state1 ){
                    for ( unsigned int state2 = 0; state2 < 20; ++state2 ){
                        if (state2!=state1){
                            if ( ( equivalencyTable(s1,state1) == 1.0 ) &&
                                 ( equivalencyTable(s2,state2) == 1.0 ) ){
                                ratesRat( state1, state2 ) += 1.0/sum;
                            }
                        }
                    }
                }
            }
        }
    }

    double reference = ( ratesRat( 0, 2 ) + ratesRat( 2, 0 ) );
    for ( unsigned int i = 0 ; i < 20; ++i ) {
        for ( unsigned int j = i + 1; j < 20; ++j ) {
            ratesRat( i, j ) = ( ratesRat( i, j ) + ratesRat( j, i ) ) /
            ( reference );
            ratesRat( j, i ) = ratesRat( i, j );
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
    return ratesRat;
}

vector < double > AaModel::retrieveEmpiricalFrequencies
( SequenceTable * sequenceTable, int modelId ) {

    vector < double > freq20( 20 );

    for ( int i = 0; i < 20; ++i ) {
        //pseudo-count
        //freq20[i] = 1.0/sequenceTable->getNumberSpecies();
        freq20[i] = 0.05;
    }

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair< string, unsigned int > > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s;

    // count the number of each amino acid in the SequenceTable, gaps are ignored.
    for ( unsigned int sequence = 0; sequence < sequenceTable->getNumberSpecies(); ++sequence ){
        for ( unsigned int site = 0; site < sequenceTable->getSequencesLength(modelId); ++site ) {
            s = getSymbolNumber( sequences( sequence, site ) );
            if (s == getNumberSymbols()){
                errorSymbol( sequenceTable, modelId, sequence, site );
            }
            double sum = 0.0;
            for ( int aminoAcid = 0; aminoAcid < 20; ++aminoAcid ) {
                if ( equivalencyTable(s,aminoAcid) == 1 ) {
                    sum += 1.0;
                }
            }
            for ( int aminoAcid = 0; aminoAcid < 20; ++aminoAcid ) {
                if ( equivalencyTable(s,aminoAcid) == 1 ) {
                    freq20[aminoAcid] += 1.0/sum;
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
        for ( int aminoAcid = 0; aminoAcid < 20; ++aminoAcid ) {
            if ( equivalencyTable(s,aminoAcid) == 1 ) {
                sum += 1.0;
            }
        }
        for ( int aminoAcid = 0; aminoAcid < 20; ++aminoAcid ) {
            if ( equivalencyTable(s,aminoAcid) == 1 ) {
                freq20[aminoAcid] += (invariantBases[site].second *
                                 sequenceTable->getNumberSpecies() / sum );
            }
        }
    }

    double sum = 0.0;
    for (int i=0; i<20; ++i) sum += freq20[i] ;

    for ( int aminoAcid = 0; aminoAcid < 20; ++aminoAcid ) {
        freq20[aminoAcid] /= sum;
    }
    return freq20;
}

string AaModel::getSymbol( unsigned int baseNumber, unsigned int ) const {
    //unrecognized symbol first
    if ( baseNumber > getNumberSymbols() ) return string( "Z" );
    if ( baseNumber < 20 ) {
        return &(aminoAcid[baseNumber]);
    }
    switch(baseNumber){
        case 20 :
            return ( string( "X" ) ); // Unknown
        case 21 :
            return ( string( "-" ) ); // Gaps
        default : //must be a codon
            return CodonParent::getSymbol(baseNumber-21);
    }
}


double AaModel::probability( unsigned int oldState, unsigned int newState, double time,
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
    for ( k = 0; k < 20; k++ ) {
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


double AaModel::diffProbability( unsigned int oldState, unsigned int newState, double time,
unsigned int category, unsigned int ) const {
    int k, cat;
    double diffProbability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) )
        return 0.0;
    cat = category - invariantCategory;
    for ( k = 0; k < 20; k++ ) {
        diffProbability += ( eigenValues[cat] [k] *
        ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    return ( diffProbability );
}


double AaModel::secondDiffProbability( unsigned int oldState, unsigned int newState, double time,
unsigned int category, unsigned int ) const {
    int k, cat;
    double diff2Probability = 0.0;
    if ( ( invariantCategory == 1 ) && ( category == 0 ) )
        return 0.0;
    cat = category - invariantCategory;
    for ( k = 0; k < 20; k++ ) {
        diff2Probability +=
        ( pow( eigenValues[cat] [k], 2.0 ) *
        ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    return ( diff2Probability );
}


void AaModel::initEquivalencyTable(){
    unsigned int i, j;
    
    equivalencyTable.resize( this->getNumberSymbols(), this->getNumberStates() );
 
    // initialization
    for (i = 0; i < this->getNumberSymbols(); ++i ){
        for (j = 0; j < this->getNumberStates(); ++j ){
            equivalencyTable( i, j ) = 0.0;
        }
    }
 
    // Set the patterns for valid amino acids states
    for ( i=0; i < this->getNumberStates(); ++i ) {
        equivalencyTable( i, i ) = 1.0;
    }
 
    // Match all amino acids (unknown) and gaps
    for( j = 0; j < getNumberStates(); ++j){
        equivalencyTable( 20, j ) = 1.0;
        equivalencyTable( 21, j ) = 1.0;
    }
 
    // Set the patterns for codons (skipped if no genetic code)
    if (!geneticCode.empty()){
        vector<unsigned int> codonMask;
        vector<unsigned int> aa;
        codonMask.reserve(geneticCode.size());
        aa.reserve(geneticCode.size());
        //for all identified codons
        for (vector< pair<string,char> >::const_iterator iter = geneticCode.begin(); iter != geneticCode.end(); ++iter ){
            unsigned int mask = (CodonParent::getSingleSymbolNumber(iter->first[0])<<8) |
                                (CodonParent::getSingleSymbolNumber(iter->first[1])<<4) |
                                CodonParent::getSingleSymbolNumber(iter->first[2]);
            codonMask.push_back(mask);
            //find corresponding aa state
            vector<char>::iterator i = find(aminoAcid.begin(), aminoAcid.end(), iter->second);
            assert(i!=aminoAcid.end());
            aa.push_back(distance(aminoAcid.begin(), i));
        }
    
        for (unsigned int k = 22; k < 22+CodonParent::getNumberSymbols(); ++k){
            unsigned int symbolMask = CodonParent::getMask(k-22);
            for ( unsigned int cod = 0; cod < codonMask.size(); ++cod ){
                if ( (symbolMask & codonMask[cod]) == codonMask[cod] ) {
                    equivalencyTable( k, aa[cod] ) = 1.0;
                }
            }
        }
    }
}

void AaModel::printParameters( ostream & outputStream ) const {
    //print common parameters
    MatrixModel::printParameters( outputStream );
    //and subsitution rate matrices
    MatrixModel::printParametersAux( outputStream );
}
