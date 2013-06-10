#include "Models/CodonModel.h"

#include <iostream>
#include <iomanip>
#include <algorithm>

#include "Sequence/SequenceTable.h"

#include "PatternDesign/Factory.h"

#include "PatternDesign/Singleton.h"
#include "Util/statlib.h"
#include "Util/matrixmath.h"
#include "Util/array2D.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"


using namespace std;

CodonModel::CodonModel(){
}

CodonModel::CodonModel( const string & registrationName ) :
        CodonParent( registrationName ){
}

CodonModel::CodonModel( ParametersSet & parameters ) : CodonParent( parameters ) {
    geneticCodeFileName = parameters.stringParameter("Genetic code file");
    readGeneticCode(geneticCodeFileName);
    matrixSize = geneticCode.size();
    initMatrix( matrixSize );
}

CodonModel::~CodonModel() {
}

string CodonModel::getName() const{
    return '(' + geneticCodeFileName + ')' + MatrixModel::getName();
}

string CodonModel::getState( unsigned int stateNumber, unsigned int ) const{
    assert(stateNumber<matrixSize);
    return geneticCode[stateNumber].first;
}

array2D< double > CodonModel::retrieveEmpiricalRates( SequenceTable * sequenceTable,
int modelId, double * propInvariant ) {

    array2D < double > ratesRat( matrixSize, matrixSize );

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair< string, unsigned int > > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s1, s2;

    for ( unsigned int i = 0 ; i < matrixSize; ++i ) {
        for ( unsigned int j = 0 ; j < matrixSize; ++j ) {
            ratesRat( i, j ) = 1e-6;
        }
    }

    for ( int site = -((int)invariantBases.size());
            site < (int)(sequences.numberColumns()); ++site ) {
        for ( unsigned int seq = 0;
              seq < sequenceTable->getNumberSpecies(); ++seq ) {
            if ( site >= 0 ){
                s1 = CodonParent::getSymbolNumber( sequences( seq, site ) );
            }
            else{
                s1 = CodonParent::getSymbolNumber(
                        invariantBases[-1-site].first );
            }
            if (s1 == CodonParent::getNumberSymbols()){
                errorSymbol( sequenceTable, modelId, seq, site );
            }
            for ( unsigned int seq2 = seq + 1;
                  seq2 < sequenceTable->getNumberSpecies(); ++seq2 ) {
                if ( site >= 0 ){
                    s2 = CodonParent::getSymbolNumber( sequences( seq2, site ) );
                }
                else{
                    s2 = s1;
                }
                if (s2 == CodonParent::getNumberSymbols()){
                    errorSymbol( sequenceTable, modelId, seq2, site );
                }
                double sum = 0.0;
                for ( unsigned int state1 = 0; state1 < matrixSize; ++state1 ){
                    for ( unsigned int state2 = 0; state2 < matrixSize; ++state2 ){
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
                for ( unsigned int state1 = 0; state1 < matrixSize; ++state1 ){
                    for ( unsigned int state2 = 0; state2 < matrixSize; ++state2 ){
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
    for ( unsigned int i = 0 ; i < matrixSize; ++i ) {
        for ( unsigned int j = i + 1; j < matrixSize; ++j ) {
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


vector < double > CodonModel::retrieveEmpiricalFrequencies( SequenceTable *
sequenceTable, int modelId ) {

    vector < double > freqVec( matrixSize );

    // Find empirical frequencies
    for ( unsigned int i = 0; i < matrixSize; ++i ) {
        //pseudo-count
        //freqVec[i] = 1.0/sequenceTable->getNumberSpecies();
        freqVec[i] = 1.0/matrixSize;
    }

    const array2D< string > & sequences = sequenceTable->getSequences( modelId );
    const vector< pair< string, unsigned int > > & invariantBases = sequenceTable->getInvariantBases( modelId );
    unsigned int s;
    double sum = 0.0;

    for ( unsigned int sequence = 0; sequence < sequenceTable->getNumberSpecies(); ++sequence ) {
        for ( unsigned int site = 0; site < sequenceTable->getSequencesLength(modelId); ++site ) {
            s = CodonParent::getSymbolNumber( sequences( sequence, site ) );
            if (s == CodonParent::getNumberSymbols()){
               errorSymbol( sequenceTable, modelId, sequence, site );
            }

			sum = 0.0;
            for ( unsigned int codon = 0; codon < matrixSize; ++codon ) {
                if ( equivalencyTable(s,codon) == 1.0 ) {
                    sum += 1.0;
                }
            }
            for ( unsigned int codon = 0; codon < matrixSize; ++codon ) {
                if ( equivalencyTable(s,codon) == 1.0 ) {
                    freqVec[codon] += 1.0/sum;
                }
            }
        }
    }
    for ( unsigned int site = 0; site < invariantBases.size(); ++site ){
        s = CodonParent::getSymbolNumber( invariantBases[site].first );
        if (s == CodonParent::getNumberSymbols()){
            errorSymbol( sequenceTable, modelId, 0, -site-1 );
        }

		sum = 0.0;
        for ( unsigned int codon = 0; codon < matrixSize; ++codon ) {
            if ( equivalencyTable(s,codon) == 1.0 ) {
                sum += 1.0;
            }
        }
        for ( unsigned int codon = 0; codon < matrixSize; ++codon ) {
            if ( equivalencyTable(s,codon) == 1.0 ) {
                freqVec[codon] += (invariantBases[site].second *
                                 sequenceTable->getNumberSpecies() / sum );
            }
        }
    }

    sum = 0.0;
    for ( unsigned int codon = 0; codon < matrixSize; ++codon ) {
        sum +=  freqVec[codon];
    }
    for ( unsigned int codon = 0; codon < matrixSize; ++codon ) {
        freqVec[codon] /= sum;
    }

    return freqVec;
}

void CodonModel::initEquivalencyTable(){

    assert(matrixSize);
    
    equivalencyTable.resize(getNumberSymbols(),matrixSize);

    // initialization
    for (unsigned int i = 0; i < getNumberSymbols(); i++ ){
        for (unsigned int j = 0; j < matrixSize; j++ ){
            equivalencyTable( i, j ) = 0.0;
        }
    }
    
    vector<unsigned int> codonMask;
    codonMask.reserve(geneticCode.size());
    //for all identified codons
    for (vector< pair<string,char> >::const_iterator iter = geneticCode.begin(); iter != geneticCode.end(); ++iter ){
        unsigned int mask = (CodonParent::getSingleSymbolNumber(iter->first[0])<<8) |
                            (CodonParent::getSingleSymbolNumber(iter->first[1])<<4) |
                            CodonParent::getSingleSymbolNumber(iter->first[2]);
        codonMask.push_back(mask);
    }
    
    
    for (unsigned int k = 0; k < CodonParent::getNumberSymbols(); ++k){
        unsigned int symbolMask = CodonParent::getMask(k);
        for ( unsigned int cod = 0; cod < codonMask.size(); ++cod ){
            //codonMask should be a set of "pure" symbols with 3 well definite nucleotides
            if ( (symbolMask & codonMask[cod]) == codonMask[cod] ) {
                equivalencyTable( k, cod ) = 1.0;
            }
        }
    }
}

double CodonModel::probability( unsigned int oldState, unsigned int newState,
                 double time, unsigned int category, unsigned int ) const {
    if ( ( invariantCategory == 1 ) && ( category == 0 ) ) {
        if ( oldState == newState ) {
            return 1.0;
        }
        else {
            return 0.0;
        }
    }

    int cat = category - invariantCategory;
    double probability = 0.0;
    for ( unsigned int k = 0; k < matrixSize; ++k )
        probability += ( ieigenMatrix[cat] ( oldState, k ) *
        eigenMatrix[cat] ( k, newState ) * exp( eigenValues[cat] [k] * time ) );

    // Just make sure
    //assert ( ( probability >= -0.000001 ) && ( probability <= 1.000001 ) );
    probability = MIN( probability, 1.0 );
    probability = MAX( probability, 0.0 );
    return ( probability );
}


double CodonModel::diffProbability( unsigned int oldState, unsigned int newState,
                   double time, unsigned int category, unsigned int ) const {
    if ( ( invariantCategory == 1 ) && ( category == 0 ) ) {
        return 0.0;
    }

    int cat = category - invariantCategory;
    double diffProbability = 0.0;

    for ( unsigned int k = 0; k < matrixSize; ++k ) {
        diffProbability += ( eigenValues[cat] [k] *
        ieigenMatrix[cat] ( oldState, k ) * eigenMatrix[cat] ( k, newState ) *
        exp( eigenValues[cat] [k] * time ) );
    }
    assert( !isnan( diffProbability ) );
    return ( diffProbability );
}


double CodonModel::secondDiffProbability( unsigned int oldState, unsigned int newState,
     double time, unsigned int category, unsigned int )  const {
    if ( ( invariantCategory == 1 ) && ( category == 0 ) ){
        return 0.0;
    }
    int cat = category - invariantCategory;
    double diff2Probability = 0.0;
    for ( unsigned int k = 0; k < matrixSize; ++k ) {
        diff2Probability += ( pow( eigenValues[cat] [k], 2.0 ) *
            ieigenMatrix[cat] ( oldState, k ) * eigenMatrix[cat] ( k, newState ) *
            exp( eigenValues[cat] [k] * time ) );
    }
    assert( !isnan( diff2Probability ) );
    return ( diff2Probability );
}


void CodonModel::printParameters( ostream & outputStream ) const {

    //print common parameters
    MatrixModel::printParameters( outputStream );

    //rely on the general method to print matrices
    MatrixModel::printParametersAux( outputStream );

}

void CodonModel::setEigenMatrix() {
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
