#include "Models/RNA7E.h"
using namespace std;

RNA7E RNA7E::prototype( "RNA7E" );
RNA7E::RNA7E( const string & registrationName ): RNA7D( registrationName ) {}
RNA7E::RNA7E( ParametersSet & parameters ): RNA7D(parameters) {}

bool RNA7E::setParameterNumbers() {
    numberRatesRatios = 1;
    numberFrequencies = 7;
	bool bp_symmetry = false;
	return bp_symmetry;
}

vector < double > RNA7E::condenseRates( array2D < double > & rates16x16 ) {
	vector < double > rates(1);
	
	double ref = ( rates16x16( 0, 1 ) + rates16x16( 1, 2 ) + rates16x16( 3, 4 ) + rates16x16( 4, 5 ) ) / 4.0;
	
	double total_mismatch_rate = 0.0;
	for ( int i = 0; i < 6; ++i ) {
		for ( int j = 6; j < 16 ; ++j ) {
			total_mismatch_rate += rates16x16( i, j );
		}
	}
	rates[0] = total_mismatch_rate / (60 * ref);
	
	return rates;
}

// Destructor
RNA7E::~RNA7E() {
}

string RNA7E::getName( void ) const {
    string name = "RNA7E" + MatrixModel::getName();
    return ( name );
}

void RNA7E::initMatrixIndex() {
    matrixIndex.resize( 7, 7 );  
    //the reference ratio is given the negative index of -2
    matrixIndex( 0, 1 ) = matrixIndex(3, 4) =
            matrixIndex( 1, 2 ) = matrixIndex( 4, 5 ) = -2;
    //forbidden double-transitions
    matrixIndex(0, 2) = matrixIndex(3, 5) = -1;
    for ( int i = 0; i < 3; ++i ) {
        for ( int j = 3; j < 6; ++j ) {
            matrixIndex(i, j) = -1;
        }
    }    
    //single-transition to/from MM
    for (unsigned int i = 0; i < 6; ++i){
        matrixIndex( i, 6 ) = 0;
    }
    //copy on the other triangle of the matrix
    for ( int i = 1; i < 7; ++i ) {
        for ( int j = 0; j < i; ++j ) {
            matrixIndex(i, j) = matrixIndex(j, i);
        }
    }
    //initialise the diagonal
    for ( int i = 0 ; i < 7 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
}

Model * RNA7E::clone( ParametersSet & parameters ) const {
    return new RNA7E( parameters );
}
