#include "Models/RNA7C.h"
using namespace std;

RNA7C RNA7C::prototype( "RNA7C" );
RNA7C::RNA7C( const string & registrationName ): RNA7A( registrationName ) {}
RNA7C::RNA7C( ParametersSet & parameters ): RNA7A(parameters) {}

bool RNA7C::setParameterNumbers() {
    numberRatesRatios = 9;
    numberFrequencies = 7;
	bool bp_symmetry = false;
	return bp_symmetry;
}

vector < double > RNA7C::condenseRates( array2D < double > & rates16x16 ) {
	vector < double > rates(9);
	double ref = rates16x16( 0, 1 );
	rates[0] = rates16x16( 1, 2 );
	rates[1] = rates16x16( 3, 4 );
	rates[2] = rates16x16( 4, 5 );
	
	for ( int i = 0; i < 6 ; ++i ) {
		double total_mismatch_rate = 0.0;
		for ( int j = 6; j < 16 ; ++j ) {
			total_mismatch_rate += rates16x16( i, j );
		}
		rates[i+3] = total_mismatch_rate /= 10;
	}
	
	for (unsigned int i = 0; i < 9; ++i){
		rates[i] /= ref;
	}
	
	return rates;
}

// Destructor
RNA7C::~RNA7C() {
}

string RNA7C::getName( void ) const {
    string name = "RNA7C" + MatrixModel::getName();
    return ( name );
}

void RNA7C::initMatrixIndex() {
    matrixIndex.resize( 7, 7 );  
    //the reference ratio is given the negative index of -2
    matrixIndex( 0, 1 ) = -2;
    //other single-transition
    matrixIndex( 1, 2 ) = 0;
    matrixIndex( 3, 4 ) = 1;
    matrixIndex( 4, 5 ) = 2;
    //single-transition to/from MM
    for (unsigned int i = 0; i < 6; ++i){
        matrixIndex( i, 6 ) = i+3;
    }
    //forbidden double-transitions
    matrixIndex(0, 2) = matrixIndex(3, 5) = -1;
    for ( int i = 0; i < 3; ++i ) {
        for ( int j = 3; j < 6; ++j ) {
            matrixIndex(i, j) = -1;
        }
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

Model * RNA7C::clone( ParametersSet & parameters ) const {
    return new RNA7C( parameters );
}
