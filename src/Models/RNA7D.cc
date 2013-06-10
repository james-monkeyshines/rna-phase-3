#include "Models/RNA7D.h"
using namespace std;

RNA7D RNA7D::prototype( "RNA7D" );
RNA7D::RNA7D( const string & registrationName ): RNA7A( registrationName ) {}
RNA7D::RNA7D( ParametersSet & parameters ): RNA7A( parameters ) {}

bool RNA7D::setParameterNumbers() {
    numberRatesRatios = 3;
    numberFrequencies = 7;
	bool bp_symmetry = false;
	return bp_symmetry;
}

vector < double > RNA7D::condenseRates( array2D < double > & rates16x16 ) {
	vector < double > rates(3);
	
	double ref = ( rates16x16( 0, 2 ) + rates16x16( 3, 5 ) ) / 2;
	
	// 0 - single transition rate ratio
	rates[0] = ( rates16x16( 0, 1 ) + rates16x16( 1, 2 ) + rates16x16( 3, 4 ) + rates16x16( 4, 5 ) );
	rates[0] /= (4 * ref);
	
	// 1 - double transversion rate ratio
	for ( int i = 0; i < 3; ++i ) {
		for ( int j = 3; j < 6; ++j ) {
			rates[1] += rates16x16( i, j );
		}
	}
	rates[1] /= (9 * ref);
	
	double total_mismatch_rate = 0.0;
	for ( int i = 0; i < 6 ; ++i ) {
		for ( int j = 6; j < 16 ; ++j ) {
			total_mismatch_rate += rates16x16( i, j );
		}
	}
	rates[2] = total_mismatch_rate /= (60 * ref);
	
	return rates;
}

// Destructor
RNA7D::~RNA7D() {
}

string RNA7D::getName( void ) const {
    string name = "RNA7D" + MatrixModel::getName();
    return ( name );
}

void RNA7D::initMatrixIndex() {
    matrixIndex.resize( 7, 7 );  
    //the reference ratio is given the negative index of -2
    matrixIndex( 0, 2 ) = matrixIndex(3, 5) = -2;
    // 0 - single transition rate ratio
    matrixIndex(0, 1) = matrixIndex(1, 2) =
        matrixIndex(3, 4) = matrixIndex(4, 5) = 0;
    // 1 - double transversion rate ratio
    for ( int i = 0; i < 3; ++i ) {
        for ( int j = 3; j < 6; ++j ) {
            matrixIndex(i, j) = 1;
        }
    }
    // 2 - mismatch rate ratio
    for ( int i = 0; i < 6; ++i ) {
        matrixIndex(i, 6) = 2;
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

Model * RNA7D::clone( ParametersSet & parameters ) const {
    return new RNA7D( parameters );
}
