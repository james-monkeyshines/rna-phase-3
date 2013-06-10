#include "Models/RNA16_alt.h"
using namespace std;

// This model is identical to RNA16 in everything but the reference ratio.
RNA16_alt RNA16_alt::prototype( "RNA16_alt" );
RNA16_alt::RNA16_alt( const string & registrationName ): RNA16( registrationName ) {}
RNA16_alt::RNA16_alt( ParametersSet & parameters ): RNA16( parameters ) {}

// Destructor
RNA16_alt::~RNA16_alt() {
}

// The models unique name
string RNA16_alt::getName( void ) const {
    return ( "RNA16_alt" + MatrixModel::getName() );
}

void RNA16_alt::initMatrixIndex() {
    matrixIndex.resize( 16, 16 );
    int index = -1;
    for ( int i = 0 ; i < 15 ; ++i ) {
        for ( int j = i + 1 ; j < 16 ; ++j ) {
            if ( ( i == 0 ) && ( j == 1 ) ) {
                //the reference ratio is given the negative index of -2
                matrixIndex(i, j) = matrixIndex(j, i) = -2;
            }
            else {
                ++index;
                matrixIndex(i, j) = matrixIndex(j, i) = index;
            }
        }
    }
    assert( index == 118 );
    //initialise the diagonal
    for ( int i = 0 ; i < 16 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
}

