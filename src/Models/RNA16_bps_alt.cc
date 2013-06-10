#include "Models/RNA16_bps_alt.h"
using namespace std;

RNA16_bps_alt RNA16_bps_alt::prototype( "RNA16_bps_alt" );
RNA16_bps_alt::RNA16_bps_alt( const string & registrationName ): RNA16_bps( registrationName ) {}
RNA16_bps_alt::RNA16_bps_alt( ParametersSet & parameters ): RNA16_bps( parameters ) {}

// Destructor
RNA16_bps_alt::~RNA16_bps_alt() {
}

// The models unique name
string RNA16_bps_alt::getName( void ) const {
    return ( "RNA16_bps_alt" + MatrixModel::getName() );
}

void RNA16_bps_alt::initMatrixIndex() {
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

