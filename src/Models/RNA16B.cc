#include "Models/RNA16B.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name RNA16B
RNA16B RNA16B::prototype( "RNA16B" );

RNA16B::RNA16B( const string & registrationName ) :
RNA16( registrationName ) {
//the private constructor is called by the prototype only
//RNA16B is a restriction of RNA16
}

RNA16B::RNA16B( ParametersSet & parameters )
: RNA16( parameters ) {
//the normal constructor is called by the prototype's clone method
//RNA16B is a restriction of RNA16
}


void RNA16B::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //The initialisation procedure for the 16B model is slightly different
    //than for the 16 model, therefore we redefine it
    //initialise the number of parameters
    numberRatesRatios = 0;
    numberFrequencies = 16;
    
    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();
    
    initMatrixIndex();

    //call the basic initialisation method
    RnaModel::initialisation( sequenceTable, modelId );
    
    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if ( sequenceTable ) {
        frequencies->initialisation( retrieveEmpiricalFrequencies( sequenceTable, modelId ) );
        //estimation of all the rate ratios (and proportion of invariant sites
        //if relevant from the empirical sequences
        ratesRatios->initialisation( 0 );
        if ( invariantCategory ) {
            retrieveEmpiricalRates( sequenceTable, modelId,
                                   & proportionInvariantSites );
        }
        else {
            //not useful
        }
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}


// Destructor
RNA16B::~RNA16B() {
}


// The models unique name
string RNA16B::getName( void ) const {
    return ( "RNA16B" + MatrixModel::getName() );
}


void RNA16B::printParameters( ostream & outputStream ) const{
    MatrixModel::printParameters( outputStream );
    outputStream << "no rate parameter with the 16B model" << endl;
    MatrixModel::printParametersAux( outputStream );
}


void RNA16B::initMatrixIndex(){
    matrixIndex.resize( 16, 16 );  
    //initialise the diagonal (-3)
    for ( int i = 0 ; i < 16 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
    //fill the matrix for double substitution (-1 for rate = 0)
    for ( int i = 0; i < 15 ; ++i ) {
        for ( int j = i+1; j < 16 ; ++j ) {
            matrixIndex(i, j) = -1;
        }
    }
    
    //the reference ratio is given the negative index of -2
    matrixIndex( 0, 1 ) = matrixIndex(2, 8) =
        matrixIndex(3, 4) = matrixIndex(5, 11) =
        matrixIndex(6, 7) = matrixIndex(6, 9) =
        matrixIndex(7, 10) = matrixIndex(9, 10) = -2;
    // second transition rate ratio C<->U
    matrixIndex( 0, 8 ) = matrixIndex( 1, 2 ) =
        matrixIndex( 3, 11 ) = matrixIndex( 4, 5 ) =
        matrixIndex( 12, 13 ) = matrixIndex( 12, 14 ) =
        matrixIndex( 13, 15 ) = matrixIndex( 14, 15 ) = -2;
    matrixIndex( 0, 13 ) = matrixIndex( 2, 9 ) =
        matrixIndex( 3, 14 ) = matrixIndex( 5, 7 ) =
        matrixIndex( 6, 8 ) = matrixIndex( 6, 11 ) =
        matrixIndex( 8, 12 ) = matrixIndex( 11, 12 ) = -2;
   matrixIndex( 0, 6 ) = matrixIndex( 0, 15 ) =
        matrixIndex( 1, 9 ) = matrixIndex( 3, 6 ) =
        matrixIndex( 3, 15 ) = matrixIndex( 4, 7 ) =
        matrixIndex( 8, 14 ) = matrixIndex( 11, 13 ) = -2;
   matrixIndex( 1, 13 ) = matrixIndex( 2, 10 ) =
        matrixIndex( 2, 12 ) = matrixIndex( 4, 14 ) =
        matrixIndex( 5, 10 ) = matrixIndex( 5, 12 ) =
        matrixIndex( 7, 8 ) = matrixIndex( 9, 11 ) = -2;
    matrixIndex( 0, 7 ) = matrixIndex( 1, 10 ) =
        matrixIndex( 1, 15 ) = matrixIndex( 2, 14 ) =
        matrixIndex( 3, 9 ) = matrixIndex( 4, 10 ) =
        matrixIndex( 4, 15 ) = matrixIndex( 5, 13 ) = -2;
    
    //symetrize the matrix
    for ( int i = 0; i < 15 ; ++i ) {
        for ( int j = i+1; j < 16 ; ++j ) {
            matrixIndex(j, i) = matrixIndex(i, j);
        }
    }
}

Model * RNA16B::clone( ParametersSet & parameters ) const {
    if( parameters.findParameter( "Number of rates ratios sets" ) ){
        if ( parameters.intParameter("Number of rates ratios sets") != 0 ){
            cerr << "The \"Number of rates ratios sets\" parameter "
                 << "cannot be used with a RNA16B model" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        parameters["Number of rates ratios sets"] = "0";
    }
    return new RNA16B( parameters );
}

