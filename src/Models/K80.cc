#include "Models/K80.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name K80
K80 K80::prototype( "K80" );

K80::K80( const string & registrationName ) : DnaModel( registrationName ){
//the private constructor is called by the prototype only
}

K80::K80( ParametersSet & parameters ) : DnaModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void K80::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 1;
    numberFrequencies = 0;
    
    //initialise the matrixIndex for correspondance between a matrix element
    //and an index in the rate ratios vector
    initMatrixIndex();
    
    //call the basic initialisation method
    DnaModel::initialisation( sequenceTable, modelId );
    
    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();
    
    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if (sequenceTable){
        frequencies->initialisation( 0 );
        //estimation of all the rate ratios (and proportion of invariant sites
        //if relevant) from the empirical sequences
        vector < double > rates5;
        if ( invariantCategory ) {
            rates5 = retrieveEmpiricalRates( sequenceTable, modelId,
            & proportionInvariantSites );
        }
        else {
            rates5 = retrieveEmpiricalRates( sequenceTable, modelId );
        }
        //the transversion rate is the only free parameter
        vector < double > rate1;
        rate1.push_back( (rates5[0]+rates5[1]+rates5[2]+rates5[4])/(2*(1.0+rates5[3])) );
        ratesRatios->initialisation( rate1 );
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();    
}
  
K80::~K80() {
}

string K80::getName( void ) const {
    return ( "K80" + MatrixModel::getName() );
}

void K80::updateEigenMatrix(){
    double sum = 0.0;
    if ( discreteGamma ) {
        for ( unsigned int i = 0 ; i < numberGammaCategories; ++i ) {
            sum = .25 + .5 * (*ratesRatios)[i][0];
            substitutionRate[i] = averageRate[i]/sum;
        }
    }
    else {
        sum = .25 + .5 * (*ratesRatios)[0][0];
        substitutionRate[0] = averageRate[0]/sum;
    }
    setEigenMatrix();
}


void K80::printParameters( ostream & outputStream ) const{
    //print common parameters
    MatrixModel::printParameters( outputStream );
    outputStream << "f[A] = .25   f[C] = .25   f[G] = .25   f[U] = .25" << endl;
    MatrixModel::printParametersAux( outputStream );
}

void K80::initMatrixIndex() {
    matrixIndex.resize( 4, 4 );
    for ( unsigned int i = 0 ; i < 4 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
    //reference (transitions)
    matrixIndex(0, 2) = matrixIndex(2, 0) = -2;
    matrixIndex(1, 3) = matrixIndex(3, 1) = -2;
    //transversion
    matrixIndex(0, 1) = matrixIndex(1, 0) = 0;
    matrixIndex(0, 3) = matrixIndex(3, 0) = 0;
    matrixIndex(1, 2) = matrixIndex(2, 1) = 0;
    matrixIndex(2, 3) = matrixIndex(3, 2) = 0;
}


void K80::setEigenMatrix() {
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
    category++ ) {
        // fill the rateMatrix (cf Hky.h)
        int rateCat = ratesRatios->ratesRatiosCat(category);
        rateMatrix[category] ( 0, 1 ) = substitutionRate[category] * .25 *
                                        (*ratesRatios)(rateCat) [0];
        rateMatrix[category] ( 0, 2 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 0, 3 ) = rateMatrix[category] ( 0, 1 ) ;
        rateMatrix[category] ( 0, 0 ) = -2 * rateMatrix[category]( 0, 1 ) -
                                        rateMatrix[category]( 0, 2 );
        rateMatrix[category] ( 1, 0 ) = rateMatrix[category] ( 0, 1 );
        rateMatrix[category] ( 1, 2 ) = rateMatrix[category] ( 0, 1 );
        rateMatrix[category] ( 1, 3 ) = rateMatrix[category] ( 0, 2 );
        rateMatrix[category] ( 1, 1 ) = rateMatrix[category] ( 0, 0 );

        rateMatrix[category] ( 2, 0 ) = rateMatrix[category] ( 0, 2 );
        rateMatrix[category] ( 2, 1 ) = rateMatrix[category] ( 0, 1 );
        rateMatrix[category] ( 2, 3 ) = rateMatrix[category] ( 0, 1 );
        rateMatrix[category] ( 2, 2 ) = rateMatrix[category] ( 0, 0 );

        rateMatrix[category] ( 3, 0 ) = rateMatrix[category] ( 0, 1 );
        rateMatrix[category] ( 3, 1 ) = rateMatrix[category] ( 0, 2 );
        rateMatrix[category] ( 3, 2 ) = rateMatrix[category] ( 0, 1 );
        rateMatrix[category] ( 3, 3 ) = rateMatrix[category] ( 0, 0 );

        for ( unsigned int i = 0; i < 4; ++i ) {
            eigenValues[category] [i] = 0.0;
            for ( unsigned int j = 0; j < 4; ++j ) {
                eigenMatrix[category] ( i, j ) = 0.0;
                ieigenMatrix[category] ( i, j ) = 0.0;
            }
        }
        
        double transitionRate = substitutionRate[category];
        double transversionRate = substitutionRate[category] * (*ratesRatios)(rateCat)[0];
        

        eigenValues[category] [0] = 0.0;
        eigenValues[category] [1] =
        -transitionRate * .5 - transversionRate * .5;
        eigenValues[category] [2] =
        -transitionRate * .5 - transversionRate * .5;
        eigenValues[category] [3] = -transversionRate;


        // Update the eigenmatrix by hand
        eigenMatrix[category] ( 0, 0 ) = .25;
        eigenMatrix[category] ( 0, 1 ) = .25;
        eigenMatrix[category] ( 0, 2 ) = .25;
        eigenMatrix[category] ( 0, 3 ) = .25;

        eigenMatrix[category] ( 1, 0 ) = 0.0;
        eigenMatrix[category] ( 1, 1 ) = -1.0;
        eigenMatrix[category] ( 1, 2 ) = 0.0;
        eigenMatrix[category] ( 1, 3 ) = 1.0;

        eigenMatrix[category] ( 2, 0 ) = -1.0;
        eigenMatrix[category] ( 2, 1 ) = 0.0;
        eigenMatrix[category] ( 2, 2 ) = 1.0;
        eigenMatrix[category] ( 2, 3 ) = 0.0;

        eigenMatrix[category] ( 3, 0 ) = -.25;
        eigenMatrix[category] ( 3, 1 ) = .25;
        eigenMatrix[category] ( 3, 2 ) = -.25;
        eigenMatrix[category] ( 3, 3 ) = .25;

        // Inverse of the eigenmatrix
        ieigenMatrix[category] ( 0, 0 ) = 1.0;
        ieigenMatrix[category] ( 0, 1 ) = 0.0;
        ieigenMatrix[category] ( 0, 2 ) = -.5;
        ieigenMatrix[category] ( 0, 3 ) = -1.0;

        ieigenMatrix[category] ( 1, 0 ) = 1.0;
        ieigenMatrix[category] ( 1, 1 ) = -.5;
        ieigenMatrix[category] ( 1, 2 ) = 0.0;
        ieigenMatrix[category] ( 1, 3 ) = 1.0;

        ieigenMatrix[category] ( 2, 0 ) = 1.0;
        ieigenMatrix[category] ( 2, 1 ) = 0.0;
        ieigenMatrix[category] ( 2, 2 ) = .5;
        ieigenMatrix[category] ( 2, 3 ) = -1.0;

        ieigenMatrix[category] ( 3, 0 ) = 1.0;
        ieigenMatrix[category] ( 3, 1 ) = .5;
        ieigenMatrix[category] ( 3, 2 ) = 0.0;
        ieigenMatrix[category] ( 3, 3 ) = 1.0;
    } //for each category (except invariant)
}


Model * K80::clone( ParametersSet & parameters ) const {
    if( parameters.findParameter( "Number of frequencies sets" ) ){
        if (parameters.intParameter("Number of frequencies sets") != 0){
            cerr << "The \"Number of frequencies sets\" parameter "
                 << "cannot be used with a K80 model" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        parameters["Number of frequencies sets"] = "0";
    }
    return new K80( parameters );
}

