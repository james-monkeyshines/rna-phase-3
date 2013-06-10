#include "Models/JC69.h"

#include <assert.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the model to the model factory with the name JC69
JC69 JC69::prototype( "JC69" );

JC69::JC69( const string & registrationName ) : DnaModel( registrationName ){
//the private constructor is called by the prototype only
}

JC69::JC69( ParametersSet & parameters ) : DnaModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void JC69::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 0;
    numberFrequencies = 0;
        
    //call the basic initialisation method
    DnaModel::initialisation( sequenceTable, modelId );
    
    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();
    
    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if (sequenceTable){
        //estimation of the proportion of invariant sites (if relevant)
        //from the empirical sequences
        frequencies->initialisation( 0 );
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
  
JC69::~JC69() {
}

string JC69::getName( void ) const {
    return ( "JC69" + MatrixModel::getName() );
}

double JC69::getExchangeability( unsigned int, unsigned int, unsigned int )
                                                                    const {
    return 1.0;
}

void JC69::updateEigenMatrix(){
    if ( discreteGamma ) {
        for ( unsigned int i = 0 ; i < numberGammaCategories; ++i ) {
            substitutionRate[i] = averageRate[i] / .750;
        }
    }
    else {
        substitutionRate[0] = averageRate[0] / .750;
    }
    setEigenMatrix();
}


void JC69::printParameters( ostream & outputStream ) const{
    MatrixModel::printParameters( outputStream );
    outputStream << "f[A] = .25    f[C] = .25    f[G] = .25    f[U] = .25    "
                 << endl;
    outputStream << "no rate parameter with the Jukes-Cantor model" << endl;
    MatrixModel::printParametersAux( outputStream );
}



void JC69::setEigenMatrix() {
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
    category++ ) {
        rateMatrix[category] ( 0, 0 ) = substitutionRate[category] * -.75;
        rateMatrix[category] ( 0, 1 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 0, 2 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 0, 3 ) = substitutionRate[category] * .25;

        rateMatrix[category] ( 1, 0 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 1, 1 ) = substitutionRate[category] * -.75;
        rateMatrix[category] ( 1, 2 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 1, 3 ) = substitutionRate[category] * .25;

        rateMatrix[category] ( 2, 0 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 2, 1 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 2, 2 ) = substitutionRate[category] * -.75;
        rateMatrix[category] ( 2, 3 ) = substitutionRate[category] * .25;

        rateMatrix[category] ( 3, 0 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 3, 1 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 3, 2 ) = substitutionRate[category] * .25;
        rateMatrix[category] ( 3, 3 ) = substitutionRate[category] * -.75;
        
        for ( int i = 0; i < 4; ++i ) {
            eigenValues[category] [i] = 0.0;
            for ( int j = 0; j < 4; ++j ) {
                eigenMatrix[category] ( i, j ) = 0.0;
                ieigenMatrix[category] ( i, j ) = 0.0;
            }
        }
        eigenValues[category] [0] = 0.0;
        eigenValues[category] [1] = -substitutionRate[category];
        eigenValues[category] [2] = -substitutionRate[category];
        eigenValues[category] [3] = -substitutionRate[category];


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

Model * JC69::clone( ParametersSet & parameters ) const {
    if( parameters.findParameter( "Number of frequencies sets" ) ){
        if (parameters.intParameter("Number of frequencies sets") != 0){
            cerr << "The \"Number of frequencies sets\" parameter "
                 << "cannot be used with a JC69 model" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        parameters["Number of frequencies sets"] = "0";
    }
    if( parameters.findParameter( "Number of rates ratios sets" ) ){
        if ( parameters.intParameter("Number of rates ratios sets") != 0 ){
            cerr << "The \"Number of rates ratios sets\" parameter "
                 << "cannot be used with a JC69 model" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        parameters["Number of rates ratios sets"] = "0";
    }
    return new JC69( parameters );
}
