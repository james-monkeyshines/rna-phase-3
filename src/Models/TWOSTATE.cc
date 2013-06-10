#include <assert.h>

#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Models/Perturbator.h"
#include "Models/TWOSTATE.h"

#include "Sequence/SequenceTable.h"

#include "Util/statlib.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"

#define SQR(n) (n*n)

using namespace std;

//register the model to the model factory with the name REV
TWOSTATE TWOSTATE::prototype( "TWOSTATE" );

TWOSTATE::TWOSTATE( const string & registrationName ) : TwoStateModel( registrationName ){
//the private constructor is called by the prototype only
}


TWOSTATE::TWOSTATE( ParametersSet & parameters ) : TwoStateModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void TWOSTATE::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 0;
    numberFrequencies = 2;
            
    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();
    
    //call the basic initialisation method
    TwoStateModel::initialisation( sequenceTable, modelId );
        
    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if (sequenceTable){
        //estimation of the model frequencies from the empirical frequencies
        frequencies->initialisation(retrieveEmpiricalFrequencies( sequenceTable, modelId ));
        ratesRatios->initialisation( 0 );
        //estimation of the proportion of invariant sites (if relevant)
        //from the empirical sequences
        if ( invariantCategory ) {
            proportionInvariantSites = getEmpiricalPropInvariant( sequenceTable, modelId );
        }
    }
    else {
        //not useful
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}  

TWOSTATE::~TWOSTATE() {
}

string TWOSTATE::getName( void ) const {
    return ( "TWOSTATE" + MatrixModel::getName() );
}


void TWOSTATE::setEigenMatrix() {
    vector< double > junk( 2 );
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
          ++category ) {
        // fill the rateMatrix (cf REV.h)
        int freqCat = frequencies->frequenciesCat(category+invariantCategory);
        rateMatrix[category] ( 0, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1]; // Pur->Pyr
        rateMatrix[category] ( 0, 0 ) = -rateMatrix[category] ( 0, 1 );

        rateMatrix[category] ( 1, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0]; // Pyr->Pur
        rateMatrix[category] ( 1, 1 ) = -rateMatrix[category] ( 1, 0 );

        for ( int i = 0; i < 2; ++i ) {
            eigenValues[category] [i] = 0.0;
            junk[i] = 0.0;
            for ( int j = 0; j < 2; ++j ) {
                eigenMatrix[category] ( i, j ) = 0.0;
                ieigenMatrix[category] ( i, j ) = 0.0;
            }
        }
        // Find the eigenValues and eigenVectors of the rate matrix
        int err = LeftEigenSystem( rateMatrix[category], eigenMatrix[category],
        eigenValues[category], junk );
        // Find the inverse of the matrix whose rows are the left
        // eigen vectors of the rate matrix
        if (err){
            cerr << "Error during matrix exponentiation" << endl;
            exit(EXIT_FAILURE);
        }
        ieigenMatrix[category] = inverse( eigenMatrix[category] );
    } //for each category (except invariant)
}

Model * TWOSTATE::clone( ParametersSet & parameters ) const {
    if( parameters.findParameter( "Number of rates ratios sets" ) ){
        if ( parameters.intParameter("Number of rates ratios sets") != 0 ){
            cerr << "The \"Number of rates ratios sets\" parameter "
                 << "cannot be used with a 2-state model" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        parameters["Number of rates ratios sets"] = "0";
    }
    return new TWOSTATE( parameters );
}
