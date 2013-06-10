#include <assert.h>

#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Models/Perturbator.h"
#include "Models/T92.h"

#include "Sequence/SequenceTable.h"

#include "Util/statlib.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"

using namespace std;

//register the model to the model factory with the name REV
T92 T92::prototype( "T92" );

T92::T92( const string & registrationName ) : DnaModel( registrationName ){
//the private constructor is called by the prototype only
}


T92::T92( ParametersSet & parameters ) : DnaModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void T92::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 1;
    numberFrequencies = 2;
    
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
        //estimation of the model frequencies from the empirical frequencies
        vector<double> freq4 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        //G+C
        vector<double> freq2;
        freq2.push_back(freq4[0] + freq4[3]);
        freq2.push_back(freq4[1] + freq4[2]);
        frequencies->initialisation( freq2 );
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
        vector<double> rate1;
        rate1.push_back( (rates5[0]+rates5[1]+rates5[2]+rates5[4])/(2*(1.0+rates5[3])) );
        ratesRatios->initialisation( rate1 );
    }
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}  

T92::~T92(){
}

string T92::getName( void ) const {
    return ( "T92" + MatrixModel::getName() );
}

double T92::getFrequency( unsigned int residue,
    unsigned int rateCategory, unsigned int ) const {
    switch (residue){
        case 0: case 3:
            return (*frequencies)[rateCategory][0]/2.0;
        case 1: case 2:
            return (*frequencies)[rateCategory][1]/2.0;
    }
    assert("wrong residue in T92:getFrequency()"==0);
    return 0.0;
}

string T92::getFrequencyState( unsigned int freqState, unsigned int ) const{
    assert(freqState<2);
    if (freqState){
        return "C+G";
    }
    else{
        return "A+T";
    }
}

void T92::initMatrixIndex() {
    matrixIndex.resize( 4, 4 );
    for ( int i = 0 ; i < 4 ; ++i ) {
        matrixIndex(i, i) = -3;
    }
    //reference
    matrixIndex(0, 2) = matrixIndex(2, 0) = -2;
    //transversion
    matrixIndex(0, 1) = matrixIndex(1, 0) = 0;
    matrixIndex(0, 3) = matrixIndex(3, 0) = 0;
    matrixIndex(1, 2) = matrixIndex(2, 1) = 0;
    matrixIndex(2, 3) = matrixIndex(3, 2) = 0;
    //second transition
    matrixIndex(1, 3) = matrixIndex(3, 1) = -2;
}


void T92::setEigenMatrix() {  
    vector< double > junk( 4 );
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
    category++ ) {
        // fill the rateMatrix (cf REV.h)
        int freqCat = frequencies->frequenciesCat(category+invariantCategory);
        int rateCat = ratesRatios->ratesRatiosCat(category);
        rateMatrix[category] ( 0, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1] * .5 *
        (*ratesRatios)(rateCat) [0]; // A->C
        rateMatrix[category] ( 0, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1] * .5;
        rateMatrix[category] ( 0, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0] * .5 *
        (*ratesRatios)(rateCat) [1]; // A->U
        rateMatrix[category] ( 0, 0 ) =
        -( ( rateMatrix[category] ) ( 0, 1 ) +
        ( rateMatrix[category] ) ( 0, 2 ) +
        ( rateMatrix[category] ) ( 0, 3 ) );
        rateMatrix[category] ( 1, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0] * .5 *
        (*ratesRatios)(rateCat) [0]; // C->A
        rateMatrix[category] ( 1, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1] * .5 *
        (*ratesRatios)(rateCat) [2]; // C->G
        rateMatrix[category] ( 1, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0] * .5 *
        (*ratesRatios)(rateCat) [3]; // C->U
        rateMatrix[category] ( 1, 1 ) =
        -( ( rateMatrix[category] ) ( 1, 0 ) +
        ( rateMatrix[category] ) ( 1, 2 ) +
        ( rateMatrix[category] ) ( 1, 3 ) );
        rateMatrix[category] ( 2, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0] * .5; // G->A
        rateMatrix[category] ( 2, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[1] * .5 *
        (*ratesRatios)(rateCat) [2]; // G->C
        rateMatrix[category] ( 2, 3 ) =
        substitutionRate[category] * (*frequencies)(freqCat)[0] * .5 *
        (*ratesRatios)(rateCat) [4]; // G->U
        rateMatrix[category] ( 2, 2 ) =
        -( ( rateMatrix[category] ) ( 2, 0 ) +
        ( rateMatrix[category] ) ( 2, 1 ) +
        ( rateMatrix[category] ) ( 2, 3 ) );
        rateMatrix[category] ( 3, 0 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [0] * .5 *
        (*ratesRatios)(rateCat) [1]; // U->A
        rateMatrix[category] ( 3, 1 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [1] * .5 *
        (*ratesRatios)(rateCat) [3]; // U->C
        rateMatrix[category] ( 3, 2 ) =
        substitutionRate[category] * (*frequencies)(freqCat) [1] * .5 *
        (*ratesRatios)(rateCat) [4]; // U->G
        rateMatrix[category] ( 3, 3 ) =
        -( ( rateMatrix[category] ) ( 3, 0 ) +
        ( rateMatrix[category] ) ( 3, 1 ) +
        ( rateMatrix[category] ) ( 3, 2 ) );

        for ( int i = 0; i < 4; ++i ) {
            eigenValues[category] [i] = 0.0;
            junk[i] = 0.0;
            for ( int j = 0; j < 4; ++j ) {
                eigenMatrix[category] ( i, j ) = 0.0;
                ieigenMatrix[category] ( i, j ) = 0.0;
            }
        }
        
        double transitionRate = substitutionRate[category];
        double transversionRate = substitutionRate[category] * (*ratesRatios)(rateCat)[0];
        
        //double purineFreq = .5;
        //double pyrimidineFreq = .5;

        eigenValues[category][0] = 0.0;
        eigenValues[category][1] = - .5 *
                (transitionRate + transversionRate);
        eigenValues[category][2] = eigenValues[category] [1];
        eigenValues[category][3] = -transversionRate;


        // Update the eigenmatrix by hand
        eigenMatrix[category] ( 0, 0 ) = (*frequencies)(freqCat)[0]*.5;
        eigenMatrix[category] ( 0, 3 ) = eigenMatrix[category] ( 0, 0 );
        eigenMatrix[category] ( 0, 1 ) = (*frequencies)(freqCat)[1]*.5;
        eigenMatrix[category] ( 0, 2 ) = eigenMatrix[category] ( 0, 1 );

        eigenMatrix[category] ( 1, 0 ) = 0.0;
        eigenMatrix[category] ( 1, 1 ) = -1.0;
        eigenMatrix[category] ( 1, 2 ) = 0.0;
        eigenMatrix[category] ( 1, 3 ) = 1.0;

        eigenMatrix[category] ( 2, 0 ) = -1.0;
        eigenMatrix[category] ( 2, 1 ) = 0.0;
        eigenMatrix[category] ( 2, 2 ) = 1.0;
        eigenMatrix[category] ( 2, 3 ) = 0.0;

        eigenMatrix[category] ( 3, 0 ) = -(*frequencies)(freqCat)[0] * .5;
        eigenMatrix[category] ( 3, 1 ) = (*frequencies)(freqCat)[1] * .5;
        eigenMatrix[category] ( 3, 2 ) = -(*frequencies)(freqCat)[1] * .5;
        eigenMatrix[category] ( 3, 3 ) = (*frequencies)(freqCat)[0] * .5;

        // Inverse of the eigenmatrix
        ieigenMatrix[category] ( 0, 0 ) = 1.0;
        ieigenMatrix[category] ( 0, 1 ) = 0.0;
        ieigenMatrix[category] ( 0, 2 ) = -(*frequencies)(freqCat)[1];
        ieigenMatrix[category] ( 0, 3 ) = -1.0;

        ieigenMatrix[category] ( 1, 0 ) = 1.0;
        ieigenMatrix[category] ( 1, 1 ) = -(*frequencies)(freqCat)[0];
        ieigenMatrix[category] ( 1, 2 ) = 0.0;
        ieigenMatrix[category] ( 1, 3 ) = 1.0;

        ieigenMatrix[category] ( 2, 0 ) = 1.0;
        ieigenMatrix[category] ( 2, 1 ) = 0.0;
        ieigenMatrix[category] ( 2, 2 ) = (*frequencies)(freqCat)[0];
        ieigenMatrix[category] ( 2, 3 ) = -1.0;

        ieigenMatrix[category] ( 3, 0 ) = 1.0;
        ieigenMatrix[category] ( 3, 1 ) = (*frequencies)(freqCat)[1];
        ieigenMatrix[category] ( 3, 2 ) = 0.0;
        ieigenMatrix[category] ( 3, 3 ) = 1.0;
    } //for each category (except invariant)
}

