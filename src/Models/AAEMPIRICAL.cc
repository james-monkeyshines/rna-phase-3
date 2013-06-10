#include <assert.h>

#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Models/Perturbator.h"
#include "Models/AAEMPIRICAL.h"

#include "Sequence/SequenceTable.h"

#include "Util/statlib.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"

using namespace std;

//register the model to the model factory with the name AAEMPIRICAL
AAEMPIRICAL AAEMPIRICAL::prototype( "AAEMPIRICAL" );

AAEMPIRICAL::AAEMPIRICAL( const string & registrationName ) : AaModel( registrationName ){
//the private constructor is called by the prototype only
}

AAEMPIRICAL::AAEMPIRICAL( ParametersSet & parameters ) : AaModel( parameters ){
    //the normal constructor is called by the prototype's clone method
    empiricalRatesFileName = parameters.stringParameter("Empirical values");
    // Attempt to open the file
    empiricalValuesFile.open( empiricalRatesFileName.c_str() );
    if ( !empiricalValuesFile.good() ) {
        cerr << "Unable to open file " << empiricalRatesFileName
             << " for reading, cannot retrieve empirical values..." << endl;
        exit(EXIT_FAILURE);
    }
        
    // initialise empirical rates from the data file
    readEmpiricalRates();
    
    plusF = false;
    
    if (parameters.findParameter("+F")){
        plusF = parameters.boolParameter("+F");
    }
    
    if(!plusF){
        readEmpiricalFrequencies();
    }
}

void AAEMPIRICAL::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 0;
    numberFrequencies = plusF ? 20 : 0;

    //call the basic initialisation method
    AaModel::initialisation( sequenceTable, modelId );

    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();

    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if (sequenceTable){
        //estimation of the model frequencies from the empirical frequencies
        if (plusF){
            vector<double> freq20 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
            frequencies->initialisation( freq20 );
        }
        else{
            // initialise empirical frequencies from the data file
            frequencies->initialisation(0);
        }
        ratesRatios->initialisation(0);
    }
    /*
    cout << endl;
    for (int i = 0; i < 20; ++i){
        cout << (*frequencies)[0][i] << "   " ;
    }
    cout << endl;
    */
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}



AAEMPIRICAL::~AAEMPIRICAL() {
}


string AAEMPIRICAL::getName( void ) const {
    string namePart;
    if ( plusF ){
        namePart = "+F";
    }
    return ( "AAEMPIRICAL(" + empiricalRatesFileName + ')' + namePart + MatrixModel::getName() );
}

double AAEMPIRICAL::getFrequency( unsigned int residue,
        unsigned int gammaCategory, unsigned int ) const {
    if(plusF){
        assert ( residue < (*frequencies)[gammaCategory].size() );
        return (*frequencies)[gammaCategory][residue];
    }
    else{
        assert((residue<getNumberStates()) && (!empiricalFrequencies.empty()));
        return empiricalFrequencies[residue];
    }
}

double AAEMPIRICAL::getExchangeability( unsigned int residue1,
    unsigned int residue2, unsigned int ) const {
    assert ( residue1 < getNumberStates() );
    assert ( residue2 < getNumberStates() );
    //error: i==j
    assert (residue1 != residue2);

    return empiricalRates(residue1, residue2);
}


void AAEMPIRICAL::readEmpiricalRates(){
    empiricalRates.resize(20, 20);

    char AaSymbols[20] = {'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};

    for (unsigned int i = 0; i < 20; ++i){
        char ch ;
        empiricalValuesFile >> ch;
        if (ch!=AaSymbols[i]){
            cerr << "The matrix raw values in the empirical file must be in that order:" << endl
                 << "A R N D C Q E G H I L K M F P S T W Y V" << endl
                 << "and lines must be LABELLED (this is a check to make sure the parsing goes ok...sorry)" << endl;
            exit(0);
        }
        for (unsigned int j = 0; j < i; ++j){
            empiricalValuesFile >> ws;
            empiricalValuesFile >> empiricalRates(i,j);
            empiricalRates(j,i) = empiricalRates(i,j);
        }
        empiricalRates(i,i) = -1.0;
    }
}

void AAEMPIRICAL::readEmpiricalFrequencies(){
    empiricalFrequencies.resize(20);
    double sum = 0.0;
    for (unsigned int i = 0; i < 20; ++i){
        empiricalValuesFile >> ws;
        empiricalValuesFile >> empiricalFrequencies[i];
        sum += empiricalFrequencies[i];
        if (empiricalFrequencies[i]==0.0){
            cerr << "Error while reading frequencies values in the empirical model file." << endl
                 << "Frequencies should not be null" << endl;
            exit(EXIT_FAILURE);
        }
    }
    if (fabs(sum-1.0)>.01){
        cerr << "Error, frequency values in the empirical model file do not sum to one." << endl;
        exit(EXIT_FAILURE);
    }
    for (unsigned int i = 0; i < 20; ++i){
        empiricalFrequencies[i]/=sum;
    }
}

void AAEMPIRICAL::setEigenMatrix() {
    for ( unsigned int category = 0; category < MAX( numberGammaCategories, 1 );
          ++category ) {

        // Initialise the rate matrix
        for ( unsigned int i = 0; i < 20 - 1; ++i ) {
            for ( unsigned int j = i + 1; j < 20; ++j ) {
                rateMatrix[category]( i, j ) =
                        getExchangeability( i, j, category ) *
                        getFrequency( j, category + invariantCategory) *
                        substitutionRate[category];
                rateMatrix[category] ( j, i ) =
                        getExchangeability( j, i, category ) *
                        getFrequency( i, category + invariantCategory) *
                        substitutionRate[category];
#ifdef DEBUG
                assert( !isnan(rateMatrix[category]( i, j )) );
                assert( !isnan(rateMatrix[category]( j, i )) );
#endif
            }
        }
        for ( unsigned int i = 0; i < 20; ++i ) {
            double sum = 0.0;
            for ( unsigned int j = 0; j < 20; ++j ) {
                if ( i != j ){
                    sum += rateMatrix[category] ( i, j );
                }
            }
            rateMatrix[category] ( i, i ) = -sum;
        }

        vector < double > junk( 20 );
        for ( unsigned int i = 0; i < 20; ++i ) {
            eigenValues[category] [i] = 0.0;
            junk[i] = 0.0;
            for ( unsigned int j = 0; j < 20; ++j ){
                eigenMatrix[category] ( i, j ) = 0.0;
                ieigenMatrix[category] ( i, j ) = 0.0;
            }
        }
        // Find the eigenvalues and eigenvectors of the rate matrix
        int err = LeftEigenSystem( rateMatrix[category], eigenMatrix[category],
        eigenValues[category], junk );
        if (err){
            cerr << "Error during matrix exponentiation" << endl;
            exit(EXIT_FAILURE);
        }
        // Find the inverse of the matrix whose rows are the left eigen
        // vectors of the rate matrix
        ieigenMatrix[category] = inverse( eigenMatrix[category] );
    }
}


Model * AAEMPIRICAL::clone( ParametersSet & parameters ) const {
    if( parameters.findParameter( "Number of rates ratios sets" ) ){
        if ( parameters.intParameter("Number of rates ratios sets") != 0 ){
            cerr << "The \"Number of rates ratios sets\" parameter "
                 << "cannot be used with a AAEMPIRICAL model" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        parameters["Number of rates ratios sets"] = "0";
    }
    if( parameters.findParameter( "Number of frequencies sets" ) ){
        if ( !parameters.findParameter("+F") ||
             !parameters.boolParameter("+F") ){
            if ( parameters.intParameter("Number of frequencies sets") != 0 ){
                cerr << "The \"Number of frequencies sets\" parameter "
                     << "cannot be used with a AAEMPIRICAL model which is not +F" << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    else{
        if ( !parameters.findParameter("+F") ||
             !parameters.boolParameter("+F") ){
            parameters["Number of frequencies sets"] = "0";
        }
    }
    return new AAEMPIRICAL( parameters );
}



