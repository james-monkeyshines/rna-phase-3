#include <assert.h>

#include <iostream>
#include <iomanip>
#include <stdio.h>

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Models/Perturbator.h"
#include "Models/YNH98.h"

#include "Sequence/SequenceTable.h"

#include "Util/array2D.h"
#include "Util/statlib.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"


#define SQR(n) (n*n)

using namespace std;

//register the model to the model factory with the name YNH98
YNH98 YNH98::prototype( "YNH98" );

YNH98::YNH98( const string & registrationName ) : CodonModel( registrationName ){
//the private constructor is called by the prototype only
}

YNH98::YNH98( ParametersSet & parameters ) : CodonModel( parameters ){
//the normal constructor is called by the prototype's clone method
}

void YNH98::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //initialise the number of parameters
    numberRatesRatios = 2; // omega and kappa
    numberFrequencies = matrixSize;

    //call the basic initialisation method
    CodonModel::initialisation( sequenceTable, modelId );

    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();

    //initialisation of the matrix which defines the model
    initDefModelMatrix();

    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if ( sequenceTable ) {
        array2D<double> arrayRates;
        vector<double> vecRates;

        if ( invariantCategory ) {
            arrayRates = retrieveEmpiricalRates( sequenceTable, modelId, & proportionInvariantSites );
        }
        else {
            arrayRates = retrieveEmpiricalRates( sequenceTable, modelId );
        }

        double sumTransitionRates = 0.0;
        double sumTransversionRates = 0.0;
        double sumSynonymous = 0.0;
        double sumNonSynonymous = 0.0;

        for (unsigned int i = 0; i < matrixSize; ++i){
            for (unsigned int j = 0; j < matrixSize; ++j){
                switch (defModelMatrix(i,j)){
                    case 1:
                        sumSynonymous += arrayRates(i,j);
                        sumTransversionRates += arrayRates(i,j);
                        break;
                    case 2:
                        sumNonSynonymous += arrayRates(i,j);
                        sumTransversionRates += arrayRates(i,j);
                        break;
                    case 3:
                        sumSynonymous += arrayRates(i,j);
                        sumTransitionRates += arrayRates(i,j);
                        break;
                    case 4:
                        sumNonSynonymous += arrayRates(i,j);
                        sumTransitionRates += arrayRates(i,j);
                        break;
                }
            }
        }

        // omega: ratio of nonsynonymous and synonymous subsitution rates
        vecRates.push_back(sumNonSynonymous/sumSynonymous);

        // kappa: ratio of transition and transversion
        vecRates.push_back( sumTransitionRates/sumTransversionRates );

        cout << "initial omega: " << vecRates[0] << " - initial kappa: " << vecRates[1] << endl;

        ratesRatios->initialisation( vecRates );

        frequencies->initialisation( retrieveEmpiricalFrequencies(sequenceTable, modelId) );
    }

    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}

YNH98::~YNH98() {
}

string YNH98::getName( void ) const {
    return ( "YNH98" + CodonModel::getName() );
}

int YNH98::isTransitionOrTransversion(char c1, char c2){
    switch (c2){
        case 'A': case 'C': case 'G': case 'U': break;
        default: cerr << "Unrecognized nucleotide in the genetic code" << endl;
                 exit(EXIT_FAILURE); return 0;
    }
    if (c1==c2){
        return 0;
    }
    switch (c1){
        case 'A': return ((c2=='G') ?  1 : -1);
        case 'G': return ((c2=='A') ?  1 : -1);
        case 'C': return ((c2=='U') ?  1 : -1);
        case 'U': return ((c2=='C') ?  1 : -1);
        default: cerr << "Unrecognized nucleotide in the genetic code" << endl;
                 exit(EXIT_FAILURE); return 0;
    }
}

void YNH98::initDefModelMatrix(){

    defModelMatrix.resize(matrixSize,matrixSize);

    for (unsigned int i = 0; i < matrixSize; ++i){
        for (unsigned int j = 0; j < matrixSize; ++j){
            if (i==j){
                defModelMatrix(i,j) = -1;
            }
            else{
                int nbTransition = 0;
                int nbTransversion = 0;
                string codon1 = getState(i,0);
                string codon2 = getState(j,0);

                for (unsigned int k = 0; k < 3; ++k){
                    int isTrans = isTransitionOrTransversion(codon1[k], codon2[k]);
                    switch (isTrans){
                        case -1: ++nbTransversion; break;
                        case  1: ++nbTransition; break;
                        default: break;
                    }
                }
                //if one and only one substitution
                if (nbTransition+nbTransversion==1){
                    bool isSynonymous = (geneticCode[i].second == geneticCode[j].second);
                    if (nbTransversion==1){
                        isSynonymous ? defModelMatrix(i,j) = 1 : defModelMatrix(i,j) = 2;
                    }
                    else{
                        isSynonymous ? defModelMatrix(i,j) = 3 : defModelMatrix(i,j) = 4;
                    }
                }
                //otherwise, instantaneous rate == 0
                else{
                    defModelMatrix(i,j) = 0;
                }
            }
        }
    }
}

double YNH98::getExchangeability( unsigned int residue1,
    unsigned int residue2, unsigned int gammaCategory ) const {
    assert ( residue1 < getNumberStates() );
    assert ( residue2 < getNumberStates() );

    switch (defModelMatrix(residue1,residue2)){
        case 0:
            return 0.0;
        case 1:
            return 1.0;
        case 2:
            return (*ratesRatios)[gammaCategory][0]; // omega
        case 3:
            return (*ratesRatios)[gammaCategory][1]; // kappa
        case 4:
            return (*ratesRatios)[gammaCategory][0] * (*ratesRatios)[gammaCategory][1]; // omega * kappa
        default:
            assert("Internal Error: YNH98::getExchangeability"==0); return 0.0;
    }
}
