#include "Models/RNA16D.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <numeric>

using namespace std;

//register the model to the model factory with the name RNA16D
RNA16D RNA16D::prototype( "RNA16D" );

RNA16D::RNA16D( const string & registrationName ) :
RNA16REVEQ( registrationName ) {
//the private constructor is called by the prototype only
//RNA16D is related to RNA16REVEQ on some aspects
}

RNA16D::RNA16D( ParametersSet & parameters )
: RNA16REVEQ( parameters ) {
//the normal constructor is called by the prototype's clone method
//RNA16D is a restriction of RNA16REVEQ
}


void RNA16D::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //The initialisation procedure for the RNA16D model is slightly different
    //than for the 16REVEQ model, therefore we redefine it
    //initialise the number of parameters
    numberRatesRatios = 3;
    numberFrequencies = 4;
    
    //init the equivalency table (correspondance symbol/states)
    initEquivalencyTable();
    
    //initialise the matrixIndex for correspondance between a matrix element
    //and an index in the rate ratios vector
    initMatrixIndex();

    //call the basic initialisation method
    RnaModel::initialisation( sequenceTable, modelId );
    
    //if a sequence table is provided initialise rate ratios and frequencies
    //accordingly
    if ( sequenceTable ) {
        array2D<double> rates16;
        vector<double> freq16;
        vector<double> freq4;
        vector<double> rates3;
        freq4.resize(4);
        //alpha=transition=1.0; beta=transversion; lambda = pair<->MM factor; theta = GU behaviour (1=equivalent(MM), lambda=equivalent(watson-crick))
        rates3.resize(3);
        
        freq16 = retrieveEmpiricalFrequencies( sequenceTable, modelId );
        freq4[0] = freq16[0]+freq16[3]+freq16[7]+freq16[8]+freq16[9]+freq16[11] + 2*freq16[6];
        freq4[1] = freq16[2]+freq16[5]+freq16[8]+freq16[11]+freq16[13]+freq16[14] + 2*freq16[12];
        freq4[2] = freq16[1]+freq16[2]+freq16[4]+freq16[5]+freq16[7]+freq16[9] + 2*freq16[10];
        freq4[3] = freq16[0]+freq16[1]+freq16[3]+freq16[4]+freq16[13]+freq16[14] + 2*freq16[15];
        double sum = accumulate(freq4.begin(),freq4.end(),0.0);
        freq4[0]/=sum; freq4[1]/=sum; freq4[2]/=sum; freq4[3]/=sum;
        frequencies->initialisation(freq4);
        //estimation of all the rate ratios (and proportion of invariant sites)
        //if relevant from the empirical sequences
        if ( invariantCategory ) {
            rates16 = retrieveEmpiricalRates( sequenceTable, modelId, & proportionInvariantSites );
        }
        else {
            rates16 = retrieveEmpiricalRates( sequenceTable, modelId );
        }
        //fill the rate vector according to the estimation
        //reference ratio A<->G, it should be divided by 8
        //but since all the rates in rates5 are initialised from
        //8 rates in rates16 it is a useless computation
        double ref = rates16( 0, 1 ) + rates16( 2, 8 ) + rates16( 3, 4 ) +
                       rates16( 5, 11 )+ rates16( 6, 7 ) + rates16( 6, 9 ) +
                       rates16( 7, 10 )+ rates16( 9, 10 ) +
                       rates16( 0, 8 ) + rates16( 1, 2 ) + rates16( 3, 11 ) +
                       rates16( 4, 5 )+ rates16( 12, 13 ) + rates16( 12, 14 ) +
                       rates16( 13, 15 )+ rates16( 14, 15 );
        rates3[0] = ( rates16( 0, 13 ) + rates16( 2, 9 ) + rates16( 3, 14 ) +
                      rates16( 5, 7 )+ rates16( 6, 8 ) + rates16( 6, 11 ) +
                      rates16( 8, 12 )+ rates16( 11, 12 ) +
                      rates16( 0, 6 ) + rates16( 0, 15 ) + rates16( 1, 9 ) +
                      rates16( 3, 6 )+ rates16( 3, 15 ) + rates16( 4, 7 ) +
                      rates16( 8, 14 )+ rates16( 11, 13 ) +
                      rates16( 1, 13 ) + rates16( 2, 10 ) + rates16( 2, 12 ) +
                      rates16( 4, 14 )+ rates16( 5, 10 ) + rates16( 5, 12 ) +
                      rates16( 7, 8 )+ rates16( 9, 11 ) +      
                      rates16( 0, 7 ) + rates16( 1, 10 ) + rates16( 1, 15 ) +
                      rates16( 2, 14 )+ rates16( 3, 9 ) + rates16( 4, 10 ) +
                      rates16( 4, 15 )+ rates16( 5, 13 ) ) / (2.0 * ref);
        //estimate kappa from pi(MM) deduce lambda and theta
        double pairMM = freq16[6]+freq16[7]+freq16[8]+freq16[9]+freq16[10]+freq16[11]+freq16[12]+freq16[13]+freq16[14]+freq16[15];
        double unpairMM = 2*freq4[0]*(freq4[0] + freq4[1] + freq4[2]) + 2*freq4[1]*(freq4[1]+freq4[3]);
        double kappa_est=pairMM/unpairMM;
        //lambda_est^2*kappa*unpairWC=pairWC
        rates3[1] = sqrt( (freq16[0]+freq16[2]+freq16[3]+freq16[5])/
                          ( kappa_est * (2*freq4[0]*freq4[3]+2*freq4[1]*freq4[2]) ) );
        rates3[2] = sqrt( (freq16[1]+freq16[4]/
                          ( kappa_est * 2 * freq4[2]*freq4[3]) ) );
        ratesRatios->initialisation(rates3);
    }
    
    //specific model check...
    if ( invariantCategory ){
        if ( (ratesRatios->getNumberRatesRatiosSets() != 1) ||
                                (frequencies->getNumberFrequenciesSets() != 1) ){
            cerr << "ERROR: Frequency and exchangeabilities parameters are closely related in " << getName() << endl
                 << "Cannot use more than 1 set of frequencies/rate ratios with +I" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        if ( ratesRatios->getNumberRatesRatiosSets() != frequencies->getNumberFrequenciesSets() ){
            cerr << "ERROR: Frequency and exchangeabilities parameters are closely related in " << getName() << endl
                 << "One must have nb frequency sets == nb exchangeability sets" << endl;
            exit(EXIT_FAILURE);
        }
        else{
            for (unsigned int i=0; i < numberRatesCategories; ++i){
                if (frequencies->frequenciesCat(i)!=ratesRatios->ratesRatiosCat(i)){
                    cerr << "ERROR: frequencies set and rate ratios set assigned to the category " << i
                         << " are different: " << frequencies->frequenciesCat(i) << "<>" << ratesRatios->ratesRatiosCat(i) << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    
    //first initialisation of the substitution rate matrix
    //(factor(s) used with the subtitution matrix to reach the specified
    //average substitution rate) and eigen system
    updateAverageRateVector();
    updateEigenMatrix();
}


// Destructor
RNA16D::~RNA16D() {
}


// The models unique name
string RNA16D::getName( void ) const {
    return ( "RNA16D" + MatrixModel::getName() );
}


double RNA16D::getExchangeability( unsigned int residue1,
    unsigned int residue2, unsigned int gammaCategory ) const {
    assert ( residue1 < getNumberStates() );
    assert ( residue2 < getNumberStates() );
    int index = matrixIndex(residue1, residue2);
    //error: i==j
    assert (index != 0xFF);
    //null
    if (index == 0x0F){
        return 0.0;
    }
    //0x06 does not exist since it is for pair<->GU and it must be a transition: 0x06|0x01==0x07
    assert ((index & 0x0F) != 0x06);
    double exch = getInvKappa(gammaCategory+invariantCategory);
    //transversion
    if (!(index & 0x01)){
        exch *= (*ratesRatios)[gammaCategory][0];
    }
    //div by theta for pair<->GU and GU<->MM ?
    if (index & 0x02){
        exch /= (*ratesRatios)[gammaCategory][2];
    }
    //div by lambda for pair<->GU and pair<->MM ?
    if (index & 0x04){
        exch /= (*ratesRatios)[gammaCategory][1];
    }
    //div by common nucleotide
    //freqSet should be = ratesSet;
    return exch/(*frequencies)[gammaCategory+invariantCategory][((index & 0xF0)>>4)];
}


void RNA16D::initMatrixIndex() {
    matrixIndex.resize( 16, 16 );  
    //initialise the diagonal (-3)
    for ( int i = 0 ; i < 16 ; ++i ) {
        matrixIndex(i, i) = 0xFF;
    }
    //fill the matrix for double substitution (-1 for rate = 0)
    for ( int i = 0; i < 15 ; ++i ) {
        for ( int j = i+1; j < 16 ; ++j ) {
            matrixIndex(i, j) = 0x0F;
        }
    }
    //transitions
    matrixIndex( 0, 1 ) = matrixIndex(2, 8) =
        matrixIndex(3, 4) = matrixIndex(5, 11) =
        matrixIndex(6, 7) = matrixIndex(6, 9) =
        matrixIndex(7, 10) = matrixIndex(9, 10) = 0x01;
    matrixIndex( 0, 8 ) = matrixIndex( 1, 2 ) =
        matrixIndex( 3, 11 ) = matrixIndex( 4, 5 ) =
        matrixIndex( 12, 13 ) = matrixIndex( 12, 14 ) =
        matrixIndex( 13, 15 ) = matrixIndex( 14, 15 ) = 0x01;
    //transversions
    matrixIndex( 0, 13 ) = matrixIndex( 2, 9 ) =
        matrixIndex( 3, 14 ) = matrixIndex( 5, 7 ) =
        matrixIndex( 6, 8 ) = matrixIndex( 6, 11 ) =
        matrixIndex( 8, 12 ) = matrixIndex( 11, 12 ) = 0x00;
   matrixIndex( 0, 6 ) = matrixIndex( 0, 15 ) =
        matrixIndex( 1, 9 ) = matrixIndex( 3, 6 ) =
        matrixIndex( 3, 15 ) = matrixIndex( 4, 7 ) =
        matrixIndex( 8, 14 ) = matrixIndex( 11, 13 ) = 0x00;
   matrixIndex( 1, 13 ) = matrixIndex( 2, 10 ) =
        matrixIndex( 2, 12 ) = matrixIndex( 4, 14 ) =
        matrixIndex( 5, 10 ) = matrixIndex( 5, 12 ) =
        matrixIndex( 7, 8 ) = matrixIndex( 9, 11 ) = 0x00;
    matrixIndex( 0, 7 ) = matrixIndex( 1, 10 ) =
        matrixIndex( 1, 15 ) = matrixIndex( 2, 14 ) =
        matrixIndex( 3, 9 ) = matrixIndex( 4, 10 ) =
        matrixIndex( 4, 15 ) = matrixIndex( 5, 13 ) = 0x00;
        
    //MM<->MM
    //....nothing more, code 0x00
    
    //pair<->MM
    //trigger the division by lambda or thetha
    string pair1;
    string pair2;
    for ( int i = 0; i < 6; ++i ) {
        pair1 = getState(i);
        for ( int j = 6; j < 16; ++j ) {
            pair2 = getState(j);
            assert(pair1.size() == 2);
            assert(pair2.size() == 2);
            //if one nucleotide is the same rate = pair/mismatch rate ratio
            if ( (pair1[0]==pair2[0]) || (pair1[1]==pair2[1]) ){
                //GU<->MM
                if ( (i == 1) || (i==4) ){
                    matrixIndex(i, j) |= 0x02;
                }
                //watson-crick<->MM
                else{
                    matrixIndex(i, j) |= 0x04;
                }
            }
        }
    }
    //trigger both division (waston-crick<->GU)
    matrixIndex( 0, 1 ) |= 0x06;
    matrixIndex(1, 2) |= 0x06;
    matrixIndex(3, 4) |= 0x06;
    matrixIndex(4, 5) |= 0x06;
    
    
    //for the REV equivalent, common echangeability parameters
    //are divided by frequencies,
    //we squeeze the frequency in the actual values, the division
    //is handled by getExchangeability(..)
    //AA<->X exchangeabilities are divided by F(A)
    //matrixIndex(0,6) += 0;
    //matrixIndex(3,6) += 0;
    //matrixIndex(6,7) += 0;
    //matrixIndex(6,8) += 0;
    //matrixIndex(6,9) += 0;
    //matrixIndex(6,11) += 0;
    //C?<->C? exchangeabilities are divided by F(C)
    matrixIndex(2,12) = matrixIndex(2,12) | 0x10;
    matrixIndex(5,12) = matrixIndex(5,12) | 0x10;
    matrixIndex(8,12) = matrixIndex(8,12) | 0x10;
    matrixIndex(11,12) = matrixIndex(11,12) | 0x10;
    matrixIndex(12,13) = matrixIndex(12,13) | 0x10;
    matrixIndex(12,14) = matrixIndex(12,14) | 0x10;
    matrixIndex(5,13) = matrixIndex(5,13) | 0x10;
    matrixIndex(11,13) = matrixIndex(11,13) | 0x10;    
    matrixIndex(2,14) = matrixIndex(2,14) | 0x10;
    matrixIndex(8,14) = matrixIndex(8,14) | 0x10;
    matrixIndex(2,8) = matrixIndex(2,8) | 0x10;
    matrixIndex(5,11) = matrixIndex(5,11) | 0x10;
    
    //G?<->G? exchangeabilities are divided by F(G)
    matrixIndex(1,10) = matrixIndex(1,10) | 0x20;
    matrixIndex(2,10) = matrixIndex(2,10) | 0x20;
    matrixIndex(4,10) = matrixIndex(4,10) | 0x20;
    matrixIndex(5,10) = matrixIndex(5,10) | 0x20;
    matrixIndex(7,10) = matrixIndex(7,10) | 0x20;
    matrixIndex(9,10) = matrixIndex(9,10) | 0x20;
    matrixIndex(1,9) = matrixIndex(1,9) | 0x20;
    matrixIndex(2,9) = matrixIndex(2,9) | 0x20;
    matrixIndex(1,2) = matrixIndex(1,2) | 0x20;
    matrixIndex(4,7) = matrixIndex(4,7) | 0x20;
    matrixIndex(5,7) = matrixIndex(5,7) | 0x20;
    matrixIndex(4,5) = matrixIndex(4,5) | 0x20;
    
    
    //U?<->U? exchangeabilities are divided by F(U)
    matrixIndex(0,15) = matrixIndex(0,15) | 0x30;
    matrixIndex(1,15) = matrixIndex(1,15) | 0x30;
    matrixIndex(3,15) = matrixIndex(3,15) | 0x30;
    matrixIndex(4,15) = matrixIndex(4,15) | 0x30;
    matrixIndex(13,15) = matrixIndex(13,15) | 0x30;
    matrixIndex(14,15) = matrixIndex(14,15) | 0x30;
    matrixIndex(3,14) = matrixIndex(3,14) | 0x30;
    matrixIndex(4,14) = matrixIndex(4,14) | 0x30;
    matrixIndex(3,4) = matrixIndex(3,4) | 0x30;
    matrixIndex(0,13) = matrixIndex(0,13) | 0x30;
    matrixIndex(1,13) = matrixIndex(1,13) | 0x30;
    matrixIndex(0,1) = matrixIndex(0,1) | 0x30;
    
    //symetrize the matrix
    for ( int i = 0; i < 15 ; ++i ) {
        for ( int j = i+1; j < 16 ; ++j ) {
            matrixIndex(j, i) = matrixIndex(i, j);
        }
    }
}

double RNA16D::getFrequency( unsigned int residue,
        unsigned int rateCategory, unsigned int ) const {
    unsigned int mixtCat = rateCategory - invariantCategory;
    double lambda2 = (*ratesRatios)[mixtCat][1] * (*ratesRatios)[mixtCat][1];
    double theta2 = (*ratesRatios)[mixtCat][2] * (*ratesRatios)[mixtCat][2];
    double kappa = 1.0/getInvKappa(rateCategory);
    switch (residue){
        case 0 : return ( kappa*lambda2 * (*frequencies)[rateCategory][0] * (*frequencies)[rateCategory][3] ); break;
        case 1 : return ( kappa*theta2 * (*frequencies)[rateCategory][2] * (*frequencies)[rateCategory][3] ); break;
        case 2 : return ( kappa*lambda2 * (*frequencies)[rateCategory][2] * (*frequencies)[rateCategory][1] ); break;
        case 3 : return ( kappa*lambda2 * (*frequencies)[rateCategory][3] * (*frequencies)[rateCategory][0] ); break;
        case 4 : return ( kappa*theta2 * (*frequencies)[rateCategory][3] * (*frequencies)[rateCategory][2] ); break;
        case 5 : return ( kappa*lambda2 * (*frequencies)[rateCategory][1] * (*frequencies)[rateCategory][2] ); break;
        case 6 : return ( kappa * (*frequencies)[rateCategory][0] * (*frequencies)[rateCategory][0] ); break;
        case 7 : return ( kappa * (*frequencies)[rateCategory][0] * (*frequencies)[rateCategory][2] ); break;
        case 8 : return ( kappa * (*frequencies)[rateCategory][0] * (*frequencies)[rateCategory][1] ); break;
        case 9 : return ( kappa * (*frequencies)[rateCategory][2] * (*frequencies)[rateCategory][0] ); break;
        case 10 : return ( kappa * (*frequencies)[rateCategory][2] * (*frequencies)[rateCategory][2] ); break;
        case 11 : return ( kappa * (*frequencies)[rateCategory][1] * (*frequencies)[rateCategory][0] ); break;
        case 12 : return ( kappa * (*frequencies)[rateCategory][1] * (*frequencies)[rateCategory][1] ); break;
        case 13 : return ( kappa * (*frequencies)[rateCategory][1] * (*frequencies)[rateCategory][3] ); break;
        case 14 : return ( kappa * (*frequencies)[rateCategory][3] * (*frequencies)[rateCategory][1] ); break;
        case 15 : return ( kappa * (*frequencies)[rateCategory][3] * (*frequencies)[rateCategory][3] ); break;
        default: assert(0); return -10000.0; break;
    }
}

double RNA16D::getInvKappa( unsigned int rateCategory ) const {
    assert(rateCategory>=invariantCategory);
    unsigned int mixtCat = rateCategory - invariantCategory;
    double lambda2 = (*ratesRatios)[mixtCat][1] * (*ratesRatios)[mixtCat][1];
    double theta2 = (*ratesRatios)[mixtCat][2] * (*ratesRatios)[mixtCat][2];
    double invKappa = 2.0 * ( lambda2 - 1.0 );
    invKappa *= ( ((*frequencies)[rateCategory][0] * (*frequencies)[rateCategory][3]) +
               ((*frequencies)[rateCategory][1] * (*frequencies)[rateCategory][2]) );
    invKappa += 2.0 * (*frequencies)[rateCategory][2] * (*frequencies)[rateCategory][3] *
                   ( theta2 - 1.0 );           
    invKappa += 1.0;
    return invKappa;
}
