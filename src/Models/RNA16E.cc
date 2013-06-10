#include "Models/RNA16E.h"

#include "Util/statlib.h"

#include <assert.h>
#include <iostream>
#include <numeric>

using namespace std;

//register the model to the model factory with the name RNA16E
RNA16E RNA16E::prototype( "RNA16E" );

RNA16E::RNA16E( const string & registrationName ) :
RNA16D( registrationName ) {
//the private constructor is called by the prototype only
//RNA16E is simplified from RNA16D
}

RNA16E::RNA16E( ParametersSet & parameters )
: RNA16D( parameters ) {
//the normal constructor is called by the prototype's clone method
//RNA16E is a restriction of RNA16D
}


void RNA16E::initialisation( SequenceTable * sequenceTable, int modelId ) {
    //The initialisation procedure for the RNA16E model is slightly different
    //than for the 16D model, therefore we redefine it
    //initialise the number of parameters
    numberRatesRatios = 2;
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
        vector<double> rates2;
        freq4.resize(4);
        //alpha=transition=1.0; beta=transversion; lambda = pair<->MM factor; GU behaviour (=1=equivalent(MM))
        rates2.resize(2);
        
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
        rates2[0] = ( rates16( 0, 13 ) + rates16( 2, 9 ) + rates16( 3, 14 ) +
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
        //estimate kappa from pi(MM) deduce lambda
        double pairMM = freq16[6]+freq16[7]+freq16[8]+freq16[9]+freq16[10]+freq16[11]+freq16[12]+freq16[13]+freq16[14]+freq16[15]+freq16[1]+freq16[4];
        double unpairMM = 2*freq4[0]*(freq4[0] + freq4[1] + freq4[2]) + 2*freq4[1]*(freq4[1]+freq4[3]) +  2*freq4[2]*freq4[3];
        double kappa_est=pairMM/unpairMM;
        //lambda_est^2*kappa*unpairWC=pairWC
        rates2[1] = sqrt( (freq16[0]+freq16[2]+freq16[3]+freq16[5])/
                          ( kappa_est * (2*freq4[0]*freq4[3]+2*freq4[1]*freq4[2]) ) );
        ratesRatios->initialisation(rates2);
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
RNA16E::~RNA16E() {
}


// The models unique name
string RNA16E::getName( void ) const {
    return ( "RNA16E" + MatrixModel::getName() );
}


double RNA16E::getExchangeability( unsigned int residue1,
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
    double exch = getInvKappa(gammaCategory+invariantCategory);
    //transversion
    if (!(index & 0x01)){
        exch *= (*ratesRatios)[gammaCategory][0];
    }
    //div by lambda for pair<->GU and pair<->MM ?
    if (index & 0x04){
        exch /= (*ratesRatios)[gammaCategory][1];
    }
    //other
    return exch/(*frequencies)[gammaCategory+invariantCategory][((index & 0xF0)>>4)];
}
    


double RNA16E::getFrequency( unsigned int residue,
        unsigned int rateCategory, unsigned int ) const {
    unsigned int mixtCat = rateCategory - invariantCategory;
    double lambda2 = (*ratesRatios)[mixtCat][1] * (*ratesRatios)[mixtCat][1];
    double kappa = 1.0/getInvKappa(rateCategory);
    switch (residue){
        case 0 : return ( kappa*lambda2 * (*frequencies)[rateCategory][0] * (*frequencies)[rateCategory][3] ); break;
        case 1 : return ( kappa * (*frequencies)[rateCategory][2] * (*frequencies)[rateCategory][3] ); break;
        case 2 : return ( kappa*lambda2 * (*frequencies)[rateCategory][2] * (*frequencies)[rateCategory][1] ); break;
        case 3 : return ( kappa*lambda2 * (*frequencies)[rateCategory][3] * (*frequencies)[rateCategory][0] ); break;
        case 4 : return ( kappa * (*frequencies)[rateCategory][3] * (*frequencies)[rateCategory][2] ); break;
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

double RNA16E::getInvKappa( unsigned int rateCategory ) const {
    assert(rateCategory>=invariantCategory);
    unsigned int mixtCat = rateCategory - invariantCategory;
    double lambda2 = (*ratesRatios)[mixtCat][1] * (*ratesRatios)[mixtCat][1];
    double invKappa = 2.0 * ( lambda2 - 1.0 );
    invKappa *= ( ((*frequencies)[rateCategory][0] * (*frequencies)[rateCategory][3]) +
               ((*frequencies)[rateCategory][1] * (*frequencies)[rateCategory][2]) );
    invKappa += 1.0;
    return invKappa;
}
