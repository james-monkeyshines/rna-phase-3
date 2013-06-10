#include "Models/MatrixModel.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

#include "Sequence/SequenceTable.h"

#include "Models/Perturbator.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/statlib.h"
#include "Util/array2D.h"
#include "Util/randombox.h"
#include "Util/ParametersSet.h"

#define SQR(n) (n*n)

using namespace std;

MatrixModel::MatrixModel(){
}

MatrixModel::MatrixModel( const string & registrationName ) {
    Singleton < Factory<Model> > & modelFactory = Singleton < Factory<Model> >::instance();
    modelFactory.subscribe( this, registrationName );
}

MatrixModel::MatrixModel( ParametersSet & parameters ) {

    perturbator = NULL;

    initInv( parameters );
    initGamma( parameters );
    numberRatesCategories = ( numberGammaCategories ? numberGammaCategories : 1 ) + invariantCategory;
    ratesRatios = new RatesRatios( parameters, this );
    frequencies = new Frequencies( parameters, this );

    initLink( parameters );
    
}


void MatrixModel::initInv( ParametersSet & parameters ){
    // invariant category, the invariant will be computed later if a
    // sequenceTable is used to initialise the Model
    if ( ( parameters.findParameter( "Invariant sites" ) ) &&
         ( parameters.boolParameter( "Invariant sites" ) ) ) {
        invariantCategory = 1;
        proportionInvariantSites = 0.2;
    }
    else {
        invariantCategory = 0;
        proportionInvariantSites = 0.0;
    }
}

void MatrixModel::initGamma( ParametersSet & parameters ){
    //discrete gamma categories
    if ( parameters.findParameter( "Discrete gamma distribution of rates" ) ) {
        discreteGamma = parameters.boolParameter( "Discrete gamma distribution of rates" );
    }
    else{
        discreteGamma = false;
    }
    // initialise number gamma categories
    if (discreteGamma){
        numberGammaCategories = parameters.intParameter( "Number of gamma categories" );
        alpha = 0.1;
        if ( numberGammaCategories < 2 ){
            cerr << "You must use a least 2 gamma categories, change your"
                 << "\"Number of gamma categories\" parameter" << endl;
            exit(EXIT_FAILURE);
        }
        // space allocation
        substitutionRate.resize( numberGammaCategories );
        averageRate.resize( numberGammaCategories );
        categoriesWeights.resize( numberGammaCategories );
        if ( (parameters.findParameter("Laguerre quadrature") ) &&
             ( parameters.boolParameter("Laguerre quadrature") ) ){
            percentileGammaApproximation = false;
        }
        else{
            percentileGammaApproximation = true;
            fill( categoriesWeights.begin(), categoriesWeights.end(),
                    1.0/numberGammaCategories );
        }
    }
    else{
        numberGammaCategories = 0;
        // space allocation
        substitutionRate.resize( 1 );    
        averageRate.resize( 1 );    
    }
}

void MatrixModel::initMatrix( unsigned int matrixSize ){
    if ( discreteGamma ){
        eigenMatrix = new array2D < double > [numberGammaCategories];
        ieigenMatrix = new array2D < double > [numberGammaCategories];
        eigenValues = new vector < double > [numberGammaCategories];
        rateMatrix = new array2D < double > [numberGammaCategories];
        for (unsigned int i = 0; i < numberGammaCategories; ++i){
            eigenMatrix[i].resize(matrixSize, matrixSize);
            ieigenMatrix[i].resize(matrixSize, matrixSize);
            eigenValues[i].resize(matrixSize);
            rateMatrix[i].resize(matrixSize, matrixSize);
        }
    }
    else{
        eigenMatrix = new array2D < double >;
        ieigenMatrix = new array2D < double >;
        eigenValues = new vector < double >;
        rateMatrix = new array2D < double >;
        eigenMatrix->resize(matrixSize, matrixSize);
        ieigenMatrix->resize(matrixSize, matrixSize);
        eigenValues->resize(matrixSize);
        rateMatrix->resize(matrixSize, matrixSize);
   }    
}
        
        
void MatrixModel::initLink( ParametersSet & parameters ){
    if ( parameters.findParameter( "Link S parameter" ) ){
        if ( (frequencies->getGammaCatAdjust() & Frequencies::LS) &&
             (ratesRatios->getGammaCatAdjust() & RatesRatios::LS) ){
            if ( parameters.boolParameter("Link S parameter") ){
                double* sparam = frequencies->linkS();
                ratesRatios->assignS( sparam );
            }
        }
        else{
            cerr << "\"Link S parameter\" option not appropriate for your model..." << endl;
            exit(EXIT_FAILURE);
        }
    }
    
    if ( parameters.findParameter( "Link X parameter" ) ){
        if ( (frequencies->getGammaCatAdjust() & Frequencies::LX) &&
             (ratesRatios->getGammaCatAdjust() & RatesRatios::LX) ){
            if ( parameters.boolParameter("Link X parameter") ){
                vector<double> * sparam = frequencies->linkX();
                ratesRatios->assignX( sparam );
            }
        }
        else{
            cerr << "\"Link X parameter\" option not appropriate for your model..." << endl;
            exit(EXIT_FAILURE);
        }
    }
}



MatrixModel::~MatrixModel() {
    if (perturbator!=NULL){
        delete perturbator;
    }
    if( frequencies ){
        delete frequencies;
    }
    if( ratesRatios ){
        delete ratesRatios;
    }
    if ( discreteGamma ){
        if ( eigenMatrix != NULL ){
            delete[] eigenMatrix;
            eigenMatrix = NULL;
        }
        if ( ieigenMatrix != NULL ){
            delete[] ieigenMatrix;
            ieigenMatrix = NULL;
        }
        if ( eigenValues != NULL ){
            delete[] eigenValues;
            eigenValues = NULL;
        }
        if ( rateMatrix != NULL ){
            delete[] rateMatrix;
            rateMatrix = NULL;
        }
    }
    else{
        if ( eigenMatrix != NULL ){
            delete eigenMatrix;
            eigenMatrix = NULL;
        }
        if ( ieigenMatrix != NULL ){
            delete ieigenMatrix;
            ieigenMatrix = NULL;
        }
        if ( eigenValues != NULL ){
            delete eigenValues;
            eigenValues = NULL;
        }
        if ( rateMatrix != NULL ){
            delete rateMatrix;
            rateMatrix = NULL;
        }
    }
}

void MatrixModel::initialisation( SequenceTable * sequenceTable, int modelId ) {
    modelId = modelId;
    //frequencies and rateRatios must have the good size
    if ( sequenceTable == NULL ) {
        frequencies->initialisation( numberFrequencies );
        ratesRatios->initialisation( numberRatesRatios );
    }
    else{
        //in this case each model must initialise itself, frequencies and
        //ratesRatios from the sequenceTable
    }
    averageSubstitutionRate = 1.0;
}


ParametersSet MatrixModel::getModelParameters() const{
  
    ParametersSet parameters("MODEL");

    char tmpString[40];
    char dblString[40];
    
    parameters["Model name"] = getName();
    frequencies->getModelParameters( parameters, this, 0 );
    ratesRatios->getModelParameters( parameters );
    
    parameters("RATESCATEGORIES")["Invariant category"] =
        (invariantCategory ? "yes" : "no");
    if (invariantCategory){
        sprintf(dblString, "%.8f",proportionInvariantSites);
        parameters("RATESCATEGORIES")["Proportion of invariant sites"] =
            dblString;
    }
    parameters("RATESCATEGORIES")["Gamma distribution of rates"] =
        (discreteGamma ? "yes" : "no");
    if ( discreteGamma ){
        sprintf(tmpString,"%d",numberGammaCategories);
        parameters("RATESCATEGORIES")["Number of gamma categories"] =
            tmpString;
        sprintf(dblString, "%.8f",alpha);
        parameters("RATESCATEGORIES")["Alpha parameter"] = dblString;
    }
    return parameters;
}


void MatrixModel::setModelParameters( ParametersSet& parameters ){

    if (getName() != parameters.stringParameter("Model name")){
        cerr << "You are trying to initialise the model " << getName() << endl;
        cerr << "with parameters for the model "
             << parameters.stringParameter("Model name") << endl;
        cerr << "Exit..." << endl;
        exit(EXIT_FAILURE);
    }
    frequencies->setModelParameters( parameters, this, 0 );
    ratesRatios->setModelParameters( parameters );
    
    if ( (invariantCategory == 1) !=
        parameters("RATESCATEGORIES").boolParameter("Invariant category") ){
        cerr << "The name of your model "
             << parameters.stringParameter("Model name") << endl;
        cerr << "does not match with the \"Invariant category\" field "
             << "in your model file" << endl;
        cerr << "Exit..." << endl;
        exit(EXIT_FAILURE);        
    }
    if (invariantCategory){
        proportionInvariantSites =
            parameters("RATESCATEGORIES").doubleParameter("Proportion of invariant sites");
    }
    if( discreteGamma !=
        parameters("RATESCATEGORIES").boolParameter("Gamma distribution of rates") ){
        cerr << "The name of your model "
             << parameters.stringParameter("Model name") << endl;
        cerr << "does not match with the \"Gamma distribution of rates\" field "
             << "in your model file" << endl;
        cerr << "Exit..." << endl;
        exit(EXIT_FAILURE);        
    }
    if ( discreteGamma ){
        if ( numberGammaCategories != (unsigned int)
            parameters("RATESCATEGORIES").intParameter("Number of gamma categories") ){
            cerr << "The name of your model "
                 << parameters.stringParameter("Model name") << endl;
            cerr << "does not match with the \"Number of gamma categories\" field "
                 << "in your model file" << endl;
            cerr << "Exit..." << endl;
            exit(EXIT_FAILURE);        
        }
        alpha = parameters("RATESCATEGORIES").doubleParameter("Alpha parameter");
    }
}


double MatrixModel::getExchangeability( unsigned int residue1,
    unsigned int residue2, unsigned int gammaCategory ) const {
    assert ( residue1 < getNumberStates() );
    assert ( residue2 < getNumberStates() );
    int index = matrixIndex(residue1, residue2);
    //error: i==j
    assert (index != -3);
    //reference
    if (index == -2){
        return 1.0;
    }
    //null
    if (index == -1){
        return 0.0;
    }
    return (*ratesRatios)[gammaCategory][index];
}
    

double MatrixModel::getAverageSubstitutionRate( unsigned int ) const {
    return averageSubstitutionRate;
}


void MatrixModel::getAllParameters( vector < double > & pars ) const {
    vector< double > p;
    
    pars.clear(); 
    frequencies->getAllParameters( p );
    pars.insert( pars.end(), p.begin(), p.end() );
    ratesRatios->getAllParameters( p );
    pars.insert( pars.end(), p.begin(), p.end() );
    
    if ( discreteGamma ) {
        pars.push_back( sqrt( alpha ) );
    }
    if ( invariantCategory == 1 ) {
        pars.push_back( sqrt( proportionInvariantSites /
                              ( 1.0 - proportionInvariantSites ) ) );
    }
    assert( pars.size() == getNumberFreeParameters() );
}


void MatrixModel::setAllParameters( const vector < double > & pars ) {

    assert( pars.size() == getNumberFreeParameters() );

    vector < double >::const_iterator iter1;
    vector < double >::const_iterator iter2;
    
    iter1 = pars.begin() + frequencies->getNumberFreeParameters();
    frequencies->setAllParameters( vector<double>( pars.begin(), iter1 ) );
    iter2 = iter1 + ratesRatios->getNumberFreeParameters();
    ratesRatios->setAllParameters( vector<double>( iter1, iter2 ) );
        
    if ( discreteGamma ) {
        alpha = SQR(*iter2);
        ++iter2;
    }
    if ( invariantCategory == 1 ) {
        proportionInvariantSites =
        SQR(*iter2) / ( 1.0 + SQR( *iter2 ) );
        ++iter2;
    }
    assert(iter2==pars.end());
}


unsigned int MatrixModel::getNumberFreeParameters() const{
    return frequencies->getNumberFreeParameters() +
           ratesRatios->getNumberFreeParameters() +
           ( discreteGamma ? 1 : 0 ) + (unsigned int)invariantCategory;
}

unsigned int MatrixModel::getNumberFreeFrequencyParameters() const{
    return frequencies->getNumberFreeParameters();
}

void MatrixModel::getOptimisableParameters( bool empiricalFreqs, vector<unsigned int> &optimisableParameters ) const{
    unsigned int nbModelParameters = getNumberFreeParameters();
	
	if ( empiricalFreqs ) {
    	unsigned int nbFreqParameters = frequencies->getNumberFreeParameters();
		optimisableParameters.resize( nbModelParameters - nbFreqParameters );
		// Frequency parameters are always the first ones, so just need to skip those.
		unsigned int j = 0;
        for ( unsigned int i = nbFreqParameters; i < nbModelParameters; ++i ){
			optimisableParameters.at(j++) = i;
		}
	}
	else{
		optimisableParameters.resize(nbModelParameters);
        for ( unsigned int i = 0; i < nbModelParameters; ++i ){
			optimisableParameters.at(i) = i;
		}
	}
}

double MatrixModel::getRateCategoryProbability( unsigned int category, unsigned int ) const {
    if ( invariantCategory == 1 ) {
        if ( category == 0 ) {
            return proportionInvariantSites;
        }
        assert ( ( category > 0 ) && ( category < getNumberRatesCategories() ) );
        if ( discreteGamma ) {
            if (percentileGammaApproximation){
                return ( 1.0 - proportionInvariantSites ) /
                       ( ( double )numberGammaCategories );
            }
            else{
                  return ( 1.0 - proportionInvariantSites ) * categoriesWeights[category-1];
            }
        }
        else {
            return ( 1.0 - proportionInvariantSites );
        }
    }
    // no invariant
    else {
        if ( discreteGamma ) {
            assert( category < numberGammaCategories );
            if (percentileGammaApproximation){
                return ( 1.0 / ( double )numberGammaCategories );
            }
            else{
                  return categoriesWeights[category];
            }
        }
        else {
            assert( category == 0 );
            return ( 1.0 );
        }
    }
}

void MatrixModel::getAllPenaltyParameters( vector<double>& params ) const{
    frequencies->getAllPenaltyParameters( params );
}
    
void MatrixModel::setAllPenaltyParameters( const vector<double>& params){
    frequencies->setAllPenaltyParameters( params );
}

unsigned int MatrixModel::getNumberPenaltyParameters() const{
    return frequencies->getNumberPenaltyParameters();
}

void MatrixModel::initialiseML( ParametersSet & parameters ){
    frequencies->initialiseML( parameters );
    //ratesRatios->initialiseML( parameters );
    updateAverageRateVector();
}
    
double MatrixModel::getLnPenalty() const{
    return frequencies->getLnPenalty( averageRate );
    //      + ratesRatios->getLnPenalty();
}

void MatrixModel::diffLnPenalty(vector<double>& gradVector) const{
    frequencies->diffLnPenalty( gradVector );
    //      + ratesRatios->diffLnPenalty();
}

void MatrixModel::initialiseMCMC( ParametersSet & parameters ) {

    perturbator = new Perturbator();
    
    frequencies->initialiseMCMC( parameters, perturbator );
    frequencies->attach( *this );
    ratesRatios->initialiseMCMC( parameters, perturbator );
    ratesRatios->attach( *this );
    
    if ( invariantCategory ) {
        if ( !parameters.findParameter("Invariant parameter, prior") ){
            parameters["Invariant parameter, prior"] = "uniform(0,1)";
        }
        perturbator->registerGaussParameter( *this, UpdateMessage::INVARIANT, "Invariant",
                              &proportionInvariantSites,
                              parameters, "Invariant parameter",
                              .05, 0.001, 0.5, "initial step" );
    }
    if ( discreteGamma ) {
        if ( !parameters.findParameter("Gamma parameter, prior") ){
            parameters["Gamma parameter, prior"] = "uniform(0,1000)";
        }
        perturbator->registerGaussParameter( *this, UpdateMessage::GAMMA, "Gamma", &alpha,
                              parameters, "Gamma parameter",
                              .2, 0.001, 5.0, "initial step" );
    }
    
    //trigger a call to frequencies->prepare() & ratesRatios->prepare()
    //necessary to save the averageRate vector if prior is VAR
    updateAverageRateVector();
}



void MatrixModel::printParameters( ostream & outputStream ) const{
    outputStream.setf(ios::fixed);
    outputStream << setprecision( 4 );
    outputStream << "Average Substitution Rate = " << averageSubstitutionRate << endl;
    if ( discreteGamma ) {
        outputStream << "Gamma shape parameter = " << alpha << endl;
    }
    if ( invariantCategory == 1 ) {
        outputStream << "Proportion of invariant sites = "
                     << proportionInvariantSites << endl;
    }
    if ( discreteGamma ) {
        vector<double> subRate;
        subRate.resize( numberGammaCategories );
        vector<double> weights;
        weights.resize( numberGammaCategories );
        double mr = averageSubstitutionRate / ( 1.0 - proportionInvariantSites );
        double beta = alpha / mr;
        if (percentileGammaApproximation){
            statlib::gamma_percentile_means( subRate, alpha, beta,
                               numberGammaCategories );
            fill( weights.begin(), weights.end(), 1.0/ numberGammaCategories );
        }
        else{
            statlib::gammaLaguerre( subRate, weights, alpha, beta,
                               numberGammaCategories );
        }
        for (unsigned int i = 0; i < numberGammaCategories; ++i){
            outputStream << "average substitution rate, gamma category " << i + 1 << ": "
                 << subRate[i] << ", weight: " << (invariantCategory ? (1.0-proportionInvariantSites) : 1.0 )
                                      * weights[i] << endl;
                    
        }
    }
    frequencies->printParameters( outputStream );
    ratesRatios->printParameters( outputStream );
    outputStream.unsetf(ios::fixed);
}


void MatrixModel::outputFreqVector( const vector<string>& labels, const vector<double>& values,
                           ostream & outputStream, unsigned int pad, unsigned int prec ) const{
    assert(labels.size()==values.size());
    outputStream << setprecision(prec);
    for ( unsigned int i = 0; i < labels.size(); ++i ){
        outputStream << "f[" << labels[i] << "]=" << values[i] << setw(pad-labels[i].size()-6-prec) << " ";
    }
    outputStream << endl;
}

void MatrixModel::outputMatrix( const vector<string>& labels, const array2D<double>& values,
                           ostream & outputStream, unsigned int pad, unsigned int prec) const{
    assert(labels.size()==values.numberColumns());
    outputStream << setprecision(prec);
    unsigned int ent=(pad*2)/3;
    for ( unsigned int i = 0; i < labels.size(); ++i ) {
        unsigned int sp2=(pad-labels[i].size())/2;
        unsigned int sp1=pad-sp2;
        if (i==0) outputStream << setw(ent) << " ";
        outputStream << setw(sp1) << labels[i] << setw(sp2) << " ";
    }
    outputStream << endl;
    for ( unsigned int i = 0; i < getNumberStates(); ++i ) {
        unsigned int sp2= ent - 1 -labels[i].size();
        outputStream << " " << labels[i] << setw(sp2) << " ";
        sp2=(pad-prec-3)/2;
        unsigned int sp1=pad-sp2;
        for ( unsigned int j = 0; j < getNumberStates(); ++j ) {
            if (i == j){
                outputStream << setw(sp1) << string(prec+2,'-') << setw(sp2) << " ";
            }
            else{
                unsigned int bef2=(pad-prec-3)/2;
                unsigned int bef1=pad-bef2;
                outputStream << setw(bef1) << values(i,j) << setw(bef2) << " ";
            }
        }
        outputStream << endl; //end line
    }
    outputStream << endl; //end rate matrix
}

void MatrixModel::printParametersAux( ostream & outputStream ) const{
  
    //IMPORTANT:    Duplicated in RNA16 => update conjointly...
    //IMPORTANT:    Skipped by some models (JC69, TwoStateModel, ...?) 
        
    //first output the sets
    outputStream.setf(ios::fixed);
    outputStream << setprecision(5);

    unsigned int nbStates = getNumberStates();
    
    vector<string> freqLabels;
    if (numberFrequencies){
        freqLabels.resize(numberFrequencies);
        vector<double> freqValues;
        freqValues.resize(numberFrequencies);
        for ( unsigned int i = 0; i < numberFrequencies; ++i ){
            freqLabels[i] = getFrequencyState( i );
        }
    
        for ( unsigned int freqCat = 0; freqCat < frequencies->getNumberFrequenciesSets(); ++freqCat ){
            if (frequencies->getNumberFrequenciesSets()>=2){
                outputStream << "set of frequencies " << freqCat + 1 << ":" << endl;
            }
            else{
                outputStream << "frequencies:" << endl;
            }
            for ( unsigned int i = 0; i < numberFrequencies; ++i ){
                freqValues[i] = (*frequencies)(freqCat)[i];
            }
            outputFreqVector( freqLabels, freqValues, outputStream, 12, 5 );
        }
    }
    
    if (numberRatesRatios){
        for ( unsigned int rateCat = 0; rateCat < ratesRatios->getNumberRatesRatiosSets(); ++rateCat ){
            if (ratesRatios->getNumberRatesRatiosSets()>=2){
                outputStream << "set of rate parameters " << rateCat + 1 << ":" << endl;
            }
            else{
                outputStream << "rate parameters:" << endl;
            }
            for ( unsigned int i = 0; i < numberRatesRatios; ++i ){
                outputStream << "r" << i+1 << "=" << (*ratesRatios)(rateCat)[i] << "  ";
            }
            outputStream << endl;
        }
        outputStream << endl;
    }
        
    //now get the equFreqs and excheangeabilities for each gamma category (+inv?)
    bool onlyOneEquFreq=true;
    vector< vector<double> > equFreqsMixture;
    equFreqsMixture.resize( MAX(numberGammaCategories,1) + invariantCategory );
    for( unsigned int mixtureCat = 0; mixtureCat < equFreqsMixture.size(); ++mixtureCat ){
        equFreqsMixture[mixtureCat].resize(nbStates);
        for (unsigned int i = 0; i < nbStates; ++i){
            equFreqsMixture[mixtureCat][i] = getFrequency( i, mixtureCat );
        }
        //turn onlyOneEquFreq to false if a difference is spotted
        if(onlyOneEquFreq&&mixtureCat){
            onlyOneEquFreq = ( equFreqsMixture[mixtureCat] == equFreqsMixture[mixtureCat-1] );
        }
    }

    
    bool onlyOneExchMat=true;
    vector< array2D<double> > exchMixture;    
    exchMixture.resize( MAX(numberGammaCategories,1) );
    for( unsigned int mixtureCat = 0; mixtureCat < exchMixture.size(); ++mixtureCat ){
        exchMixture[mixtureCat].resize( nbStates, nbStates );
        for (unsigned int i = 0; i < nbStates; ++i){
            for (unsigned int j = 0; j < nbStates; ++j){
                if (i!=j){
                    exchMixture[mixtureCat](i,j) = getExchangeability( i, j, mixtureCat );
                }
                else{
                    exchMixture[mixtureCat](i,j) = -20000.0;
                }
            }
        }
        //turn onlyOneExchMat to false if a difference is spotted
        if(onlyOneExchMat&&mixtureCat){
            onlyOneExchMat = ( exchMixture[mixtureCat] == exchMixture[mixtureCat-1] );
        }
    }


    vector<string> labels;
    labels.resize(nbStates);
    for (unsigned int i = 0; i < labels.size(); ++i){
        labels[i]=getState(i);
    }
    
    //in case frequencies and state frequencies are "different parameters"
    //(eg Galtier-Gouy98) or RNA 6-state models with 4 frequencies, turn on "different"
    bool different = false;
    if(numberFrequencies){
        different = (labels.size()!=freqLabels.size());
        //also check that names are the same
        for (unsigned int i = 0; !different && (i < labels.size()); ++i){
            different = (labels[i]!=freqLabels[i]);
        }
    }

    
    if (onlyOneEquFreq && different){
        outputStream << "equilibrium frequencies:" << endl;
        outputFreqVector( labels, equFreqsMixture[0], outputStream, 12, 4 );
    }
    if (onlyOneExchMat && numberRatesRatios){
        outputStream << "rate ratios matrix:" << endl;
        outputMatrix( labels, exchMixture[0], outputStream, 9, 4 );
    }
    if (!onlyOneEquFreq && invariantCategory){
        outputStream << "frequency parameters, invariant sites:" << endl;
        outputFreqVector( labels, equFreqsMixture[0], outputStream, 12, 4 );
    }
    outputStream << "--------------------------" << endl;
    array2D<double> transitionMatrix;
    transitionMatrix.resize( nbStates, nbStates );
    if ( onlyOneEquFreq && onlyOneExchMat ){
        outputStream << "transition rate matrix (rescaled so that "
                     << "its average substitution rate is 1.0):" << endl;
        double sum = 0.0;
        for ( unsigned int i = 0; i < nbStates - 1; ++i ) {
            for ( unsigned int j = i + 1; j < nbStates; ++j ) {
                 sum += 2.0 * equFreqsMixture[0][i] *
                              equFreqsMixture[0][j] *
                              exchMixture[0](i,j);
            }
        }
        for ( unsigned int i = 0; i < nbStates; ++i ) {
            for ( unsigned int j = 0; j < nbStates; ++j ) {
                if(i!=j){
                    transitionMatrix(i,j) = equFreqsMixture[0][j] * exchMixture[0](i,j) / sum;
                }
            }
        }
        outputMatrix( labels, transitionMatrix, outputStream, 9, 4 );
        outputStream << "------------------------------" << endl;
    }
    else{
        for ( unsigned int gammaCat = 0; gammaCat < MAX(numberGammaCategories,1); ++ gammaCat ) {
            if(!onlyOneEquFreq){
                outputStream << "equilibrium frequencies:" << endl;
                outputFreqVector( labels, equFreqsMixture[gammaCat+invariantCategory], outputStream, 12, 4 );
            }
            if (!onlyOneExchMat ){
                outputStream << "rate ratios matrix:";
                outputMatrix( labels, exchMixture[gammaCat], outputStream, 9, 4 );
            }
            outputStream << "transition rate matrix (rescaled so that "
                         << "its average substitution rate is 1.0):" << endl;
            double sum = 0.0;
            for ( unsigned int i = 0; i < nbStates - 1; ++i ) {
                for ( unsigned int j = i + 1; j < nbStates; ++j ) {
                     sum += 2.0 * equFreqsMixture[gammaCat+invariantCategory][i] *
                                  equFreqsMixture[gammaCat+invariantCategory][j] *
                                  exchMixture[gammaCat](i,j);
                }
            }
            for ( unsigned int i = 0; i < nbStates; ++i ) {
                for ( unsigned int j = 0; j < nbStates; ++j ) {
                    if(i!=j){
                        transitionMatrix(i,j) = equFreqsMixture[gammaCat+invariantCategory][j] * exchMixture[gammaCat](i,j) / sum;
                    }
                }
            }
            outputMatrix( labels, transitionMatrix, outputStream, 9, 4 );
            outputStream << "------------------------------" << endl;
        }
    }
    outputStream << endl;
}
  
  
void MatrixModel::printPerturbationParameters( ostream & outputStream ) {
    assert( perturbator );
    perturbator->printPerturbationParameters( outputStream );
}


double MatrixModel::perturb() {
#ifdef DEBUG2  
    cout << "Matrix model " << this << " receives a perturb order, transmit to perturbator" << endl;
#endif
    return perturbator->perturb();
}


bool MatrixModel::validatePerturbation( bool validation ) {
#ifdef DEBUG2
    cout << "Matrix model " << this << " receives a validation " << validation << ", transmit to perturbator" << endl;
#endif
    bool res = perturbator->validatePerturbation( validation );
    return res;    
}

void MatrixModel::stopBurn(){
    perturbator->stopBurn();
}


void MatrixModel::getAllPerturbationParameters( vector<double>& params ) const{
    perturbator->getAllPerturbationParameters( params );
}

void MatrixModel::setAllPerturbationParameters( const vector<double>& params){
    perturbator->setAllPerturbationParameters( params );
}

unsigned int MatrixModel::getNumberPerturbationParameters() const{
    return perturbator->getNumberPerturbationParameters();
}


void MatrixModel::getAllPriorParameters( vector<double>& params ) const{
    perturbator->getAllPriorParameters( params );
    vector< double > p;
    frequencies->getAllPriorParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
    ratesRatios->getAllPriorParameters( p );
    params.insert( params.end(), p.begin(), p.end() );
    assert( params.size() == getNumberPriorParameters() );
}

void MatrixModel::setAllPriorParameters( const vector<double>& params){
    vector<double>::const_iterator iter = params.begin();

    vector<double> p;
    unsigned int nbp = perturbator->getNumberPriorParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    perturbator->setAllPriorParameters( p );
    nbp = frequencies->getNumberPriorParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    frequencies->setAllPriorParameters( p );
    nbp = ratesRatios->getNumberPriorParameters();
    p.resize(nbp);
    for ( unsigned int i = 0; i < nbp; ++i ){
        p[i] = *iter;
        ++iter;
    }
    ratesRatios->setAllPriorParameters( p );
    assert( iter == params.end() );
}

unsigned int MatrixModel::getNumberPriorParameters() const{
    return perturbator->getNumberPriorParameters() +
            frequencies->getNumberPriorParameters() +
            ratesRatios->getNumberPriorParameters();
}

double MatrixModel::getLnPrior() const{
#ifdef DEBUG2
    cout << "model prior on frequencies: " << frequencies->getLnPrior(averageRate) << endl;
    cout << "model prior on rates: " << ratesRatios->getLnPrior() << endl;
    cout << "model prior on parameters: " << perturbator->getLnPrior() << endl;
#endif    
    return (perturbator->getLnPrior()+frequencies->getLnPrior(averageRate)+ratesRatios->getLnPrior());
}



void MatrixModel::printLine( ostream & outputStream ) {
    outputStream.setf(ios::fixed);
    outputStream << setprecision( 4 );
    if ( invariantCategory ){
        outputStream << setw( 10 ) << proportionInvariantSites << ' ';
    }
    if (discreteGamma){
        outputStream << setw( 14 ) << alpha << ' ';
    }
    frequencies->printLine( outputStream );
    ratesRatios->printLine( outputStream );    
    outputStream.unsetf(ios::fixed);
}


void MatrixModel::printPriorLine( ostream & outputStream ) {
    vector< double > priors;
    getAllPriorParameters( priors );
    outputStream.setf(ios::fixed);
    outputStream << setprecision( 4 );
    for (unsigned int i = 0; i < priors.size(); ++i){
        outputStream << setw( 14 ) << priors[i] << ' ';
    }
    outputStream.unsetf(ios::fixed);
}


void MatrixModel::fromLine( const vector<double>& parameters ) {
    assert( parameters.size() == getNumberLineParameters() );

    vector<double>::const_iterator iter1 = parameters.begin();
    if ( invariantCategory ){
        proportionInvariantSites = *iter1;
        ++iter1;
    }
    if (discreteGamma){
        alpha = *iter1;
        ++iter1;
    }
    vector<double>::const_iterator iter2 = iter1 + frequencies->getNumberLineParameters();
    frequencies->fromLine( vector<double>( iter1, iter2 ) );
    iter1 = iter2 + ratesRatios->getNumberLineParameters();
    ratesRatios->fromLine( vector<double>( iter2, iter1 ) );    
    
    updateAverageRateVector();
    updateEigenMatrix();
}

unsigned int MatrixModel::getNumberLineParameters() const{
    return frequencies->getNumberLineParameters() +
           ratesRatios->getNumberLineParameters() +
           invariantCategory + (discreteGamma ? 1 : 0 );
}

void MatrixModel::errorSymbol( SequenceTable * sequenceTable,
                               unsigned int modelId, unsigned int speciesId,
                               int site){
    vector < vector< unsigned int > > pos;
    string symbol;
    //invariant ??
    if (site < 0){
        symbol = sequenceTable->getInvariantBases( modelId )[-site-1].first;
    }
    //usual
    else{
        symbol = sequenceTable->getSequences( modelId )( speciesId, site );
        cerr << endl << "Species: " << sequenceTable->species[speciesId] << ", ";
    }
    
    cerr << "Unrecognized symbol: " << symbol << ", pos: ";
    pos = sequenceTable->retrieveInitialNucleotides( modelId, site );        
    for ( vector < vector< unsigned int > >::iterator iter = pos.begin();
         iter != pos.end(); ++iter ){
        for ( vector< unsigned int >::iterator iter2 = (*iter).begin();
              iter2 != (*iter).end(); ++iter2 ){
            if (iter2 != (*iter).begin() ){
                cout << ", ";
            }
            cout << *iter2 + 1;
        }
        cout << "; ";
    }
    cout << endl;
    exit(EXIT_FAILURE);
}

void MatrixModel::validChange(){
    updateAverageRateVector();
    updateEigenMatrix();
}

void MatrixModel::updateAverageRateVector(){
    if ( discreteGamma ) {
        double beta = alpha * ( 1.0 - proportionInvariantSites )/averageSubstitutionRate;
        if (percentileGammaApproximation){
            statlib::gamma_percentile_means( averageRate, alpha, beta,
            numberGammaCategories );
        }
        else{
            statlib::gammaLaguerre( averageRate, categoriesWeights, alpha, beta,
                               numberGammaCategories );
        }
        if (invariantCategory){
            vector< double > rates;
            rates.push_back( 0 );
            rates.insert( rates.end(), averageRate.begin(), averageRate.end() );
            frequencies->prepare( rates );
        }
        else{
            frequencies->prepare( averageRate );
        }
        ratesRatios->prepare( averageRate );
    }
    else{
        averageRate[0] = averageSubstitutionRate/( 1.0 - proportionInvariantSites );
    }
}

void MatrixModel::setAverageSubstitutionRate( double newAverageRate, unsigned int ){
    assert( newAverageRate > 0.0 );
    //if ( newAverageRate < 1e-8 ) {
	//	newAverageRate = 1e-8;
	//}
    averageSubstitutionRate = newAverageRate;
}

void MatrixModel::updateEigenMatrix(){

    //BEWARE, function redefined in many descendants
    //TwoStateModel, GG, JC69, K80, RNA7A,RNA7D
    
    double sum = 0.0;

    if ( discreteGamma ) {
        for ( unsigned int gammaCategory = 0 ;
              gammaCategory < numberGammaCategories; ++gammaCategory ) {
            int freqCat = gammaCategory + invariantCategory;
            sum = 0.0;
            for ( unsigned int i = 0; i < getNumberStates() -1; ++i ) {
                for ( unsigned int j = i + 1; j < getNumberStates(); ++j ) {
                    sum += 2.0 * getFrequency( i, freqCat ) *
                                 getFrequency( j, freqCat ) *
                                 getExchangeability( i, j, gammaCategory );
                }
            }
            substitutionRate[gammaCategory] = averageRate[gammaCategory]/sum;
        }
    }
    else{
        sum = 0.0;
        int freqCat = invariantCategory;
        for ( unsigned int i = 0; i < getNumberStates() - 1; ++i ) {
            for ( unsigned int j = i + 1; j < getNumberStates(); ++j ) {
                sum += 2.0 * getFrequency( i, freqCat ) *
                             getFrequency( j, freqCat ) *
                             getExchangeability( i, j, 0 );
            }
        }
        substitutionRate[0] = averageRate[0] / sum;
    }
    setEigenMatrix();
}

string MatrixModel::getName() const{
    char num[3];
    string name;
    if ( discreteGamma ) {
        sprintf( num, "%i", numberGammaCategories );
        name = name + string(" + dG") + string( num );
    }
    if ( invariantCategory == 1 ){
        name = name + string( " + I" );
    }
    //if linked s or linked x
    if ( frequencies->getGammaCatAdjust() & Frequencies::LINK ){
        if ( frequencies->getGammaCatAdjust() & Frequencies::LS ){
            name = name + string(" + s");
        }
        else{
            sprintf( num, "%i", frequencies->getNumberFrequenciesXSets() );
            name = name + string(" + x") + string(num);
        }
        return name;
    }
    //in other case
    if (frequencies->getNumberFrequenciesSets() > 1){
        if ( frequencies->getGammaCatAdjust() & Frequencies::L ){
            if ( frequencies->getGammaCatAdjust() & Frequencies::LX ){
                sprintf( num, "%i", frequencies->getNumberFrequenciesXSets() );
                name = name + string(" + xf") + string(num);
            }
            else{
                if ( frequencies->getGammaCatAdjust() & Frequencies::LS ){
                    name = name + string(" + sf");
                }
                else{
                    name = name + string(" + lf");
                }
            }
        }
        else{
            sprintf( num, "%i", frequencies->getNumberFrequenciesSets() );
            name = name + string(" + f") + string(num);
        }
    }
    if (ratesRatios->getNumberRatesRatiosSets() > 1){
        if ( ratesRatios->getGammaCatAdjust() & RatesRatios::L ){
            if ( ratesRatios->getGammaCatAdjust() & RatesRatios::LX ){
                sprintf( num, "%i", ratesRatios->getNumberRatesRatiosXSets() );
                name = name + string(" + xr") + string(num);
            }
            else{
                if ( ratesRatios->getGammaCatAdjust() & RatesRatios::LS ){
                    name = name + string(" + sr");
                }
                else{
                    name = name + string(" + lr");
                }
            }
        }
        else{
            sprintf( num, "%i", ratesRatios->getNumberRatesRatiosSets() );
            name = name + string(" + r") + string(num);
        }
    }
    return name;
}

void MatrixModel::update(UpdateMessage* subject){
    clear();
    
    setType( UpdateMessage::MODEL_TYPE );
#ifdef DEBUG3
    cout << "MatrixModel:received an update message ";
#endif
    if (!subject->hasFlag(PRIOR_FLAG)){
#ifdef DEBUG3
        cout << " without PRIOR flag" << endl;
#endif
        if ( subject->hasFlag(GAMMA) ||
             subject->hasFlag(INVARIANT) ||
             subject->hasFlag(AVERAGE_RATE) ){
            updateAverageRateVector();
        }
        updateEigenMatrix();
    }
    else{
#ifdef DEBUG3
        cout << " with PRIOR flag" << endl;
#endif
        setFlag(MODEL_PRIOR_FLAG);
    }
#ifdef DEBUG3
    cout << "processed, now transmit" << endl;
#endif
    notify();
}


