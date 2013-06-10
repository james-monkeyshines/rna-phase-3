#include "Models/Frequencies.h"

#include "Util/ParametersSet.h"
#include "Util/matrixmath.h"
#include "Models/MatrixModel.h"
#include "Models/Perturbator.h"
#include "Models/PerturbatorGaussParameter.h"
#include "Models/PerturbatorFrequenciesParameter.h"

#include <iostream>
#include <iomanip>
#include <float.h>
#include <numeric>

#define SQR(n) (n*n)

Frequencies::Frequencies( ParametersSet& parameters, MatrixModel* model ){
    initPrimitive( parameters, model->getNumberGammaCategories(),
            model->getInvariant() );
}

void Frequencies::initPrimitive(ParametersSet& parameters, unsigned int numberGammaCategories, unsigned int invariant ){

    this->invariant = invariant;
    this->numberGammaCategories = numberGammaCategories;

    //initialisation: no adjustment of frequencies by default
    gammaCatAdjust = 0;
    variationType = UNIFORM;
    numberFrequenciesXSets = 0;
    sF = NULL;
    xF = NULL;
    frequenciesX = NULL;

    //default value for the number of frequencies categories equals to 1
    //model without frequencies (numberFrequencies = 0)
    //must make sure this value is filled with 0 before construction
    if ( parameters.findParameter( "Number of frequencies sets" ) ) {
        numberFrequenciesSets = parameters.intParameter( "Number of frequencies sets" );
    }
    else{
        numberFrequenciesSets = 1;
    }

    string frequenciesAdjustment;
    if ( parameters.findParameter( "Frequencies adjustment" ) ){
        if (numberFrequenciesSets==0){
            cerr << "No \"Frequencies adjustment\" possible for a model without free frequencies parameters" << endl;
            exit(EXIT_FAILURE);
        }
        if (numberFrequenciesSets>1){
            cerr << "Do not use \"Number of frequencies sets\" and \"Frequencies adjustment\" simultaneously" << endl;
            exit(EXIT_FAILURE);
        }
        frequenciesAdjustment = parameters.stringParameter("Frequencies adjustment");
        if ( frequenciesAdjustment == "linear" ){
            gammaCatAdjust = L;
        }
        else if ( frequenciesAdjustment=="linearS" ){
            gammaCatAdjust = LS | L;
            if ( ( !parameters.findParameter( "Link S parameter" ) ) ||
                 ( !parameters.boolParameter("Link S parameter") ) ){
                sF = new double;
                *sF = 0.5;
            }
            else{
                gammaCatAdjust = gammaCatAdjust | LINK;
            }
        }
        else if ( frequenciesAdjustment=="linearX" ){
            gammaCatAdjust = LX | L;
            numberFrequenciesXSets = parameters.intParameter("Number of X frequencies sets");
            if (numberFrequenciesXSets < 1){
                cerr << "Invalid number of X frequencies sets" << endl;
                exit(EXIT_FAILURE);
            }
            frequenciesX = new vector< double >[numberFrequenciesXSets];
            if  ( !parameters.findParameter("Link X parameter") ||
                  !parameters.boolParameter("Link X parameter") ){
                xF = new vector<double>;
                xF->resize( numberFrequenciesXSets );
                if ( gammaCatAdjust & LX ){
                    for ( unsigned int i = 0; i < numberFrequenciesXSets; ++i ){
                        (*xF)[i] = .5;
                    }
                }
            }
            else{
                gammaCatAdjust = gammaCatAdjust | LINK;
            }
        }
        else{
            cerr << "Unrecognized value \"" << frequenciesAdjustment << "\" for the parameters "
                 << "\"Frequencies adjustment\"... abort" << endl;
            exit(EXIT_FAILURE);
        }
    }

    //consistency checking at that point, in case of identifiability issue
    unsigned int numberRatesCategories;

    if (numberGammaCategories){
        numberRatesCategories = numberGammaCategories + invariant;
        if ( numberFrequenciesXSets > numberRatesCategories ){
            cerr << "Invalid number of frequencies sets (" << numberFrequenciesXSets
                 << "). The expected number is in [1.." << numberRatesCategories << "]." << endl;
            exit(EXIT_FAILURE);
        }
        // if no invariant and 2 gamma category and LSF, sF not identifiable
        if ( ( numberRatesCategories == 2 ) && (gammaCatAdjust & LS) ){
            cerr << "Identifiability issue, with 2 set of frequencies (for the two gamma categories),"
                 << " you cannot evaluate sF" << endl;
            exit(EXIT_FAILURE);
        }
        if ( gammaCatAdjust & LX ){
            if ( numberFrequenciesXSets > numberRatesCategories ){
                cerr << "Identifiability issue, \"number of X rates categories\" > \"Number of rates categories\"" << endl;
                exit(EXIT_FAILURE);
            }
            if ( numberFrequenciesXSets > (numberRatesCategories - 2)/2 ){
                cerr << "WARNING: You probably should not be using more than "
                     << (numberRatesCategories - 2)/2 << " X frequencies categories." << endl;
            }
        }
    }
    else{
        numberRatesCategories = 1 + invariant;
        //no linear frequencies without +I
        if (gammaCatAdjust & L){
            if ( !invariant ){
                 cerr << "Cannot have a \"Frequencies adjustment\"  with a model without gamma" << endl
                      << "distribution or invariant. Exit..." << endl;
                 exit(EXIT_FAILURE);
            }
            else{
                // no linearS/X frequencies allowed without gamma model.
                if ( (gammaCatAdjust & LS) || (gammaCatAdjust & LX) ){
                    cerr << "Identifiability issue, a model +I without rate variation across sites cannot" << endl
                         << "be used with \"Frequencies adjustment = " << frequenciesAdjustment << '\"' << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    if ( numberFrequenciesSets > numberRatesCategories ){
        cerr << "Error, not enough categories in your mixture model (" << numberRatesCategories
             << ") to accomodate " << numberFrequenciesSets << " sets of frequency parameters" << endl;
        exit(EXIT_FAILURE);             
    }
/** Set up frequencies without the discrete gamma model **************************************************/
    if ( !numberGammaCategories ) {
        // if the model uses some frequencies
        if ( numberFrequenciesSets != 0 ){
            frequenciesCategory.resize( numberRatesCategories );
            if ( !invariant ){
                assert(numberFrequenciesSets == 1); //check should be already done
                frequenciesCategory[0] = 0;
            }
            // two allowed categories otherwise (but not compulsory, it is user choice)
            else{
                assert( (numberFrequenciesSets == 1) ||(numberFrequenciesSets == 2) ); //check should be already done
                //if a linear frequency adjustment is used act as if the user asked for two categories
                if ( gammaCatAdjust & L ) {
                    numberFrequenciesSets = 2;
                }
                frequenciesCategory[0] = 0;
                frequenciesCategory[1] = numberFrequenciesSets-1;
            }
            frequencies = new vector<double>[numberFrequenciesSets];
        }
        else{  //no frequencies in the model
            frequenciesCategory.resize(0);
        }
    }
/** Set up frequencies with the discrete gamma model **************************************************/
    else {
        // if the model uses some frequencies
        if ( numberFrequenciesSets != 0 ){
            frequenciesCategory.resize( numberRatesCategories );
            if ( gammaCatAdjust & L ){
                numberFrequenciesSets = numberRatesCategories;
            }
            //if no indication provided
            //initialise the frequencyCategory vector used to match a rate category
            //with a matrix
            if ( !parameters.findParameter( "Frequencies for category 1" ) ) {
                // allow only two simple case only one set of frequencies for all the
                // rate categories(number of gamma categories + invariant category )
                // or one set for each category
                if ( ( numberFrequenciesSets != 1 ) &&
                     ( numberFrequenciesSets != numberRatesCategories ) ){
                    cerr << "Since \"Number of frequencies sets\" is not equal to " << numberRatesCategories
                         << ". You have to specify the set to be used with each category." << endl;
                    if (!invariant){
                        cerr << "Use the fields \"Frequencies for category 1\", \"Frequencies for category 2\", ..." << endl;
                    }
                    else{
                        cerr << "Use the fields \"Frequencies for invariant category\", \"Frequencies for category 1\", "
                             << "\"Frequencies for category 2\", ..." << endl;
                    }
                    exit(EXIT_FAILURE);
                }
                if ( numberFrequenciesSets == 1 ) {
                    for ( unsigned int i = 0; i < numberRatesCategories; ++i ) {
                        this->frequenciesCategory[i] = 0;
                    }
                }
                if ( numberFrequenciesSets == numberRatesCategories ) {
                    for (unsigned  int i = 0; i < numberRatesCategories; ++i ) {
                        this->frequenciesCategory[i] = i;
                    }
                }
            }
            else {
                if (gammaCatAdjust & L){
                    cerr << "Do not specify the frequencies set for each gamma category "
                         << "when using the field \"Frequencies adjustment\"" << endl;
                    exit(EXIT_FAILURE);
                }
                if (invariant){
                    frequenciesCategory[0] = parameters.intParameter( "Frequencies for invariant category" ) - 1;
                    if( (frequenciesCategory[0]<0) && ((unsigned int)frequenciesCategory[0]>=numberFrequenciesSets) ){
                        cerr << "ERROR ! Cannot assign rate matrice " << frequenciesCategory[0] + 1
                             << "for the invariant rate category " << endl;
                        exit(EXIT_FAILURE);
                    }
                    if ( frequenciesCategory[0] != 0 ){
                        cerr << "WARNING ! Unusual allocation for the frequencies of"
                             << "the invariant category, 1 is usually used" << endl;
                    }
                }
                int cat = 0;
	            char label[35];
                sprintf( label, "Frequencies for category ");
                char* numberPos = label;
                while (*(++numberPos));
                for ( unsigned int i = 0; i < numberGammaCategories; ++i ) {
                    sprintf( numberPos, "%d", i + 1 );
                    frequenciesCategory[i+invariant] = parameters.intParameter( label )-1;
                    if ( frequenciesCategory[i+invariant] < cat ) {
                        cerr << "WARNING ! Unusual allocation of frequencies for "
                             << "the rate category " << i+1 << endl;
                    }
                    cat = frequenciesCategory[i+invariant];
                    if ( ( cat < 0 ) || ( (unsigned int)cat >= numberFrequenciesSets ) ) {
                        cerr << "ERROR ! assigning frequencies set " << cat+1
                             << "for the rate category " << i + 1 << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
            frequencies = new vector<double>[numberFrequenciesSets];
        }
        else{  //no frequencies in the model
            frequenciesCategory.resize(0);
        }
    } // end with discreteGamma
    lnPrior = 0.0;
}


void Frequencies::initialisation( unsigned int numberFrequencies ){
    vector< double > freq;
    freq.resize(numberFrequencies);
    for ( unsigned int i = 0; i <  numberFrequencies; ++i ) {
        freq[i] = 1.0 / numberFrequencies;
    }
    initialisation( freq );
}

void Frequencies::initialisation( const vector<double>& freq ){
    this->numberFrequencies = freq.size();
    assert( !freq.size()||
            (fabs(accumulate(freq.begin(),freq.end(),0.0)-1.0) < .001));
    if ( numberFrequenciesSets ){
        for ( unsigned int i = 0; i < numberFrequenciesSets; ++i ) {
            frequencies[i] = freq;
        }
        if (gammaCatAdjust & LX){
            for ( unsigned int i = 0; i < numberFrequenciesXSets; ++i ) {
                frequenciesX[i] = freq;
            }
        }
    }
    else{
       // frequencies = NULL;
       // frequencies = new vector<double>;
       // *frequencies = freq;
    }
}


void Frequencies::updateLinear(){
    if (gammaCatAdjust & LX){
        unsigned int x = 0;
        double highRate = savedRates[0];
        double lowRate = highRate;
        vector<double>* highFreq = &(frequencies[0]);
        vector<double>* lowFreq = highFreq;
        double ratio;
        for ( unsigned int i = 1; i < numberFrequenciesSets-1; ++i ){
            while( savedRates[i] >= highRate ){
                lowRate = highRate;
                lowFreq = highFreq;
                if ( x == numberFrequenciesXSets ){
                    highRate = savedRates.back();
                    highFreq = &(frequencies[numberFrequenciesSets-1]);
                }
                else{
                    highRate = (*xF)[x];
                    for (unsigned int j = x+1; j < numberFrequenciesXSets; ++j){
                        highRate = highRate * (*xF)[j];
                    }
                    highRate = highRate * (savedRates.back()-savedRates.front()) + savedRates.front();
                    highFreq = &(frequenciesX[x]);
                    ++x;
                }
            }
            ratio = (savedRates[i]-lowRate)/(highRate - lowRate);
            for ( unsigned int j = 0; j < numberFrequencies; ++j ){
                frequencies[i][j] = (*lowFreq)[j] + ( (*highFreq)[j] - (*lowFreq)[j] ) * ratio;
            }
        }
    }
    //else linearS
    else{
        assert(gammaCatAdjust & L);
        //fill the frequencies between the lowest and the highest category
        double delta = savedRates.back() - savedRates.front();
        double ratio;
        for ( unsigned int i = 1; i < numberFrequenciesSets-1; ++i ){
            ratio = (savedRates[i]-savedRates.front())/delta;
            //if linearS
            if (gammaCatAdjust & LS){
                ratio = pow( ratio, *sF );
            }
            for ( unsigned int j = 0; j < numberFrequencies; ++j ){
                frequencies[i][j] = frequencies[0][j]+(frequencies[numberFrequenciesSets-1][j]-frequencies[0][j])*ratio;
            }
        }
    }
}

void Frequencies::prepare( const vector<double> & rates ){
    if ( (gammaCatAdjust & L) || (variationType == GP) ){
        savedRates = rates;
        if (gammaCatAdjust & L){
            updateLinear();
        }
        if (variationType == GP){
            computeCovMatrix();
            lnPrior = computeLnPrior();
        }
    }
}

double* Frequencies::linkS(){
    assert( gammaCatAdjust & LINK );
    assert( !(gammaCatAdjust & LINKN) );
    assert( !(gammaCatAdjust & LINKY) );
    assert( gammaCatAdjust & LS );
    gammaCatAdjust = gammaCatAdjust | LINKY;

    sF = new double;
    *sF = .5;
    return sF;
}

void Frequencies::assignS( double* sF ){
    assert( gammaCatAdjust & LINK );
    assert( !(gammaCatAdjust & LINKN) );
    assert( !(gammaCatAdjust & LINKY) );
    assert( gammaCatAdjust & LS );

    gammaCatAdjust = gammaCatAdjust | LINKN;
    this->sF = sF;
}

vector<double>* Frequencies::linkX(){
    assert( gammaCatAdjust & LINK );
    assert( !(gammaCatAdjust & LINKN) );
    assert( !(gammaCatAdjust & LINKY) );
    assert( gammaCatAdjust & LX );

    gammaCatAdjust = gammaCatAdjust | LINKY;
    xF = new vector<double>;
    xF->resize( numberFrequenciesXSets );
    for ( unsigned int i = 0; i < numberFrequenciesXSets; ++i ){
        (*xF)[i] = .5;
    }
    return xF;
}

void Frequencies::assignX( vector<double>* xF ){
    assert( gammaCatAdjust & LINK );
    assert( !(gammaCatAdjust & LINKN) );
    assert( !(gammaCatAdjust & LINKY) );
    assert( gammaCatAdjust & LX );

    if( xF->size() != numberFrequenciesXSets ){
        cerr << "Cannot link two X parameters, sizes are different "
             << numberFrequenciesXSets << "<>" << xF->size() << endl;
        exit(EXIT_FAILURE);
    }

    gammaCatAdjust = gammaCatAdjust | LINKN;
    this->xF = xF;
}

void Frequencies::getModelParameters( ParametersSet & parameters, const Model* model, unsigned int symbolCategory ) const{
    ParametersSet* currentCategory;
    char tmpString[40];
    char dblString[40];

    if (numberFrequenciesSets>=1){
        if (numberFrequenciesSets>=2){
            sprintf(tmpString,"%d",numberFrequenciesSets);
            parameters("FREQUENCIES")["Number of frequencies sets"] = tmpString;
            if (gammaCatAdjust & LS){
                sprintf(dblString, "%.8f", *sF);
                parameters("FREQUENCIES")["Variation parameter s"] = dblString;
            }
            unsigned int increment = (gammaCatAdjust&L) ? numberFrequenciesSets - 1 : 1;
            for ( unsigned int i = 0; i < numberFrequenciesSets; i += increment ){
                sprintf(tmpString,"FREQUENCIESCATEGORY%d",i+1);
                currentCategory = &parameters("FREQUENCIES")(tmpString);
                for (unsigned int j = 0; j < numberFrequencies; ++j){
                    sprintf(dblString, "%.8f",frequencies[i][j]);
                    (*currentCategory)["F("+model->getFrequencyState(j,symbolCategory)+")"] = dblString;
                }
            }
            if (gammaCatAdjust & LX){
                double ratio = 1.0;
                for ( int i = numberFrequenciesXSets-1; i >= 0; --i ){
                    sprintf(tmpString,"FREQUENCIESCATEGORYX%d",i+1);
                    currentCategory = &parameters("FREQUENCIES")(tmpString);
                    ratio = ratio * (*xF)[i];
                    sprintf(dblString, "%.8f", ratio);
                    (*currentCategory)["Variation parameter x"] = dblString;
                    for (unsigned int j = 0; j < numberFrequencies; ++j){
                        sprintf(dblString, "%.8f",frequenciesX[i][j]);
                        (*currentCategory)["F("+model->getFrequencyState(j,symbolCategory)+")"] = dblString;
                    }
                }
            }
            //store the matching between rate category and set of frequencies
            //in user defined matching case
            if (!gammaCatAdjust&L){
                for ( unsigned int i = 0; i < numberFrequenciesSets; ++i ){
                    sprintf(tmpString,"RATESCATEGORY%d",i+1);
                    currentCategory = &parameters("RATESCATEGORIES")(tmpString);
                    sprintf(dblString,"%d",frequenciesCategory[i]+1);
                   (*currentCategory)["Frequencies set for the category"] = dblString;
                }
            }
        }
        else{
            for (unsigned int j = 0; j < numberFrequencies; ++j){
                sprintf(dblString, "%.8f",frequencies[0][j]);
                parameters("FREQUENCIES")["F("+model->getFrequencyState(j,symbolCategory)+")"] = dblString;
            }
        }
    }

}

void Frequencies::setModelParameters( ParametersSet & parameters, Model* model, unsigned int symbolCategory ){
    ParametersSet* currentCategory;
    char tmpString[40];

    if ( numberFrequenciesSets >= 1 ){
        if ( numberFrequenciesSets >= 2 ){
            if( numberFrequenciesSets !=
                   (unsigned int)parameters("FREQUENCIES").intParameter("Number of frequencies sets") ){
                cerr << "The number of frequencies sets does not match in the model file." << endl;
                exit(EXIT_FAILURE);
            }
            if (gammaCatAdjust & LS){
                *sF = parameters("FREQUENCIES").doubleParameter("Variation parameter s");
            }
            unsigned int increment = (gammaCatAdjust&L) ? numberFrequenciesSets - 1 : 1;
            for ( unsigned int i = 0; i < numberFrequenciesSets; i += increment ){
                sprintf(tmpString,"FREQUENCIESCATEGORY%d",i+1);
                currentCategory = &parameters("FREQUENCIES")(tmpString);
                double sum = 0.0;
                for (unsigned int j = 0; j < numberFrequencies; ++j){
                    frequencies[i][j] = currentCategory->doubleParameter("F("+model->getFrequencyState(j,symbolCategory)+")");
                    sum += frequencies[i][j];
                }
                if(fabs(sum-1.0)>=0.02){
                    cerr << "WARNING: frequencies do not sum to one in your model file ("
                         << sum << ')' << endl;
                }
                for (unsigned int j = 0; j < numberFrequencies; ++j){
                    frequencies[i][j]/=sum;
                }
            }
            if (gammaCatAdjust & LX){
                for ( unsigned int i = 0; i < numberFrequenciesXSets; ++i ){
                    sprintf(tmpString,"FREQUENCIESCATEGORYX%d",i+1);
                    currentCategory = &parameters("FREQUENCIES")(tmpString);
                    double sum = 0.0;
                    (*xF)[i] = currentCategory->doubleParameter("Variation parameter x");
                    for (unsigned int j = 0; j < numberFrequencies; ++j){
                        frequenciesX[i][j] = currentCategory->doubleParameter("F("+model->getFrequencyState(j,symbolCategory)+")");
                        sum += frequenciesX[i][j];
                    }
                    if(fabs(sum-1.0)>=0.02){
                        cerr << "WARNING: frequencies do not sum to one in your model file ("
                         << sum << ')' << endl;
                    }
                    for (unsigned int j = 0; j < numberFrequencies; ++j){
                        frequenciesX[i][j]/=sum;
                    }
                }
                for ( unsigned int i = 0; i < numberFrequenciesXSets-1; ++i ){
                    (*xF)[i] = (*xF)[i] / (*xF)[i+1];
                }
            }
            if (!gammaCatAdjust&L){
                for ( unsigned int i = 0; i < numberFrequenciesSets; ++i ){
                    sprintf(tmpString,"RATESCATEGORY%d",i+1);
                    currentCategory = &parameters("RATESCATEGORIES")(tmpString);
                    frequenciesCategory[i] = currentCategory->intParameter("Frequencies set for the category") - 1;
                    assert( (frequenciesCategory[i]>=0) &&
                                  ((unsigned int)frequenciesCategory[i]<numberFrequenciesSets) );
                }
            }
       }
        else{
            double sum = 0.0;
            for (unsigned int j = 0; j < numberFrequencies; ++j){
                frequencies[0][j] = parameters("FREQUENCIES").doubleParameter("F("+model->getFrequencyState(j,symbolCategory)+")");
                sum += frequencies[0][j];
            }
            if(fabs(sum-1.0)>=0.02){
                cerr << "WARNING: frequencies do not sum to one in your model file ("
                     << sum << ')' << endl;
            }
            for (unsigned int j = 0; j < numberFrequencies; ++j){
                frequencies[0][j]/=sum;
            }
        }
    }
}

void Frequencies::printParameters( ostream & outputStream ) const{
    if( (!(gammaCatAdjust & LS) && !(gammaCatAdjust & LX)) || (gammaCatAdjust & LINKN) ){
        return;
    }
    if (!(gammaCatAdjust & LINK)){
        outputStream << "Frequencies, ";
    }
    //if S
    if (gammaCatAdjust & LS){
        outputStream << "S variation parameter = " << *sF << endl;
    }
    //if X
    if (gammaCatAdjust & LX){
        outputStream << "X variation parameters = ";
        vector< double > x = *xF;
        for (unsigned int i = numberFrequenciesXSets-1; i > 0 ; --i){
            x[i-1]=x[i]*x[i-1];
        }
        for (unsigned int i = 0; i < numberFrequenciesXSets; ++i){
            if (i) {
                outputStream << ", ";
            }
            outputStream << x[i];
        }
        outputStream << endl;
    }
}


void Frequencies::registerVariationParameter( ParametersSet & parameters, Perturbator* perturbator ){
    char label[50];
    char* numberPos;

    double extrMin, extrMax;

    if ( !(gammaCatAdjust & LINK) ){
        if (gammaCatAdjust & LS){
            if (!parameters.findParameter("Frequencies variation parameter S, prior")){
                parameters["Frequencies variation parameter S, prior"] = "uniform(0,1.0)";
            }
            PerturbatorGaussParameter* pert =
                perturbator->registerGaussParameter( *this, UpdateMessage::FREQ,
                           "Frequencies variation parameter S", sF,
                           parameters, "Frequencies variation parameter S",
                           .05, .001, 0.5, "initial step" );
            pert->getExtrema( extrMin, extrMax );
            if (extrMin<0.0){
                cerr << "invalid prior for the frequencies variation parameter S: "
                     << parameters.stringParameter("Frequencies variation parameter S, prior") << endl;
                exit(EXIT_FAILURE);
            }
        }
        if (gammaCatAdjust & LX){
            sprintf( label, "Frequencies variation parameter X");
            string genericLabel = string(label);
            if (!parameters.findParameter(genericLabel+", prior")){
                parameters[genericLabel+", prior"] = "uniform(0,1.0)";
            }
            numberPos = label;
            while (*(++numberPos));
            for (unsigned int i = 0; i < numberFrequenciesXSets; ++i){
                sprintf( numberPos, "%d, ", i+1 );
                PerturbatorGaussParameter* pert =
                    perturbator->registerGaussParameter( *this, UpdateMessage::FREQ, label, &((*xF)[i]),
                               parameters, genericLabel,
                               .05, .001, 0.5, "initial step" );
                pert->getExtrema( extrMin, extrMax );
                if ( (extrMin<0.0) || (extrMax>1.0) ){
                    cerr << "invalid prior for the " << label << ": "
                         << parameters[string(label)+", prior"]
                         << ". Should be bounded between 0.0 and 1.0" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    else{
        if (gammaCatAdjust & LINKY){
            if (gammaCatAdjust & LS){
                if (!parameters.findParameter("Variation parameter S, prior")){
                    parameters["Variation parameter S, prior"] = "uniform(0,1.0)";
                }
                PerturbatorGaussParameter* pert =
                    perturbator->registerGaussParameter( *this, UpdateMessage::FREQ,
                               "Variation parameter S", sF,
                               parameters, "Variation parameter S",
                               .05, .001, 0.5, "initial step" );
                pert->getExtrema( extrMin, extrMax );
                if (extrMin<0.0){
                    cerr << "invalid prior for the variation parameter S: "
                         << parameters["Variation parameter S, prior"] << endl;
                    exit(EXIT_FAILURE);
                }
            }
            if (gammaCatAdjust & LX){
                sprintf ( label, "Variation parameter X");
                numberPos = label;
                while (*(++numberPos));
                for (unsigned int i = 0; i < numberFrequenciesXSets; ++i){
                    sprintf( numberPos, "%d", i+1 );
                    if (!parameters.findParameter(string(label)+", prior")){
                        parameters[string(label)+", prior"] = "uniform(0,1.0)";
                    }
                    PerturbatorGaussParameter* pert =
                        perturbator->registerGaussParameter( *this, UpdateMessage::FREQ,
                                   label, &((*xF)[i]), parameters, label,
                                   .05, .001, 0.5, "initial step");
                    pert->getExtrema( extrMin, extrMax );
                    if ( (extrMin<0.0) || (extrMax>1.0) ){
                        cerr << "invalid prior for the " << label << ": "
                             << parameters[string(label)+", prior"]
                             << ". Should be bounded between 0.0 and 1.0" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    }
}



void Frequencies::getAllParameters( vector < double > & pars ) const {
    double sum;
    unsigned int increment;
    int index = -1;

    pars.resize( getNumberFreeParameters() );
    if ( (gammaCatAdjust & LS) && (!(gammaCatAdjust & LINK)||(gammaCatAdjust & LINKY)) ){
        pars[++index] = sqrt( *sF / ( 1.0 - *sF ) );
    }
    if ( (gammaCatAdjust & LX) && (!(gammaCatAdjust & LINK)||(gammaCatAdjust & LINKY)) ){
       //the variables optimized are bounded in [0..1]
        for (unsigned int i = 0; i < numberFrequenciesXSets; ++i){
            pars[++index] = sqrt( (*xF)[i] / ( 1.0 - (*xF)[i] ) );
        }
    }
    //we use a trick to limit the code when linear/linearS/linearX frequencies are used
    increment = (gammaCatAdjust & L) ? numberFrequenciesSets - 1 : 1;
    for ( unsigned int freqCat = 0; freqCat < numberFrequenciesSets; freqCat += increment) {
        sum = ( 1.0 - frequencies[freqCat][numberFrequencies - 1] ) /
        frequencies[freqCat][numberFrequencies - 1];
#ifdef DEBUG3
        assert( !isnan( sum ) );
        assert( !isinf( sum ) );
#endif
        for ( unsigned int freq = 0; freq < numberFrequencies - 1; ++freq ) {
            pars[++index] = sqrt( ( 1.0 + sum ) * frequencies[freqCat] [freq] );
        }
    }
    if (gammaCatAdjust & LX){
        for ( unsigned int freqCat = 0; freqCat < numberFrequenciesXSets; ++freqCat ) {
            sum = ( 1.0 - frequenciesX[freqCat] [numberFrequencies - 1] ) /
                                frequenciesX[freqCat] [numberFrequencies - 1];
#ifdef DEBUG3
            assert( !isnan( sum ) );
            assert( !isinf( sum ) );
#endif
            for ( unsigned int freq = 0; freq < numberFrequencies - 1; ++freq ) {
                pars[++index] = sqrt( ( 1.0 + sum ) * frequenciesX[freqCat] [freq] );
            }
        }
    }
    assert(index+1==(int)pars.size());
}

void Frequencies::setAllParameters(const vector < double > & pars ){
    double sum;
    int index = 0;

    assert(pars.size()==getNumberFreeParameters());
    if ( (gammaCatAdjust & LS) && (!(gammaCatAdjust & LINK)||(gammaCatAdjust & LINKY)) ){
        *sF = SQR(pars[index]) / ( 1.0 + SQR( pars[index] ) );
        ++index;
    }
    if ( (gammaCatAdjust & LX) && (!(gammaCatAdjust & LINK)||(gammaCatAdjust & LINKY)) ){
       //the variables optimized are bounded in [0..1]
        for (unsigned int i = 0; i < numberFrequenciesXSets; ++i){
            (*xF)[i] = SQR( pars[index] ) / ( 1.0 + SQR(pars[index]) ) ;
            ++index;
        }
    }
    //we use a trick to limit the code when linear frequencies are used
    unsigned int increment = (gammaCatAdjust & L) ? numberFrequenciesSets - 1 : 1;
    for ( unsigned int freqCat = 0; freqCat < numberFrequenciesSets; freqCat += increment) {
        sum = 1.0;
        for ( unsigned int i = 0; i < numberFrequencies - 1; ++i ) {
            frequencies[freqCat][i] = SQR(pars[index]);
            ++index;
            sum += frequencies[freqCat][i];
        }
        for ( unsigned int i = 0; i < numberFrequencies - 1; ++i ) {
            frequencies[freqCat] [i] /= sum;
#ifdef DEBUG3
            assert( !isnan( frequencies[freqCat] [i] ) );
            assert( !isinf( frequencies[freqCat] [i] ) );
#endif
        }
        frequencies[freqCat] [numberFrequencies - 1] = 1.0 / sum;
    }
    if (gammaCatAdjust & LX){
        for ( unsigned int freqCat = 0; freqCat < numberFrequenciesXSets; ++freqCat ) {
            sum = 1.0;
            for ( unsigned int i = 0; i < numberFrequencies - 1; ++i ) {
                frequenciesX[freqCat] [i] = SQR(pars[index]);
                ++index;
                sum += frequenciesX[freqCat] [i];
            }
            for ( unsigned int i = 0; i < numberFrequencies - 1; ++i ) {
                frequenciesX[freqCat][i] /= sum;
#ifdef DEBUG3
                assert( !isnan( frequenciesX[freqCat] [i] ) );
                assert( !isinf( frequenciesX[freqCat] [i] ) );
#endif
            }
            frequenciesX[freqCat][numberFrequencies - 1] = 1.0 / sum;
        }
    }
    assert( index == (int)pars.size() );
}


unsigned int Frequencies::getNumberFreeParameters() const{
    unsigned int numberFreeParameters = ( numberFrequencies - 1 ) *
            ( ((gammaCatAdjust & LX) ? numberFrequenciesXSets : 0) +
              ((gammaCatAdjust & L) ? 2 : numberFrequenciesSets) ) ;
    if ( (gammaCatAdjust & L) &&
         (!(gammaCatAdjust&LINK) || (gammaCatAdjust&LINKY)) ){
        if ( gammaCatAdjust & LS ){
            numberFreeParameters+=1;
        }
        if ( gammaCatAdjust & LX ){
            numberFreeParameters+=numberFrequenciesXSets;
        }
    }
    return numberFreeParameters;
}

void Frequencies::printLine( ostream & outputStream ){
    outputStream << setprecision( 4 );
    if ( (gammaCatAdjust & L) &&
         (!(gammaCatAdjust&LINK) || (gammaCatAdjust&LINKY)) ){
        if ( gammaCatAdjust & LS ){
            outputStream << setw( 9 ) << *sF << ' ';
        }
        if ( gammaCatAdjust & LX ){
            vector< double > x = (*xF);
            for( unsigned int i = numberFrequenciesXSets-1; i >0; --i ){
                x[i-1] = x[i-1]*x[i];
            }
            for( unsigned int i = 0; i < numberFrequenciesXSets; ++i ){
                outputStream << setw( 9 ) << x[i] << ' ';
            }
        }
    }
    unsigned int increment = (gammaCatAdjust & L) ? numberFrequenciesSets - 1 : 1;
    for (unsigned int freqCat = 0; freqCat < numberFrequenciesSets; freqCat += increment){
        for ( unsigned int i = 0; i < numberFrequencies; ++i ) {
            outputStream << setw( 9 ) << frequencies[freqCat] [i] << ' ';
        }
    }
    if ( gammaCatAdjust & LX ){
        for (unsigned int freqCat = 0; freqCat < numberFrequenciesXSets; ++freqCat ){
            for ( unsigned int i = 0; i < numberFrequencies; ++i ) {
                outputStream << setw( 9 ) << frequenciesX[freqCat] [i] << ' ';
            }
        }
    }
}


void Frequencies::fromLine( const vector<double>& parameters ){
    int index = -1;
    if ( (gammaCatAdjust & L) &&
         (!(gammaCatAdjust&LINK) || (gammaCatAdjust&LINKY)) ){
        if ( gammaCatAdjust & LS ){
            *sF = parameters[++index];
        }
        if ( gammaCatAdjust & LX ){
            for( unsigned int i = 0; i < numberFrequenciesXSets; ++i ){
                (*xF)[i] = parameters[++index];
            }
            for( unsigned int i = 0; i < numberFrequenciesXSets-1; ++i ){
                (*xF)[i] = (*xF)[i] / (*xF)[i+1];
            }
        }
    }
    double sum;
    unsigned int increment = (gammaCatAdjust & L) ? numberFrequenciesSets - 1 : 1;
    for (unsigned int freqCat = 0; freqCat < numberFrequenciesSets; freqCat += increment){
        sum = 0.0;
        for ( unsigned int i = 0; i < numberFrequencies; ++i ) {
			// Small pseudocount, to prevent later matrix operations failing.
            frequencies[freqCat] [i] = parameters[++index] + 0.000001;
            sum += frequencies[freqCat][i];
        }
        if ( fabs(sum-1.0) > 0.004 ){
            cerr << "WARNING, frequencies do not sum to 1 while "
                 << "reading parameters from a line (" <<sum << ')' << endl;
        }
        for ( unsigned int i = 0; i < numberFrequencies; ++i ) {
            frequencies[freqCat][i] /= sum;
	cout << "frequencies[freqCat][i] = [" << freqCat << "][" << i << "] = " << frequencies[freqCat][i] << endl;
        }
    }
    if ( gammaCatAdjust & LX ){
        for (unsigned int freqCat = 0; freqCat < numberFrequenciesXSets; ++freqCat ){
            sum = 0.0;
            for ( unsigned int i = 0; i < numberFrequencies; ++i ) {
                frequenciesX[freqCat][i] = parameters[++index];
                sum += frequenciesX[freqCat][i];
            }
            if ( fabs(sum-1.0) > 0.004 ){
                cerr << "WARNING, frequencies do not sum to 1 (" << sum
                     << ") while reading parameters from a line" << endl;
            }
            for ( unsigned int i = 0; i < numberFrequencies; ++i ) {
                frequenciesX[freqCat][i] /= sum;
            }
        }
    }
    assert(index == (int)parameters.size()-1);
}

unsigned int Frequencies::getNumberLineParameters() const{
    unsigned int numberLineParameters = numberFrequencies *
            ( ((gammaCatAdjust & LX) ? numberFrequenciesXSets : 0) +
              ((gammaCatAdjust & L) ? 2 : numberFrequenciesSets) );
    if ( (gammaCatAdjust & L) &&
         (!(gammaCatAdjust&LINK) || (gammaCatAdjust&LINKY)) ){
        if ( gammaCatAdjust & LS ){
            numberLineParameters+=1;
        }
        if ( gammaCatAdjust & LX ){
            numberLineParameters+=numberFrequenciesXSets;
        }
    }
    return numberLineParameters;
}



void Frequencies::initialiseMCMC( ParametersSet & parameters, Perturbator* perturbator ){
    char label[50];

    //prepare the update message (base class) sent to the model for an update
    setType( UpdateMessage::PARAM_TYPE );
    setFlag( UpdateMessage::FREQ );

    if (numberFrequenciesSets != 0){
        //register S or X variation parameter(s)
        registerVariationParameter( parameters, perturbator );
        //register frequencies
        char* numberPos = label;
        if (numberFrequenciesSets!=1){
            sprintf ( label, "Frequencies, Set ");
            while (*(++numberPos));
        }
        else{
            sprintf ( label, "Frequencies");
        }
        unsigned int increment = (gammaCatAdjust & L) ? numberFrequenciesSets - 1 : 1;
        for ( unsigned int i = 0; i < numberFrequenciesSets; i += increment){
            if (numberFrequenciesSets!=1){
                sprintf( numberPos, "%d", i+1 );
            }
            if (!parameters.findParameter("Frequencies, prior")){
                parameters["Frequencies, prior"] = "dirichlet(1.0)";
            }
            perturbator->registerFrequencies( *this, UpdateMessage::FREQ, label, &(frequencies[i]),
                   parameters, "Frequencies", 1000.0,
                   "initial Dirichlet tuning parameter");
        }
        //register X frequencies
        if (gammaCatAdjust & LX){
            sprintf ( label, "Frequencies, Set X");
            char* numberPos = label;
            while (*(++numberPos));
            for ( unsigned int i = 0; i < numberFrequenciesXSets; ++i){
                sprintf( numberPos, "%d", i+1 );
                if (!parameters.findParameter("Frequencies X, prior")){
                    if (!parameters.findParameter("Frequencies, prior")){
                        parameters["Frequencies X, prior"] = "dirichlet(1.0)";
                    }
                     parameters["Frequencies X, prior"] = parameters.stringParameter("Frequencies, prior");
                }
                perturbator->registerFrequencies( *this, UpdateMessage::FREQ, label, &(frequenciesX[i]),
                       parameters, "Frequencies X", 1000.0,
                       "initial Dirichlet tuning parameter");
            }
        }
    } // end if numberFrequencies != 0
    //more than one set? smoothing prior....
    if(numberFrequenciesSets>1){
        string stringPrior;
        variationType = NB_TYPE;
        double bound1, bound2, mean1;
        string name = "Frequencies variation prior";
        stringPrior = parameters.stringParameter(name);
        PriorField variationPrior(stringPrior);
        if (variationPrior.name()=="gp"){
            if (variationPrior.getNumberParameters() != 6){
                cerr << "Error, six fields are expected with the \"gp\" frequencies variation prior" << endl;
                exit(EXIT_FAILURE);
            }
            priorParams.resize(numberFrequencies);
            hyperPriorsPert.resize(numberFrequencies);
            hyperPriors.resize(numberFrequencies);
            //priorParams.resize(1);
            //priorParams[0].resize(7);
            //hyperPriorsPert.resize(1);
            //hyperPriors.resize(1);
            //hyperPriors[0].clear();
            //hyperPriorsPert[0].clear();
            invCovMatrix.resize(numberFrequencies);
            for (unsigned int i = 0; i < numberFrequencies; ++i){
                priorParams[i].resize(5);
                hyperPriors[i].clear();
                hyperPriorsPert[i].clear();
                invCovMatrix[i].resize(numberFrequenciesSets,numberFrequenciesSets);
            }
            logdet.resize(numberFrequencies);
            sumAct.resize(numberFrequenciesSets);
            hyperSumAct.clear();
            activationMatrix.resize(numberFrequencies,numberFrequenciesSets);
            bound1 = 0;
            bound2 = +DBL_MAX;
            create( perturbator, name+" theta1 hyperparameter", variationPrior, 0, bound1, bound2, mean1,
                    parameters);
            bound1 = 0;
            bound2 = +DBL_MAX;
            create( perturbator, name+" scale hyperparameter", variationPrior, 1, bound1, bound2, mean1,
                    parameters);
            //theta2 adds a constant component to the regression function.
            //this constant component has a gaussian prior mean 0, variance theta2
            bound1 = 0;
            bound2 = +DBL_MAX;
            create( perturbator, name+" theta2 hyperparameter", variationPrior, 2, bound1, bound2, mean1,
                    parameters);
            //linear trend adds a linear component to the regression function.
            bound1 = -DBL_MAX;
            bound2 = +DBL_MAX;
            create( perturbator, name+" linear trend hyperparameter", variationPrior, 3, bound1, bound2, mean1,
                    parameters);
            bound1 = 0;
            bound2 = +DBL_MAX;
            create( perturbator, name+" jitter hyperparameter", variationPrior, 4, bound1, bound2, mean1,
                    parameters);
            for (unsigned int i =0;i<numberFrequenciesSets; ++i){
                bound1 = -DBL_MAX;
                bound2 = +DBL_MAX;
                createSumAct( perturbator, name+" activation sum hyperparameter", variationPrior, i, bound1, bound2, mean1,
                              parameters);
            }
            variationType = GP;
        }
        if (variationPrior.name()=="uniform"){
            if (variationPrior.getNumberParameters()){
                cerr << "invalid prior for frequencies: " << variationPrior.toString() << endl;
                exit(EXIT_FAILURE);
            }
            variationType = UNIFORM;
        }
        if (variationType==NB_TYPE){
            cerr << "unrecognized prior for frequencies variation: " << variationPrior.name() << endl;
            exit(EXIT_FAILURE);
        }
    }
}

void Frequencies::initialiseML( ParametersSet & parameters ){
    if(numberFrequenciesSets>1){
        variationType = NB_TYPE;
        string name = "Frequencies variation penalty";
        string stringPenalty = parameters.stringParameter(name);
        PriorField variationPrior(stringPenalty);
        if (variationPrior.name()=="gp"){
            if (variationPrior.getNumberParameters() != 6){
                cerr << "Error, six fields are expected with the \"gp\" frequencies variation prior" << endl;
                exit(EXIT_FAILURE);
            }
            priorParams.resize(numberFrequencies);
            hyperPriors.resize(numberFrequencies);
            //priorParams.resize(1);
            //priorParams[0].resize(7);
            //hyperPriors.resize(1);
            //hyperPriors[0].clear();
            invCovMatrix.resize(numberFrequencies);
            for (unsigned int i = 0; i < numberFrequencies; ++i){
                priorParams[i].resize(5);
                hyperPriors[i].clear();
                invCovMatrix[i].resize(numberFrequenciesSets,numberFrequenciesSets);
            }
            logdet.resize(numberFrequencies);
            sumAct.resize(numberFrequenciesSets);
            hyperSumAct.clear();
            activationMatrix.resize(numberFrequencies,numberFrequenciesSets);
            for ( unsigned int i = 0; i < 5; ++i ){
                if (variationPrior.isConstant(i)){
                    for (unsigned int freq = 0; freq < numberFrequencies; ++freq){
                        priorParams[freq][i] = variationPrior.getValue(i);
                    }
                    //priorParams[0][i] = variationPrior.getValue(i);
                }
                else{
                    if ( variationPrior.getHyperPriorField(i).name()!="var" ){
                        cerr << "Error, constant number or var(x) expected in the penalized approach" << endl;
                        exit(EXIT_FAILURE);
                    }
                    if (variationPrior.getHyperPriorField(i).getNumberParameters()!=1){
                        cerr << "Only one parameter is expected in " << variationPrior.getHyperPriorField(i).toString() << endl;
                        exit(EXIT_FAILURE);
                    }
                    for (unsigned int freq = 0; freq < numberFrequencies; ++freq){
                        priorParams[freq][i] = variationPrior.getHyperPriorField(i).getValue(0);
                        hyperPriors[freq].push_back(i);
                    }
                    //priorParams[0][i] = variationPrior.getHyperPriorField(i).getValue(0);
                    //hyperPriors[0].push_back(i);
                }
            }
            for (unsigned int i =0;i<numberFrequenciesSets; ++i){
                if (variationPrior.isConstant(5)){
                    sumAct[i] = variationPrior.getValue(5);
                }
                else{
                    if ( variationPrior.getHyperPriorField(5).name()!="var" ){
                        cerr << "Error, constant number or var(x) expected in the penalized approach" << endl;
                        exit(EXIT_FAILURE);
                    }
                    if (variationPrior.getHyperPriorField(5).getNumberParameters()!=1){
                        cerr << "Only one parameter is expected in " << variationPrior.getHyperPriorField(i).toString() << endl;
                        exit(EXIT_FAILURE);
                    }
                    sumAct[i] = variationPrior.getHyperPriorField(5).getValue(0);
                    hyperSumAct.push_back(i);
                }
            }
            variationType = GP;
        }
        if (variationPrior.name()=="uniform"){
            if (variationPrior.getNumberParameters()){
                cerr << "invalid prior for frequencies: " << variationPrior.toString() << endl;
                exit(EXIT_FAILURE);
            }
            variationType = UNIFORM;
            lnPrior = 0;
        }
        if (variationType==NB_TYPE){
            cerr << "unrecognized penalty for frequencies variation: " << variationPrior.name() << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        variationType = UNIFORM;
        lnPrior = 0;
    }
}


double Frequencies::getLnPenalty( const vector<double> & rates ) {
    if (variationType == GP){
        assert(rates==savedRates);
        array2D<double> testCov = invCovMatrix[1];
        computeCovMatrix();
        assert( invCovMatrix[1] == testCov );
        assert( lnPrior == computeLnPrior() );
    }
    if (lnPrior>100000000.0){
        cerr << "WARNING: rates / freq = " << endl;
        assert(rates.size()==numberFrequenciesSets);
        for ( unsigned int i=0; i < numberFrequenciesSets; ++i ){
            cerr << i << ":   subst. rate=" << rates[i] << ",    freq=[";
            for ( unsigned int j = 0; j < numberFrequencies; ++j ){
                cerr << frequencies[i][j] << ',';
            }
            cerr << ']' << endl;
        }
    }
    return lnPrior;
}

void Frequencies::diffLnPenalty( vector<double>& gradVector ) const{
  //warning, parameters given to the optimizer via getAllPenaltyParameters are : sqr(O0), sqr(O1), sqr(O2), sqr(O3), sqr(O4)
    gradVector.clear();
    if(variationType==GP){
        vector<double> logRates;
        logRates.resize(savedRates.size());
        for (unsigned int k = 0; k< savedRates.size(); ++k){
            if (savedRates[k]<1e-14){
                logRates[k]=log(1e-14);
            }
            else{
                logRates[k]=log(savedRates[k]);
            }
        }

        array2D<double> temp;
        temp.resize(savedRates.size(),savedRates.size());
        array2D<double> prod;
        prod.resize(savedRates.size(),savedRates.size());
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            if (logdet[freq]>=FLT_MAX){
                for (unsigned int paramId = 0; paramId < hyperPriors[freq].size(); ++paramId){
                    gradVector.push_back(0.0);
                }
            }
            else{
                // dL/dO = -.5 * trace(C-1 * dC/dO) + .5 * freq' * (C-1 * dC/dO) * C-1 * freq
                for (unsigned int paramId = 0; paramId < hyperPriors[freq].size(); ++paramId){
                    for(unsigned int lineTemp = 0; lineTemp < savedRates.size(); ++lineTemp){
                        for(unsigned int colTemp = 0; colTemp < savedRates.size(); ++colTemp){
                            temp(lineTemp,colTemp) = 0.0;
                        }
                    }
                    for(unsigned int line = 0; line < savedRates.size(); ++line){
                        if (hyperPriors[freq][paramId]==4){
                            for(unsigned int col = 0; col < savedRates.size(); ++col){
                                temp(line,col) = invCovMatrix[freq](line,col) * 4 * sqrt(priorParams[freq][4]) * priorParams[freq][4];
                            }
                        }
                        else{
                            for(unsigned int col = 0; col < savedRates.size(); ++col){
                                double dCdO = 0.0;
                                double sqr = 0.0;
                                switch(hyperPriors[freq][paramId]){
                                  case 0:
                                    sqr = logRates[col] - logRates[line];
                                    sqr = sqr*sqr;
                                    dCdO = 4 * sqrt(priorParams[freq][0]) * priorParams[freq][0] * exp(-sqr*(priorParams[freq][1] * priorParams[freq][1]));
                                  break;
                                  case 1:
                                    sqr = logRates[col] - logRates[line];
                                    sqr = sqr*sqr;
                                    dCdO = priorParams[freq][0] * priorParams[freq][0];
                                    dCdO = dCdO * -sqr * 4 * sqrt(priorParams[freq][1]) * priorParams[freq][1];
                                    dCdO = dCdO * exp(-sqr*(priorParams[freq][1] * priorParams[freq][1]));
                                  break;
                                  case 2:
                                    dCdO = 4 * sqrt(priorParams[freq][2]) * priorParams[freq][2];
                                  break;
                                  case 3:
                                    dCdO = 4 * logRates[col]*logRates[line] * sqrt(priorParams[freq][3]) * priorParams[freq][3];
                                  break;
                                  default:
                                    assert(0);
                                }
                                for(unsigned int lineTemp = 0; lineTemp < savedRates.size(); ++lineTemp){
                                    temp(lineTemp,col) += dCdO * invCovMatrix[freq](lineTemp,line);
                                }
                            }
                        }
                    }

                    double res = 0.0;
                    for(unsigned int i = 0; i < savedRates.size(); ++i){
                        res += temp(i,i);
                    }
                    res = -.5 * res;
                    prod = operator*(temp,invCovMatrix[freq]);
                    for (unsigned int i =  0; i < savedRates.size(); ++i){
                        double tmp = 0.0;
                        for (unsigned int j =  0; j < savedRates.size(); ++j){
                            tmp += prod(i,j)*activationMatrix(freq,j);
                        }
                        res += .5 * activationMatrix(freq,i) * tmp;
                    }
                    gradVector.push_back(res);
                }
            }
        }
    }
    //assert(gradVector.size()==getNumberPenaltyParameters());
}

void Frequencies::create( Perturbator* perturbator,
                      const string& name, const PriorField& prior,
                      unsigned int priorField, double& min, double& max, double& mean,
                      ParametersSet& parameters ){
    assert(priorField<priorParams[0].size());
    if (prior.isConstant(priorField)){
        if ( ( prior.getValue(priorField) >= min ) &&
             ( prior.getValue(priorField) <= max ) ){
            min = prior.getValue(priorField);
            max = min;
            mean = min;
        }
        else{
            cerr << name << " invalid: " << prior.getHyperPriorField(priorField).toString() << endl;
            exit(EXIT_FAILURE);
        }
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            priorParams[freq][priorField]=mean;
        }
    }
    else{
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            hyperPriors[freq].push_back(priorField);
            PerturbatorGaussParameter* pert =
                perturbator->registerGaussParameter( *this, UpdateMessage::HYPER_PARAM, name,
                    &(priorParams[freq][priorField]),
                    prior.getHyperPriorField(priorField), parameters, name,
                    .1, .001, 10.0, "initial step" );
            hyperPriorsPert[freq].push_back(pert);
            pert->getExtrema( min, max );
            mean = pert->getMean();
            priorParams[freq][priorField]=mean;
            //this call is vital to refresh the prior computed for the old value
            pert->invalidate();
        }
    }
}

void Frequencies::createSumAct( Perturbator* perturbator,
                      const string& name, const PriorField& prior,
                      unsigned int priorField, double& min, double& max, double& mean,
                      ParametersSet& parameters){
    assert(prior.getNumberParameters()==6);
    if (prior.isConstant(5)){
        if ( ( prior.getValue(5) >= min ) &&
             ( prior.getValue(5) <= max ) ){
            min = prior.getValue(5);
            max = min;
            mean = min;
        }
        else{
            cerr << name << "invalid: " << prior.toString() << endl;
            exit(EXIT_FAILURE);
        }
        sumAct[priorField]=mean;
    }
    else{
        hyperSumAct.push_back(priorField);
        PerturbatorGaussParameter* pert =
            perturbator->registerGaussParameter( *this, UpdateMessage::HYPER_PARAM, name,
                &(sumAct[priorField]),
                prior.getHyperPriorField(5), parameters, name,
                .3, .001, 10.0, "initial step" );
        hyperSumActPert.push_back(pert);
        pert->getExtrema( min, max );
        mean = pert->getMean();
        sumAct[priorField]=mean;
        //this call is vital to refresh the prior computed for the old value
        pert->invalidate();
    }
}

void Frequencies::getAllPriorParameters( vector<double>& params ) const{
    if(variationType==GP){
        params.resize(getNumberPriorParameters());
        vector<double>::iterator iter = params.begin();
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            for (unsigned int i = 0; i< hyperPriors[freq].size(); ++i){
                *iter = priorParams[freq][hyperPriors[freq][i]];
                ++iter;
            }
        }
        for (unsigned int i = 0; i< hyperSumAct.size(); ++i){
            assert(hyperSumAct[i]==i);
            *iter = sumAct[hyperSumAct[i]];
            ++iter;
        }
    }
    else{
        params.clear();
    }
}

void Frequencies::setAllPriorParameters( const vector<double>& params){
    if(variationType==GP){
        vector<double>::const_iterator iter = params.begin();
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            for (unsigned int i = 0; i< hyperPriors[freq].size(); ++i){
                priorParams[freq][hyperPriors[freq][i]] = *iter;
                ++iter;
                hyperPriorsPert[freq][i]->invalidate();
            }
        }
        for (unsigned int i = 0; i< hyperSumAct.size(); ++i){
            sumAct[hyperSumAct[i]] = *iter;
            hyperSumActPert[i]->invalidate();
            ++iter;
        }
        computeCovMatrix();
        lnPrior = computeLnPrior();
    }
    else{
        assert(params.size()==0);
    }
}

unsigned int Frequencies::getNumberPriorParameters() const{
    if(variationType==GP){
        unsigned int nb = hyperSumAct.size();
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            nb += hyperPriors[freq].size();
        }
        return nb;
    }
    else{
        return 0;
    }
}

double Frequencies::getLnPrior( const vector<double> & rates ){
    if (variationType == GP){
        //to be removed, leave it a while for safety...
        assert(savedRates==rates);
        array2D<double> testCov = invCovMatrix[1];
        computeCovMatrix();
        assert( invCovMatrix[1] == testCov );
        assert( lnPrior == computeLnPrior() );
    }
    return lnPrior;
}

void Frequencies::getAllPenaltyParameters( vector<double>& params ) const{
    if(variationType==GP){
        params.resize(getNumberPenaltyParameters());
        vector<double>::iterator iter = params.begin();
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            for (unsigned int i = 0; i< hyperPriors[freq].size(); ++i){
                *iter = sqrt(priorParams[freq][hyperPriors[freq][i]]);
                ++iter;
            }
        }
        for (unsigned int i = 0; i< hyperSumAct.size(); ++i){
            assert(hyperSumAct[i]==i);
            *iter = sumAct[hyperSumAct[i]];
            ++iter;
        }
    }
    else{
        params.clear();
    }
}

void Frequencies::setAllPenaltyParameters( const vector<double>& params ){
    assert(params.size()==getNumberPenaltyParameters());
    if(variationType==GP){
        vector<double>::const_iterator iter = params.begin();
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            for (unsigned int i = 0; i< hyperPriors[freq].size(); ++i){
                priorParams[freq][hyperPriors[freq][i]] = SQR(*iter);
                ++iter;
            }
        }
        for (unsigned int i = 0; i< hyperSumAct.size(); ++i){
            sumAct[hyperSumAct[i]] = *iter;
            ++iter;
        }
    }
    else{
        assert(params.size()==0);
    }
}


unsigned int Frequencies::getNumberPenaltyParameters() const{
    if(variationType==GP){
        unsigned int nb = hyperSumAct.size();
        for ( unsigned int freq = 0; freq < numberFrequencies; ++freq ){
            nb += hyperPriors[freq].size();
        }
        return nb;
    }
    else{
        return 0;
    }
}


void Frequencies::computeCovMatrix(){
    assert(savedRates.size()==numberFrequenciesSets);
    vector<double> logRates;
    logRates.resize(savedRates.size());
    for (unsigned int k = 0; k< savedRates.size(); ++k){
        if (savedRates[k]<1e-14){
            logRates[k]=log(1e-14);
        }
        else{
            logRates[k]=log(savedRates[k]);
        }
    }
    
    array2D<double> sumCov(savedRates.size(),savedRates.size());
    for (unsigned int k = 0; k< savedRates.size(); ++k){
        for (unsigned int i = 0; i < savedRates.size(); ++i){
            sumCov(i,k)=0.0;
        }
    }
    //for each frequency
    for (unsigned int freq = 0; freq < numberFrequencies; ++freq){
        double pr0sqr = priorParams[freq][0] * priorParams[freq][0];
        double pr1sqr = priorParams[freq][1] * priorParams[freq][1];
        if(!pr1sqr){
            cerr << "Error: scale prior == 0.0" << endl;
            pr1sqr = 1e-14;
        }
        double pr2sqr = priorParams[freq][2] * priorParams[freq][2];
        double pr3sqr = priorParams[freq][3] * priorParams[freq][3];
        double pr4sqr = priorParams[freq][4] * priorParams[freq][4];
        //fill invCovMatrix with the values of the covariance matrix
        for (unsigned int k = 0; k< savedRates.size(); ++k){
            for (unsigned int i = 0; i< k; ++i){
                double sqr = logRates[k] - logRates[i];
                sqr = sqr*sqr;
                invCovMatrix[freq](i,k) = pr0sqr *exp(-sqr*pr1sqr) + pr2sqr + pr3sqr*logRates[k]*logRates[i];
                sumCov(i,k) += invCovMatrix[freq](i,k); //add to the sum
                invCovMatrix[freq](k,i) = invCovMatrix[freq](i,k);
            }
            //(k,k) diagonal element add the jitter
            invCovMatrix[freq](k,k) = pr0sqr + pr2sqr + pr3sqr*logRates[k]*logRates[k] + pr4sqr;
            sumCov(k,k) += invCovMatrix[freq](k,k); //add to the sum
            assert(!isnan(invCovMatrix[freq](k,k)));
            assert(!isinf(invCovMatrix[freq](k,k)));
        }
        
        logdet[freq] = inverseSPD2( invCovMatrix[freq], savedRates.size() );
        if(isnan(logdet[freq])  || isinf(logdet[freq]) || (logdet[freq]<-10000.0) ){
            cerr << "Err, thresholding the matrix determinant" << endl;
            logdet[freq]=FLT_MAX;
        }
    }
    //finish filling sumCov
    for ( unsigned int k = 0; k< savedRates.size(); ++k ){
        for (unsigned int i = 0; i< k; ++i){
            sumCov(k,i) = sumCov(i,k);
        }
    }
    logdetsum = inverseSPD2( sumCov, savedRates.size() );
    if(isnan(logdetsum)  || isinf(logdetsum) || (logdetsum<-10000.0) ){
        cerr << "Err, thresholding the matrix determinant" << endl;
        logdetsum=-FLT_MAX;
    }
}


double Frequencies::computeLnPrior(){
    double ret = 0.0;
    if (variationType==GP){
        //f(C) = exp(a(C)) / sum ( exp(a(C')) )
        //we compute the log of the denominator:      log(  sum( exp(a(C')) ) ) knowing that sum(a(C')) = sumAct
        assert(savedRates.size()==numberFrequenciesSets);
        for (unsigned int rateCat = 0; rateCat< numberFrequenciesSets; ++rateCat){
            double logSumExpActivation = 0.0;
            for (unsigned int j = 0; j < numberFrequencies; ++j){
                logSumExpActivation -= log ( frequencies[rateCat][j] );
            }
            logSumExpActivation = (logSumExpActivation+sumAct[rateCat])/numberFrequencies;
            for (unsigned int j = 0; j < numberFrequencies; ++j){
                activationMatrix(j,rateCat) = logSumExpActivation + log(frequencies[rateCat][j]);
            }
        }

        //see Neal (1997) and Gowri-Shankar's PhD thesis for the maths.
        //-log(Z) = - .5 * numberCategories*(numberFreq-1)*log(2 pi)
        //          - .5 * sum_{i=numberFreq} logdet(C_i)
        //          + .5 * logdet( sum_{i=numberFreq} C_i )
        ret -= (double)(numberFrequenciesSets * (numberFrequencies-1) ) * M_LN_SQRT_2PI;
        ret += .5*logdetsum;
        //now prob(f(C)) can be computed according to prob(a(C))
        for (unsigned int freq = 0; freq < numberFrequencies; ++freq){
            if (logdet[freq]>=FLT_MAX){
                return -FLT_MAX;
            }
            ret -= .5 * logdet[freq];
            for (unsigned int rateCat = 0; rateCat< numberFrequenciesSets; ++rateCat){
                double mom = 0.0;
                for (unsigned int r = 0; r< numberFrequenciesSets; ++r){
                    mom += invCovMatrix[freq](rateCat,r) * activationMatrix(freq,r);
                }
                ret -= .5 * mom * activationMatrix(freq,rateCat);
            }
        }
    }
    assert(!isnan(ret));
    assert(!isinf(ret));
    return ret;
}

void Frequencies::update(UpdateMessage* subject){
    assert(subject->hasType(UpdateMessage::PARAM_TYPE));
    if (subject->hasFlag(UpdateMessage::HYPER_PARAM)){
#ifdef DEBUG2
        cout << "frequencies received a hyperparam message" << endl;
#endif
        if (!subject->hasFlag(UpdateMessage::PRIOR_FLAG)){
            setFlag( PRIOR_FLAG );
            computeCovMatrix();
            lnPrior = computeLnPrior();
            notify();
            unsetFlag( PRIOR_FLAG );
        }
    }
    else{
#ifdef DEBUG2
        cout << "frequencies received a modif update message" << endl;
#endif
        if (gammaCatAdjust & L) updateLinear();
        lnPrior = computeLnPrior();
        notify();
    }
}
