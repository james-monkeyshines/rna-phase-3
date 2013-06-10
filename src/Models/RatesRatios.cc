#include "Models/RatesRatios.h"

#include <iostream>
#include <iomanip>

#include "Util/ParametersSet.h"
#include "Models/Perturbator.h"
#include "Models/PerturbatorGaussParameter.h"
#include "Models/MatrixModel.h"

#define SQR(n) (n*n)

RatesRatios::RatesRatios( ParametersSet& parameters, MatrixModel * model ){
    initPrimitive( parameters, model->getNumberGammaCategories(),
                        model->getInvariant() );
}


void RatesRatios::initPrimitive( ParametersSet& parameters,
          unsigned int numberGammaCategories, unsigned int invariant ){

    this->invariant = invariant;
    this->numberGammaCategories = numberGammaCategories;

    //initialisation: no adjustment of rates by default
    gammaCatAdjust = 0;
    numberRatesRatiosXSets = 0;
    sR = NULL;
    xR = NULL;
    ratesRatiosX = NULL;

    //default value for the number of rates matrices equals to 1
    //model without rates (numberRatesRatiosPerSet = 0)
    //MUST make sure this value is filled with 0 before construction
    if ( parameters.findParameter( "Number of rates ratios sets" ) ){
        numberRatesRatiosSets = parameters.intParameter( "Number of rates ratios sets" );
    }
    else{
        numberRatesRatiosSets = 1;
    }

    string ratesAdjustment;
    if ( parameters.findParameter( "Rates adjustment" ) ){
        if (numberRatesRatiosSets==0){
            cerr << "No \"Rates adjustment\" possible for a model without free rates parameters" << endl;
            exit(EXIT_FAILURE);
        }
        if (numberRatesRatiosSets>1){
            cerr << "Do not use \"Number of rates ratios sets\" and \"Rates adjustment\" simultaneously" << endl;
            exit(EXIT_FAILURE);
        }
        ratesAdjustment = parameters.stringParameter("Rates adjustment");
        if ( ratesAdjustment == "linear" ){
            gammaCatAdjust = L;
        }
        else if ( ratesAdjustment=="linearS" ){
            gammaCatAdjust = LS | L;
            if ( ( !parameters.findParameter( "Link S parameter" ) ) ||
                 ( !parameters.boolParameter("Link S parameter") ) ){
                sR = new double;
                *sR = 0.5;
            }
            else{
                gammaCatAdjust = gammaCatAdjust | LINK;
            }
        }
        else if ( ratesAdjustment=="linearX" ){
            gammaCatAdjust = LX | L;
            numberRatesRatiosXSets = parameters.intParameter("Number of X rates ratios sets");
            if (numberRatesRatiosXSets < 1){
                cerr << "Invalid number of X rates sets" << endl;
                exit(EXIT_FAILURE);
            }
            if ( ( !parameters.findParameter( "Link X parameter" ) ) ||
                 ( !parameters.boolParameter("Link X parameter") ) ){
                xR = new vector<double>;
                xR->resize( numberRatesRatiosXSets );
                ratesRatiosX = new vector< double >[numberRatesRatiosXSets];
                if ( gammaCatAdjust & LX ){
                    for ( unsigned int i = 0; i < numberRatesRatiosXSets; ++i ){
                        (*xR)[i] = .5;
                    }
                }
            }
            else{
                gammaCatAdjust = gammaCatAdjust | LINK;
            }
        }
        else{
            cerr << "Unrecognized value \"" << ratesAdjustment << "\" for the parameters "
                 << "\"Rates adjustment\"... abort" << endl;
            exit(EXIT_FAILURE);
        }
    }

    //consistency checking at that point, in case of identifiability issue
    if (numberGammaCategories){
        //at most numberGammaCategories rate matrices
        if (numberRatesRatiosSets > numberGammaCategories){
            cerr << "Invalid number of rates ratios sets (" << numberRatesRatiosSets
                 << "). The expected number is in [1.." << numberGammaCategories << "]." << endl;
            exit(EXIT_FAILURE);
        }
        if ( ( numberGammaCategories == 2 ) && (gammaCatAdjust & LS)
             && !(gammaCatAdjust & LINK) ){
                cerr << "Identifiability issue, with 2 rates matrices (for the two gamma categories),"
                     << " you cannot evaluate sR" << endl;
                exit(EXIT_FAILURE);
        }
        //if unlink X, Xrates > numberGammaCategories is not allowed, Xrates > (numberGammaCategories - 2)/2 not encouraged
        if ( gammaCatAdjust & LX ){
            if ( ( numberRatesRatiosXSets > numberGammaCategories )  && !(gammaCatAdjust & LINK) ){
                cerr << "Identifiability issue, \"number of X rates categories\" > \"Number of gamma categories\"" << endl;
                exit(EXIT_FAILURE);
            }
            if ( ( numberRatesRatiosXSets > (numberGammaCategories - 2)/2 ) && !(gammaCatAdjust & LINK) ){
                cerr << "WARNING: You probably should not be using more than "
                     << (numberGammaCategories - 2)/2 << " X rates categories." << endl;
            }
        }
    }
    else{
        //only one rate matrix allowed if we do not use a gamma
        //distribution
        if ( numberRatesRatiosSets > 1 ){
            cerr << "No more than 1 rates matrix without a discrete gamma model. Exit..." << endl;
            exit(EXIT_FAILURE);
        }
        // no linear rates, linearS/X rates allowed without gamma model.
        if ( gammaCatAdjust & L ){
            cerr << "You must use a discrete gamma model when turning on \"Rates adjustment\". Exit..." << endl;
            exit(EXIT_FAILURE);
        }
    }

/** Set up rate ratios without the discrete gamma model **************************************************/
    if ( !numberGammaCategories ) {
        // if the model uses some substitution rates
        if (numberRatesRatiosSets != 0){
            assert(numberRatesRatiosSets==1); //check already done
            ratesCategory.resize( 1 );
            ratesCategory[0] = 0;
            ratesRatios = new vector< double >[1];
        }
        else{
            ratesCategory.resize( 0 );
        }
    }
/** Set up rate ratios with the discrete gamma model **************************************************/
    else {
        // if the model uses some frequencies
        if ( numberRatesRatiosSets != 0 ){
            ratesCategory.resize( numberGammaCategories );
            if ( gammaCatAdjust & L ){
                numberRatesRatiosSets = numberGammaCategories;
            }
            //if no indication provided initialise the ratesCategory vector
            //used to match a rate category with a matrix
            if ( !parameters.findParameter( "Rates for category 1" ) ){
                //allow only two simple case only one matrice for all the categories
                // or one matrice for each category
                if( ( numberRatesRatiosSets != 1 ) && ( numberRatesRatiosSets != numberGammaCategories ) ){
                    cerr << "Since \"Number of rates ratios sets\" is not equal to " << numberGammaCategories
                         << ". You have to specify the set to be used with each gamma category." << endl;
                    cerr << "Use the fields \"Rates for category 1\", \"Rates for category 2\", ..." << endl;
                    exit(EXIT_FAILURE);
                }
                if ( numberRatesRatiosSets == 1 ) {
                    for ( unsigned int i = 0; i < numberGammaCategories; ++i ) {
                        this->ratesCategory[i] = 0;
                    }
                }
                if ( numberRatesRatiosSets == numberGammaCategories ) {
                    for ( unsigned int i = 0; i < numberGammaCategories; ++i ) {
                        this->ratesCategory[i] = i;
                    }
                }
            }
            else {
                if (gammaCatAdjust & L){
                    cerr << "Do not specify the rates ratios sets for each gamma "
                         << "category when using the field \"Rates adjustment\"" << endl;
                    exit(EXIT_FAILURE);
                }
                int mat = 0;
	            char label[25];
                sprintf( label, "Rates for category ");
                char* numberPos = label;
                while (*(++numberPos));
                for ( unsigned int i = 0; i < numberGammaCategories; ++i ) {
                    sprintf( numberPos, "%d", i + 1);
                    ratesCategory[i] = parameters.intParameter( label )-1;
                    if ( ratesCategory[i] < mat ) {
                        cout << "WARNING ! Strange allocation of a rate matrix for "
                        << "the rate category " << i+1 << endl;
                    }
                    mat = ratesCategory[i];
                    if ( ( mat < 0 ) || ( (unsigned int)mat >= numberRatesRatiosSets ) ) {
                        cerr << "ERROR ! assigning rates set " << mat+1
                        << "for the rate category " << i + 1 << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
            ratesRatios = new vector< double >[numberRatesRatiosSets];
        }
        else{
            ratesCategory.resize( 0 );
        }
    } // end with discreteGamma
    lnPrior = 0.0;
}

void RatesRatios::initialisation( unsigned int numberRatesRatios ){
    vector< double > rates;
    rates.resize(numberRatesRatios);
    for ( unsigned int i = 0; i <  numberRatesRatios; ++i ) {
        rates[i] = 1.0;
    }
    initialisation( rates );
}

void RatesRatios::initialisation( const vector<double>& rates ){
    this->numberRatesRatios = rates.size();
    if (numberRatesRatiosSets){
        for ( unsigned int i = 0; i < numberRatesRatiosSets; ++i ) {
            ratesRatios[i] = rates;
        }
        if (gammaCatAdjust & LX){
            for ( unsigned int i = 0; i < numberRatesRatiosXSets; ++i ) {
                ratesRatiosX[i] = rates;
            }
        }
    }
    else{
   //     ratesRatios = NULL;
   //     ratesRatios = new vector<double>;
   //     *ratesRatios = rates;
    }
}


void RatesRatios::updateLinear(){
    if (gammaCatAdjust & LX){
        unsigned int x = 0;
        double highRate = savedRates[0];
        double lowRate = highRate;
        vector<double>* highRatesRatios = &(ratesRatios[0]);
        vector<double>* lowRatesRatios = highRatesRatios;
        double ratio;
        for ( unsigned int i = 1; i < numberRatesRatiosSets-1; ++i ){
            while ( savedRates[i] >= highRate ){
                lowRate = highRate;
                lowRatesRatios = highRatesRatios;
                if ( x == numberRatesRatiosXSets ){
                    highRate = savedRates.back();
                    highRatesRatios = &(ratesRatios[numberRatesRatiosSets-1]);
                }
                else{
                    highRate = (*xR)[x];
                    for (unsigned int j = x+1; j < numberRatesRatiosXSets; ++j){
                        highRate = highRate * (*xR)[j];
                    }
                    highRate = highRate * (savedRates.back()-savedRates.front()) + savedRates.front();
                    highRatesRatios = &(ratesRatiosX[x]);
                    ++x;
                }
            }
            ratio = (savedRates[i]-lowRate)/(highRate - lowRate);
            for ( unsigned int j = 0; j < numberRatesRatios; ++j ){
                ratesRatios[i][j] = (*lowRatesRatios)[j] + ( (*highRatesRatios)[j] - (*lowRatesRatios)[j] ) * ratio;
            }
        }
    }
    //else linearS
    else{
        //fill the frequencies between the lowest and the highest category
        double delta = savedRates.back() - savedRates.front();
        double ratio;
        for ( unsigned int i = 1; i < numberRatesRatiosSets-1; ++i ){
            ratio = (savedRates[i]-savedRates.front())/delta;
            //if linearS
            if (gammaCatAdjust & LS){
                ratio = pow( ratio, *sR );
            }
            for ( unsigned int j = 0; j < numberRatesRatios; ++j ){
                ratesRatios[i][j] = ratesRatios[0][j]+(ratesRatios[numberRatesRatiosSets-1][j]-ratesRatios[0][j])*ratio;
            }
        }
    }
}

void RatesRatios::prepare( const vector<double> & rates ){
    savedRates = rates;
    if (gammaCatAdjust & L) updateLinear();
//    computeCovMatrix();
//    lnPrior = computeLnPrior();
}

double* RatesRatios::linkS(){
    assert( gammaCatAdjust & LINK );
    assert( !(gammaCatAdjust & LINKN) );
    assert( !(gammaCatAdjust & LINKY) );
    assert( gammaCatAdjust & LS );
    gammaCatAdjust = gammaCatAdjust | LINKY;

    sR = new double;
    *sR = .5;
    return sR;
}

void RatesRatios::assignS( double* sR ){
    assert( gammaCatAdjust & LINK );
    assert( !(gammaCatAdjust & LINKN) );
    assert( !(gammaCatAdjust & LINKY) );
    assert( gammaCatAdjust & LS );

    gammaCatAdjust = gammaCatAdjust | LINKN;
    this->sR = sR;
}

vector<double>* RatesRatios::linkX(){
    assert( gammaCatAdjust & LINK );
    assert( !(gammaCatAdjust & LINKN) );
    assert( !(gammaCatAdjust & LINKY) );
    assert( gammaCatAdjust & LX );

    gammaCatAdjust = gammaCatAdjust | LINKY;
    xR = new vector<double>;
    xR->resize( numberRatesRatiosXSets );
    ratesRatiosX = new vector< double >[numberRatesRatiosXSets];
    for ( unsigned int i = 0; i < numberRatesRatiosXSets; ++i ){
        (*xR)[i] = .5;
    }
    return xR;
}

void RatesRatios::assignX( vector<double>* xR ){
    assert( gammaCatAdjust & LINK );
    assert( !(gammaCatAdjust & LINKN) );
    assert( !(gammaCatAdjust & LINKY) );
    assert( gammaCatAdjust & LX );

    if( xR->size() != numberRatesRatiosXSets ){
        cerr << "Cannot link two X parameters, sizes are different "
             << numberRatesRatiosXSets << "<>" << xR->size() << endl;
        exit(EXIT_FAILURE);
    }
    gammaCatAdjust = gammaCatAdjust | LINKN;
    this->xR = xR;
    ratesRatiosX = new vector< double >[numberRatesRatiosXSets];
}


void RatesRatios::getModelParameters( ParametersSet & parameters ) const{
    ParametersSet* currentCategory;
    char tmpString[40];
    char dblString[40];

    if (numberRatesRatiosSets>=1){
        if (numberRatesRatiosSets>=2){
            sprintf(tmpString,"%d",numberRatesRatiosSets);
            parameters("RATESRATIOS")["Number of rates ratios sets"] = tmpString;
            if (gammaCatAdjust & LS){
                sprintf(dblString, "%.8f", *sR);
                parameters("RATESRATIOS")["Variation parameter s"] = dblString;
            }
            unsigned int increment = (gammaCatAdjust&L) ? numberRatesRatiosSets - 1 : 1;
            for ( unsigned int i = 0; i < numberRatesRatiosSets; i += increment ){
                sprintf(tmpString,"RATESRATIOSCATEGORY%d",i+1);
                currentCategory = &parameters("RATESRATIOS")(tmpString);
                for (unsigned int j = 0; j < numberRatesRatios; ++j){
                    sprintf(dblString, "%.8f",ratesRatios[i][j]);
                    sprintf(tmpString,"Rate ratio %d",j+1);
                    (*currentCategory)[tmpString] = dblString;
                }
            }
            if (gammaCatAdjust & LX){
                double ratio = 1.0;
                for ( int i = numberRatesRatiosXSets-1; i >= 0; --i ){
                    sprintf(tmpString,"RATESRATIOSCATEGORYX%d",i+1);
                    currentCategory = &parameters("RATESRATIOS")(tmpString);
                    ratio = ratio * (*xR)[i];
                    sprintf(dblString, "%.8f", ratio);
                    (*currentCategory)["Variation parameter x"] = dblString;
                    for (unsigned int j = 0; j < numberRatesRatios; ++j){
                        sprintf(dblString, "%.8f",ratesRatiosX[i][j]);
                        sprintf(tmpString,"Rate ratio %d",j+1);
                        (*currentCategory)[tmpString] = dblString;
                    }
                }
            }
            if (!gammaCatAdjust&L){
                for ( unsigned int i = 0; i < ratesCategory.size(); ++i ){
                    sprintf(tmpString,"RATESCATEGORY%d",i+1+invariant);
                    currentCategory = &parameters("RATESCATEGORIES")(tmpString);
                    sprintf(dblString,"%d",ratesCategory[i]+1);
                    (*currentCategory)["Rates ratios set for the category"] = dblString;
                }
            }
        }
        else{
            for (unsigned int j = 0; j < numberRatesRatios; ++j){
                sprintf(dblString, "%.8f",ratesRatios[0][j]);
                sprintf(tmpString,"Rate ratio %d",j+1);
                parameters("RATESRATIOS")[tmpString] = dblString;
            }
        }
    }
}

void RatesRatios::setModelParameters( ParametersSet & parameters ){
    ParametersSet* currentCategory;
    char tmpString[40];

    if ( numberRatesRatiosSets >= 1){
        if ( numberRatesRatiosSets >= 2 ){
            if ( numberRatesRatiosSets != (unsigned int)
                parameters("RATESRATIOS").intParameter("Number of rates ratios sets") ){
                cerr << "The number of rates ratios sets does not match in the model file." << endl;
                exit(EXIT_FAILURE);
            }
            if (gammaCatAdjust & LS){
                *sR = parameters("RATESRATIOS").doubleParameter("Variation parameter s");
            }
            unsigned int increment = (gammaCatAdjust&L) ? numberRatesRatiosSets - 1 : 1;
            for ( unsigned int i = 0; i < numberRatesRatiosSets; i += increment ){
                sprintf(tmpString,"RATESRATIOSCATEGORY%d",i+1);
                currentCategory = &parameters("RATESRATIOS")(tmpString);
                for (unsigned int j = 0; j < numberRatesRatios; ++j){
                    sprintf(tmpString,"Rate ratio %d",j+1);
                    ratesRatios[i][j] = currentCategory->doubleParameter(tmpString);
                }
            }
            if (gammaCatAdjust & LX){
                for ( unsigned int i = 0; i < numberRatesRatiosXSets; ++i ){
                    sprintf(tmpString,"RATESRATIOSCATEGORYX%d",i+1);
                    currentCategory = &parameters("RATESRATIOS")(tmpString);
                    (*xR)[i] = currentCategory->doubleParameter("Variation parameter x");
                    for (unsigned int j = 0; j < numberRatesRatios; ++j){
                        sprintf(tmpString,"Rate ratio %d",j+1);
                        ratesRatiosX[i][j] = currentCategory->doubleParameter(tmpString);
                    }
                }
                for ( unsigned int i = 0; i < numberRatesRatiosXSets-1; ++i ){
                    (*xR)[i] = (*xR)[i] / (*xR)[i+1];
                }
            }
            if (!gammaCatAdjust&L){
                for ( unsigned int i = 0; i < ratesCategory.size(); ++i ){
                    sprintf(tmpString,"RATESCATEGORY%d",i+1+invariant);
                    currentCategory = &parameters("RATESCATEGORIES")(tmpString);
                    ratesCategory[i] = currentCategory->intParameter("Rates ratios set for the category")  - 1;
                    assert( (ratesCategory[i]>=0) &&
                                  ((unsigned int)ratesCategory[i]<numberRatesRatiosSets) );
                }
            }
       }
        else{
            for (unsigned int j = 0; j < numberRatesRatios; ++j){
                sprintf(tmpString,"Rate ratio %d",j+1);
                ratesRatios[0][j] = parameters("RATESRATIOS").doubleParameter(tmpString);
            }
        }
    }
}


void RatesRatios::printParameters( ostream & outputStream ) const{
    if( (!(gammaCatAdjust & LS) && !(gammaCatAdjust & LX)) || (gammaCatAdjust & LINKN) ){
        return;
    }
    if (!(gammaCatAdjust & LINK)){
        outputStream << "Rates ratios, ";
    }
    //if S
    if (gammaCatAdjust & LS){
        outputStream << "S variation parameter = " << *sR << endl;
    }
    //if X
    if (gammaCatAdjust & LX){
        outputStream << "X variation parameters = ";
        vector< double > x = *xR;
        for (unsigned int i = numberRatesRatiosXSets-1; i > 0; --i){
            x[i-1]=x[i]*x[i-1];
        }
        for (unsigned int i = 0; i < numberRatesRatiosXSets; ++i){
            if (i) {
                outputStream << ", ";
            }
            outputStream << x[i];
        }
        outputStream << endl;
    }
}

void RatesRatios::registerVariationParameter( ParametersSet & parameters, Perturbator* perturbator ){
    char label[70];
    char* numberPos;

    double extrMin, extrMax;
    //register S or X variation parameter(s)
    if ( !(gammaCatAdjust & LINK) ){
        if (gammaCatAdjust & LS){
            if (!parameters.findParameter("Rates ratios variation parameter S, prior")){
                parameters["Rates ratios variation parameter S, prior"] = "uniform(0,1.0)";
            }
            PerturbatorGaussParameter* pert =
                perturbator->registerGaussParameter( *this, UpdateMessage::RATE,
                           "Rates ratios variation parameter S", sR,
                           parameters, "Rates ratios variation parameter S",
                           .05, .001, 0.5, "initial step" );
            pert->getExtrema( extrMin, extrMax );
            if (extrMin<0.0){
                cerr << "invalid prior for the rates variation parameter S: "
                     << parameters["Rates ratios variation parameter S, prior"] << endl;
                exit(EXIT_FAILURE);
            }
        }
        if (gammaCatAdjust & LX){
            sprintf ( label, "Rates Ratios variation parameter X");
            numberPos = label;
            while (*(++numberPos));
            for (unsigned int i = 0; i < numberRatesRatiosXSets; ++i){
                sprintf( numberPos, "%d", i+1 );
                if (!parameters.findParameter(string(label)+", prior")){
                    parameters[string(label)+", prior"] = "uniform(0,1.0)";
                }
                PerturbatorGaussParameter* pert =
                    perturbator->registerGaussParameter( *this, UpdateMessage::RATE, label, &((*xR)[i]),
                               parameters, "Rates ratios variation parameter X",
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
                    perturbator->registerGaussParameter( *this, UpdateMessage::RATE,
                               "Variation parameter S", sR,
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
                for (unsigned int i = 0; i < numberRatesRatiosXSets; ++i){
                    sprintf( numberPos, "%d", i+1 );
                    if (!parameters.findParameter(string(label)+", prior")){
                        parameters[string(label)+", prior"] = "uniform(0,1.0)";
                    }
                    PerturbatorGaussParameter* pert =
                        perturbator->registerGaussParameter( *this, UpdateMessage::RATE,
                                   label, &((*xR)[i]), parameters, label,
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


void RatesRatios::initialiseMCMC( ParametersSet & parameters, Perturbator* perturbator ){
    char label[70];

    double extrMin, extrMax;

    //prepare the update message (base class) sent to the model for an update
    setType( UpdateMessage::PARAM_TYPE );
    setFlag( UpdateMessage::RATE );

    if (numberRatesRatiosSets != 0){
        registerVariationParameter( parameters, perturbator );
        //register rates
        if (!parameters.findParameter("Rate ratios, prior")){
            parameters["Rate ratios, prior"] = "uniform(.0005,2000.0)";
        }
        char* numberPos = label;
        char* numberPos2 = label;
        if (numberRatesRatiosSets!=1){
            sprintf ( label, "Rate ratios, Set ");
            while (*(++numberPos));
        }
        else{
            sprintf ( label, "Rate ratios no ");
            while (*(++numberPos2));
        }
        unsigned int increment = (gammaCatAdjust & L) ? numberRatesRatiosSets - 1 : 1;
        for ( unsigned int i = 0; i < numberRatesRatiosSets; i += increment){
            if (numberRatesRatiosSets!=1){
                sprintf( numberPos, "%d, no ", i+1 );
                numberPos2 = numberPos;
                while (*(++numberPos2));
            }
            for ( unsigned int j = 0; j < numberRatesRatios; ++j ){
                sprintf( numberPos2, "%d", j+1 );
                PerturbatorGaussParameter* pert =
                    perturbator->registerGaussParameter( *this, UpdateMessage::RATE, label,
                           &(ratesRatios[i][j]), parameters, "Rate ratios",
                           .2, .005, 10.0, "initial step" );
                pert->getExtrema( extrMin, extrMax );
                if ( extrMin<0.0 ){
                    cerr << "invalid prior for the rate ratios: "
                         << parameters["Rate ratios, prior"]
                         << ". Should not allow negative values" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
        //register X rates
        if (gammaCatAdjust & LX){
            sprintf ( label, "Rates Ratios, Set X");
            char* numberPos = label;
            while (*(++numberPos));
            for ( unsigned int i = 0; i < numberRatesRatiosXSets; ++i){
                sprintf( numberPos, "%d, no ", i+1 );
                numberPos2 = numberPos;
                while (*(++numberPos2));
                for ( unsigned int j = 0; j < numberRatesRatios; ++j ){
                    sprintf( numberPos2, "%d", j+1 );
                    PerturbatorGaussParameter* pert =
                        perturbator->registerGaussParameter( *this, UpdateMessage::RATE, label,
                               &(ratesRatios[i][j]), parameters, "Rate ratios",
                               .2, .005, 10.0, "initial step" );
                    pert->getExtrema( extrMin, extrMax );
                    if ( extrMin<0.0 ){
                        cerr << "invalid prior for the rate ratios: "
                             << parameters["Rate ratios, prior"]
                             << ". Should not allow negative values" << endl;
                        exit(EXIT_FAILURE);
                    }
                }
            }
        }
    } // end if numberRatesPrio != 0
    //more than one set? smoothing prior....
    if(numberRatesRatiosSets>1){
    }
}


void RatesRatios::getAllParameters( vector < double > & pars ) const {
    unsigned int increment;
    int index = -1;

    pars.resize( getNumberFreeParameters() );
    if ( (gammaCatAdjust & LS) && (!(gammaCatAdjust & LINK)||(gammaCatAdjust & LINKY)) ){
        pars[++index] = sqrt( *sR / ( 1.0 - *sR ) );
    }
    if ( (gammaCatAdjust & LX) && (!(gammaCatAdjust & LINK)||(gammaCatAdjust & LINKY)) ){
       //the variables optimized are x bounded in [0..1]
        for (unsigned int i = 0; i < numberRatesRatiosXSets; ++i){
            pars[++index] = sqrt( (*xR)[i] / ( 1.0 - (*xR)[i] ) );
        }
    }
    //we use a trick to limit the code when linear frequencies are used
    increment = (gammaCatAdjust & L) ? numberRatesRatiosSets - 1 : 1;
    for ( unsigned int rateMat = 0; rateMat < numberRatesRatiosSets; rateMat += increment) {
        for ( unsigned int ratio = 0; ratio < numberRatesRatios; ++ratio ) {
#ifdef DEBUG3
            assert( !isinf( ratesRatios[rateMat] [ratio] ) );
            assert( !isnan( ratesRatios[rateMat] [ratio] ) );
            assert( ratesRatios[rateMat] [ratio] >= 0.0 );
#endif
            pars[++index] = sqrt( ratesRatios[rateMat] [ratio] );
        }
    }
    if (gammaCatAdjust & LX){
        for ( unsigned int rateMat = 0; rateMat < numberRatesRatiosXSets; ++rateMat ) {
            for ( unsigned int ratio = 0; ratio < numberRatesRatios; ++ratio ) {
#ifdef DEBUG3
                assert( !isinf( ratesRatiosX[rateMat] [ratio] ) );
                assert( !isnan( ratesRatiosX[rateMat] [ratio] ) );
                assert( ratesRatiosX[rateMat] [ratio] >= 0.0 );
#endif
                pars[++index] = sqrt( ratesRatiosX[rateMat][ratio] );
            }
        }
    }
    assert(index+1==(int)pars.size());
}


void RatesRatios::setAllParameters(const vector < double > & pars ){
    int index = 0;

    if ( (gammaCatAdjust & LS) && (!(gammaCatAdjust & LINK)||(gammaCatAdjust & LINKY)) ){
        *sR = SQR(pars[index]) / ( 1.0 + SQR( pars[index] ) );
        ++index;
    }
    if ( (gammaCatAdjust & LX) && (!(gammaCatAdjust & LINK)||(gammaCatAdjust & LINKY)) ){
        for (unsigned int i = 0; i < numberRatesRatiosXSets; ++i){
            (*xR)[i] = SQR( pars[index] ) / ( 1.0 + SQR(pars[index]) );
            ++index;
        }
    }

    //we use a trick to limit the code when linear frequencies are used
    unsigned int increment = (gammaCatAdjust & L) ? numberRatesRatiosSets - 1 : 1;
    for ( unsigned int rateMat = 0; rateMat < numberRatesRatiosSets; rateMat += increment ){
        for ( unsigned int ratio = 0; ratio < numberRatesRatios; ++ratio ) {
            ratesRatios[rateMat] [ratio] = ( pars[index] * pars[index] );
#ifdef DEBUG3
            assert( !isnan( ratesRatios[rateMat] [ratio] ) );
            assert( !isinf( ratesRatios[rateMat] [ratio] ) );
#endif
            ++index;
        }
    }
    if (gammaCatAdjust & LX){
        for ( unsigned int rateMat = 0; rateMat < numberRatesRatiosXSets; ++rateMat ){
            for ( unsigned int ratio = 0; ratio < numberRatesRatios; ++ratio ) {
                ratesRatiosX[rateMat] [ratio] = ( pars[index] * pars[index] );
#ifdef DEBUG3
                assert( !isnan( ratesRatiosX[rateMat] [ratio] ) );
                assert( !isinf( ratesRatiosX[rateMat] [ratio] ) );
#endif
                ++index;
            }
        }
    }
    assert( index == (int)pars.size() );
}


unsigned int RatesRatios::getNumberFreeParameters() const{
    unsigned int numberFreeParameters = numberRatesRatios *
            ( ((gammaCatAdjust & LX) ? numberRatesRatiosXSets : 0) +
              ((gammaCatAdjust & L) ? 2 : numberRatesRatiosSets) );
    if ( (gammaCatAdjust & L) &&
         (!(gammaCatAdjust&LINK) || (gammaCatAdjust&LINKY)) ){
        if ( gammaCatAdjust & LS ){
            numberFreeParameters+=1;
        }
        if ( gammaCatAdjust & LX ){
            numberFreeParameters+=numberRatesRatiosXSets;
        }
    }
    return numberFreeParameters;
}

void RatesRatios::printLine( ostream & outputStream ){
    outputStream << setprecision( 4 );
    if ( (gammaCatAdjust & L) &&
         (!(gammaCatAdjust&LINK) || (gammaCatAdjust&LINKY)) ){
        if ( gammaCatAdjust & LS ){
            outputStream << setw( 9 ) << *sR << ' ';
        }
        if ( gammaCatAdjust & LX ){
            vector< double > x = (*xR);
            for( unsigned int i = numberRatesRatiosXSets-1; i >0; --i ){
                x[i-1] = x[i-1]*x[i];
            }
            for( unsigned int i = 0; i < numberRatesRatiosXSets; ++i ){
                outputStream << setw( 9 ) << x[i] << ' ';
            }
        }
    }
    unsigned int increment = (gammaCatAdjust & L) ? numberRatesRatiosSets - 1 : 1;
    for ( unsigned int rateMat = 0 ; rateMat < numberRatesRatiosSets; rateMat += increment ) {
        for ( unsigned int rate = 0 ; rate < numberRatesRatios; ++rate ) {
            outputStream << setw( 9 ) << ratesRatios[rateMat] [rate] << ' ';
        }
    }
    if ( gammaCatAdjust & LX ){
        for ( unsigned int rateMat = 0 ; rateMat < numberRatesRatiosXSets; ++rateMat ) {
            for ( unsigned int rate = 0 ; rate < numberRatesRatios; ++rate ) {
                outputStream << setw( 9 ) << ratesRatiosX[rateMat] [rate] << ' ';
            }
        }
    }
}

void RatesRatios::fromLine( const vector<double>& parameters ){
    int index = -1;
    if ( (gammaCatAdjust & L) &&
         (!(gammaCatAdjust&LINK) || (gammaCatAdjust&LINKY)) ){
        if ( gammaCatAdjust & LS ){
            *sR = parameters[++index];
        }
        if ( gammaCatAdjust & LX ){
            for( unsigned int i = 0; i < numberRatesRatiosXSets; ++i ){
                (*xR)[i] = parameters[++index];
            }
            for( unsigned int i = 0; i < numberRatesRatiosXSets-1; ++i ){
                (*xR)[i] = (*xR)[i] / (*xR)[i+1];
            }
        }
    }
    unsigned int increment = (gammaCatAdjust & L) ? numberRatesRatiosSets - 1 : 1;
    for ( unsigned int rateMat = 0 ; rateMat < numberRatesRatiosSets; rateMat += increment ) {
        for ( unsigned int rate = 0 ; rate < numberRatesRatios; ++rate ) {
            ratesRatios[rateMat][rate]  = parameters[++index];
        }
    }
    if ( gammaCatAdjust & LX ){
        for ( unsigned int rateMat = 0 ; rateMat < numberRatesRatiosXSets; ++rateMat ) {
            for ( unsigned int rate = 0 ; rate < numberRatesRatios; ++rate ) {
                ratesRatiosX[rateMat][rate]  = parameters[++index];
            }
        }
    }
    assert(index == (int)parameters.size()-1);
}

unsigned int RatesRatios::getNumberLineParameters() const{
    unsigned int numberLineParameters = numberRatesRatios *
            ( ((gammaCatAdjust & LX) ? numberRatesRatiosXSets : 0) +
              ((gammaCatAdjust & L) ? 2 : numberRatesRatiosSets) );
    if ( (gammaCatAdjust & L) &&
         (!(gammaCatAdjust&LINK) || (gammaCatAdjust&LINKY)) ){
        if ( gammaCatAdjust & LS ){
            ++numberLineParameters;
        }
        if ( gammaCatAdjust & LX ){
            numberLineParameters+=numberRatesRatiosXSets;
        }
    }
    return numberLineParameters;
}

void RatesRatios::getAllPriorParameters( vector<double>& params ) const{
    params.resize(getNumberPriorParameters());
}

void RatesRatios::setAllPriorParameters( const vector<double>& ){
    //assert(params.size()==0);
}

unsigned int RatesRatios::getNumberPriorParameters() const{
    return 0;
}


double RatesRatios::getLnPrior() const{
    return lnPrior;
}

double RatesRatios::computeLnPrior(){
    double ret = 0.0;
    if (priorType==VAR){
        assert(0.0);
    }
    return ret;
}

void RatesRatios::update(UpdateMessage* subject){
    assert(subject->hasType(UpdateMessage::PARAM_TYPE));
    if (subject->hasFlag(UpdateMessage::HYPER_PARAM)){
#ifdef DEBUG2
        cout << "ratesratios received a hyperparam message" << endl;
#endif
        assert(0);
        if (!subject->hasFlag(UpdateMessage::PRIOR_FLAG)){
            setFlag( PRIOR_FLAG );
            //computeCovMatrix();
            assert(0);
            lnPrior = computeLnPrior();
            notify();
            unsetFlag( PRIOR_FLAG );
        }
    }
    else{
#ifdef DEBUG2
        cout << "ratesratios received a modif update message" << endl;
#endif
        if (gammaCatAdjust & L) updateLinear();
        lnPrior = computeLnPrior();
        notify();
    }
}

