#include "Likelihood.h"


#include <iostream>
#include <iomanip>
#include <algorithm>
#include <ctype.h>
#include <functional>


#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"
#include "Util/randombox.h"

#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include "Models/Model.h"
#include "Tree/LikelihoodTree.h"

using namespace std;



void Likelihood::initObject(){
    // The tree returns the ancestral states and rates in its own format.
    // Sites are grouped by categories and invariants are lumped together
    // at the end of the sequence, we have to translate the tree
    // format into the 'initial' format to fit with the sequences
    tree->getRateCategory( rateCat );
    tree->getAncestralStates( ancestralStates );
    tree->getSiteSpecificLikelihood( siteSpecificLik );

    initialSequenceLength = table->getInitialLength();

    numberCategories = table->getNumberCategories();
    assert(rateCat.size()==numberCategories);
    assert(ancestralStates.front().second.size()==numberCategories);

    length.resize(numberCategories);
    invLength.resize(numberCategories);
    nonInvLength.resize(numberCategories);
    numberStates.resize(numberCategories);
    numberRates.resize(numberCategories);

    for ( unsigned int symbolCat = 0; symbolCat < numberCategories; ++symbolCat ){
        invLength[symbolCat] = table->getInvariantsLength( symbolCat );
        nonInvLength[symbolCat] = table->getSequencesLength( symbolCat );
        length[symbolCat] = invLength[symbolCat] + nonInvLength[symbolCat];
        assert(ancestralStates.front().second[symbolCat].numberRows() == length[symbolCat]);
        assert(rateCat[symbolCat].numberRows() == length[symbolCat]);

        numberStates[symbolCat] = model->getNumberStates( symbolCat );
        numberRates[symbolCat] = model->getNumberRatesCategories( symbolCat );
        assert(numberStates[symbolCat]==ancestralStates.front().second[symbolCat].numberColumns());
        assert(numberRates[symbolCat]==rateCat[symbolCat].numberColumns());
    }
}


double Likelihood::getStateFrequencies( const vector< bool > & speciesMask,
                                        unsigned int symbolCategory,
                                        const vector< bool >& rateCategoriesMask,
                                        vector<double>& stateFreq){

    stateFreq.resize(numberStates[symbolCategory]);
    fill(stateFreq.begin(), stateFreq.end(), 0.0);

    //the sequence table is responsible for the reordering of the sequence
    //by category. table->retrieveInitialNucleotides(category, siteIndex)
    //will return the initial nucleotides position.

    //load the contemporary sequences
    const array2D< string > & seq = table->getSequences( symbolCategory );
    const vector< pair< string, unsigned int > > inv =
            table->getInvariantBases( symbolCategory );

    //prepare the vector to store the composition in contemporary sequence
    //at this site
    vector< double > freqSite;
    freqSite.resize( numberStates[symbolCategory] );


    double postProbRateMask = 0.0;
    double total = 0.0;

    unsigned int countSpecies = count( speciesMask.begin(), speciesMask.end() , true );
    for ( unsigned int site = 0; site < length[symbolCategory]; ++site ){
        //we need to know what this site correspond to in the initial
        //sequence: fill initpos
        //initPos[] is a vector which contains 1 (or more when it was
        //an invariant) set of nucleotides associated with the state.
        //factor is equal to the actual number of symbols represented
        //by this site in the sequence used by the tree
        //(invariant sites are lumped into one site)
        vector < vector< unsigned int > > initPos;
        vector < unsigned int > initRaw;
        unsigned int factor;
        if ( site < nonInvLength[symbolCategory] ){
            initPos = table->retrieveInitialNucleotides( symbolCategory, site );
            factor = 1;
        }
        else{ //invariant site
            int invId = site-nonInvLength[symbolCategory];
            initPos = table->retrieveInitialNucleotides( symbolCategory, -1 - invId );
            factor = inv[invId].second;
        }

        //now fill the rateComposition vector, we have to retrieve the
        //contemporary composition at this site and add that to the
        //current values
        fill (freqSite.begin(), freqSite.end(), 0);
        //for variable sites
        if ( site < nonInvLength[symbolCategory] ){
            for (unsigned int species = 0; species < seq.numberRows(); ++species){
                if (speciesMask[species]){
                    unsigned int nb = 0;
                    for (unsigned int state = 0; state < numberStates[symbolCategory]; ++state){
                        if ( model->getEquivalencyTable(symbolCategory)(
                            model->getSymbolNumber( seq(species, site), symbolCategory ),
                                state ) == 1 ){
                            ++nb;
                        }
                    }
                    for (unsigned int state = 0; state < numberStates[symbolCategory]; ++state){
                        if ( model->getEquivalencyTable(symbolCategory)(
                             model->getSymbolNumber( seq(species, site), symbolCategory ),
                             state ) == 1 ){
                            freqSite[state] = freqSite[state] + 1.0/(double)nb;
                        }
                    }
                }
            }
            for (unsigned int state = 0; state < numberStates[symbolCategory]; ++state){
                freqSite[state] /= (double)countSpecies;
            }
        }
        //for invariant sites it is actually easier
        else{
            unsigned int nb = 0;
            for (unsigned int state = 0; state < numberStates[symbolCategory]; ++state){
                if ( model->getEquivalencyTable(symbolCategory)(
                        model->getSymbolNumber( inv[site-nonInvLength[symbolCategory]].first, symbolCategory ),
                            state ) == 1 ){
                    ++nb;
                }
            }
            for (unsigned int state = 0; state < numberStates[symbolCategory]; ++state){
                if ( model->getEquivalencyTable(symbolCategory)(
                     model->getSymbolNumber( inv[site-nonInvLength[symbolCategory]].first, symbolCategory ),
                     state ) == 1 ){
                    freqSite[state] = 1.0/(double)nb;
                }
            }
        }

        //and add this site to the frequencies
        for( unsigned int rate=0; rate < numberRates[symbolCategory];  ++rate ){
            total += rateCat[symbolCategory](site, rate) * (double)factor;
            //if in the rate mask
            if (rateCategoriesMask[rate]){
                postProbRateMask += rateCat[symbolCategory](site, rate) * (double)factor;
                for (unsigned int state = 0; state < numberStates[symbolCategory]; ++state){
                    stateFreq[state] +=
                            rateCat[symbolCategory](site, rate) * freqSite[state] * (double)factor;
                }
            }
        }
    } // for all site in the category
    unsigned int invNumber = 0;
    for ( unsigned int j = 0; j < invLength[symbolCategory]; ++j ){
        invNumber += inv[j].second;
    }
    assert(round(total)==nonInvLength[symbolCategory] + invNumber);
    total = 0.0;
    for (unsigned int state = 0; state < numberStates[symbolCategory]; ++state){
        total += stateFreq[state];
        stateFreq[state] /= postProbRateMask;
    }
    return (postProbRateMask / (double)(nonInvLength[symbolCategory] + invNumber) );
}


void Likelihood::getMaxProb( const array2D<double> & probMatrix, vector<unsigned int> & maxProb ){
    maxProb.resize(probMatrix.numberRows());
    for( unsigned int i = 0; i < probMatrix.numberRows(); ++i ){
        maxProb[i]=0;
        double best = probMatrix(i,0);
        for ( unsigned int j = 1; j < probMatrix.numberColumns(); ++j ){
            if ( probMatrix(i,j) > best){
                maxProb[i]=j;
                best = probMatrix(i,j);
            }
        }
    }
}


void Likelihood::getAncestralSequences(vector< vector< vector<double> > > & probReconstruction, vector< pair<InferenceNode*, string> > & ancestralSequences){
    //prepare the ancestral string and fill it with gaps
    unsigned int numberInternalNodes = ancestralStates.size();
    ancestralSequences.resize(numberInternalNodes);
    probReconstruction.resize(numberInternalNodes);
    assert(ancestralStates.size()==numberInternalNodes);
    list< pair<InferenceNode*, vector< array2D<double> > > >::iterator iter = ancestralStates.begin();
    for ( unsigned int i = 0; i < numberInternalNodes; ++i ){
        ancestralSequences[i].first = iter->first;
        ++iter;
        ancestralSequences[i].second.resize(initialSequenceLength);
        fill(ancestralSequences[i].second.begin(), ancestralSequences[i].second.end(), '-');
        probReconstruction[i].resize(initialSequenceLength);
    }
    for ( unsigned int symbolCat = 0; symbolCat < numberCategories; ++symbolCat ){
        for ( unsigned int site = 0; site < length[symbolCat]; ++site ){
            //we need to know what this site correspond to in the initial
            //sequence: fill initpos
            //initPos[] is a vector which contains 1 (or more when it was
            //an invariant) set of nucleotides associated with the state.
            //factor is equal to the actual number of symbols represented
            //by this site in the sequence used by the tree
            //(invariant sites are lumped into one site)
            vector < vector< unsigned int > > initPos;
            if ( site < nonInvLength[symbolCat] ){
                initPos = table->retrieveInitialNucleotides( symbolCat, site );
            }
            else{ //invariant site
                int invId = site-nonInvLength[symbolCat];
                initPos = table->retrieveInitialNucleotides( symbolCat, -1 - invId );
            }

            //initPos can be common to all internal nodes but what follows cannot
            vector<unsigned int> bestPostProbState;
            string stateSymbol;
            list< pair<InferenceNode*, vector< array2D<double> > > >::iterator iter = ancestralStates.begin();
            for ( unsigned int i = 0; i < numberInternalNodes; ++i ){
                getMaxProb( iter->second[symbolCat], bestPostProbState );
                /*bestPostProbState is the index of the state for the model
                  map this value to the symbol */
                stateSymbol = model->getState(bestPostProbState[site],symbolCat);
                array2D<double> & probStatesMatrix = iter->second[symbolCat];
                /**** FILL ancestralSequence and probReconstruction with this site */
                //for all the sets in initPos we also have to write at the
                //corresponding position in the string the nucleotides
                for ( unsigned int repeat = 0; repeat < initPos.size(); ++repeat ){
                    assert( stateSymbol.size() == initPos[repeat].size() );
                    //we assume that the symbol used for the state correspond
                    //to the nucleotides... state name MUST be of the correct size
                    //mismatch pairs will be assigned the two nucleotides M and M
                    //with & state models
                    for (unsigned int nuc = 0; nuc < initPos[repeat].size(); ++nuc){
                        ancestralSequences[i].second[initPos[repeat][nuc]] = stateSymbol[nuc];
                        vector<double> & probStates = probReconstruction[i][initPos[repeat][nuc]];
                        probStates.resize(numberStates[symbolCat]);
                        for ( unsigned int state = 0; state < numberStates[symbolCat]; ++state ){
                            probStates[state] = probStatesMatrix(site, state);
                        }
                    }
                }
                ++iter;
           } //for all internal node
        } // for all site in the category
    } //for all category
    for ( unsigned int i = 0; i < numberInternalNodes; ++i ){
        if ( find( ancestralSequences[i].second.begin(), ancestralSequences[i].second.end(), '-' ) !=
                      ancestralSequences[i].second.end() ){
            cerr << "Sorry, error during the ancestral sequences reconstruction" << endl;
            exit(EXIT_FAILURE);
        }
    }
}



void Likelihood::getRateCategory( vector< vector<double> > & probReconstruction, vector<unsigned int>& ratesMAP){
    //prepare the result and fill it with gaps
    ratesMAP.resize(initialSequenceLength);
    probReconstruction.resize(initialSequenceLength);
    fill( ratesMAP.begin(), ratesMAP.end(), (unsigned int) -1 );

    for ( unsigned int symbolCat = 0; symbolCat < numberCategories; ++symbolCat ){
        vector<unsigned int> bestPostProbRate;
        getMaxProb( rateCat[symbolCat], bestPostProbRate );
        for ( unsigned int site = 0; site < length[symbolCat]; ++site ){
            //we need to know what this site correspond to in the initial
            //sequence: fill initpos
            //initPos[] is a vector which contains 1 (or more when it was
            //an invariant) set of nucleotides associated with the state.
            //factor is equal to the actual number of symbols represented
            //by this site in the sequence used by the tree
            //(invariant sites are lumped into one site)
            vector < vector< unsigned int > > initPos;
            if ( site < nonInvLength[symbolCat] ){
                initPos = table->retrieveInitialNucleotides( symbolCat, site );
            }
            else{ //invariant site
                int invId = site-nonInvLength[symbolCat];
                initPos = table->retrieveInitialNucleotides( symbolCat, -1 - invId );
            }

            /**** FILL rate with this site */
            //for all the sets in initPos we also have to write at the
            //corresponding position in the rate vector
            for (unsigned int repeat = 0; repeat < initPos.size(); ++repeat){
                for (unsigned int nuc = 0; nuc < initPos[repeat].size(); ++nuc){
                    ratesMAP[initPos[repeat][nuc]] = bestPostProbRate[site];
                    vector<double> & probRates = probReconstruction[initPos[repeat][nuc]];
                    probRates.resize( numberRates[symbolCat] );
                    for (unsigned int rate = 0; rate < numberRates[symbolCat]; ++rate){
                        probRates[rate] = rateCat[symbolCat](site, rate);
                    }
                }
            }
        } // for all site in the category
    } //for all category
    if ( find( ratesMAP.begin(), ratesMAP.end(), (unsigned int) -1 ) != ratesMAP.end() ){
        cerr << "Sorry, error during the evaluation of posterior probablities of the site-specific evolutionary rate" << endl;
        exit(EXIT_FAILURE);
    }
}

int Likelihood::run( int argc, char* argv[]){

    Phase::run( argc, argv );
    
    if ( argc != 2 ) {
        cerr << "usage : " << argv[0] << " control_file" << endl;
        exit(EXIT_FAILURE);
    }

    // Create a FileParser from the control-file
    FileParser * fileParser = new FileParser( argv[1] );

    // Retrieve the control-file parameters
    ParametersSet* ret = fileParser->retrieveParametersSet();
    delete fileParser;
    if (!ret){
        cerr << "invalid control file " << argv[1] << endl;
        exit(EXIT_FAILURE);
    }
    ParametersSet parameters = *ret;
    delete ret;


    /** ************************************************************************
     * Data File Initialization
     ************************************************************************ */

    // Read in molecular data and remove unknown characters from the data matrix
    string dataFile = parameters( "DATAFILE" ).stringParameter( "Data file" );
    table = new SequenceTable( parameters( "DATAFILE" ) );

    /** ************************************************************************
     * Model Initialization
     ************************************************************************ */
    Singleton < Factory<Model> > & modelFactory =
            Singleton < Factory<Model> >::instance();
    model = modelFactory.create( parameters( "MODEL" ).stringParameter( "Model" ),
         parameters( "MODEL" ) );
    model->initialisation( NULL );
    // Create a FileParser from the model parameters-file
    FileParser * modelFileParser = new FileParser( parameters("MODEL").stringParameter("Model parameters file") );
    // Retrieve the model parameters
    ParametersSet* modelParameters = modelFileParser->retrieveParametersSet();
    delete modelFileParser;
    if (!modelParameters){
        cerr << "Invalid model parameters file" << endl;
        exit(EXIT_FAILURE);
    }
    //initialise the model with them
    model->setModelParameters(*modelParameters);
    model->validChange();
    delete modelParameters;

    /** ************************************************************************
     * Tree Initialization
     ************************************************************************ */
    string userTreeFileName = parameters("TREE").stringParameter( "Tree file" );
    ifstream userTreeFile( userTreeFileName.c_str() );
    string stringTree;
    FileParser::readTree(userTreeFile, stringTree);
    if ( stringTree.length() == 0 ){
        cerr << "Error while reading the tree in the file : "
             << userTreeFileName << endl;
        exit(EXIT_FAILURE);
    }
    tree = new LikelihoodTree( stringTree );
    
    cout << tree->toString() << endl;

    tree->loadDataAndModel(table, model);

    model->printParameters( cout );

    cout.setf(ios::fixed);
    cout << setprecision(4);
    cout << "likelihood = " << tree->loglikelihood() << endl;

    initObject();

    if ( parameters.findParameter("Ancestral sequences") ){
        /* output the ancestral sequences */
        ofstream output;
        FileParser::confirmOpenFile( output,  parameters.stringParameter("Ancestral sequences"), true );
        vector< vector< vector<double> > > probReconstruction;
        vector< pair<InferenceNode*, string> > ancestralSequences;
        getAncestralSequences( probReconstruction, ancestralSequences );
        output << "Ancestral sequences..." << ancestralSequences[0].second.size()
             << " nucleotides" << endl;
        for ( unsigned int node = 0; node < ancestralSequences.size(); ++node ){
            output << '#' << ancestralSequences[node].first->toString(true) << endl;
            output << ancestralSequences[node].second << endl;
            /* output BPP for each nucleotide */
            output << "#Posterior probabilities" << endl;
            output.setf(ios::fixed);
            output << setprecision(5);
            for ( unsigned int i = 0; i < probReconstruction[node].size(); ++i ){
                output << setw(10) << i+1 << "   ";
                for ( unsigned int j = 0; j < probReconstruction[node][i].size(); ++j ){
                    output << setw(10) << probReconstruction[node][i][j] << ' ';
                }
                output << endl;
            }
        }
    }
    
    if ( parameters.findParameter("Site-specific loglikelihood") ){
        //output the rate vector
        ofstream output;
        FileParser::confirmOpenFile( output,  parameters.stringParameter("Site-specific loglikelihood"), true );
        output.setf(ios::fixed);
        output << setprecision(5);
        for ( unsigned int i = 0; i < siteSpecificLik.size(); ++i ){
            output << "#site-specific loglikelihoods, category/model " << i+1 << " (warning, order is meaningless)" << endl;
            for ( unsigned int site = 0; site < siteSpecificLik[i].size(); ++site ){
                output << setw(9) << siteSpecificLik[i][site] << ' ';
            }
            output << endl << endl;
        }
    }

    
    if ( parameters.findParameter("Site-specific substitution rates") ){
        //output the rate vector
        ofstream output;
        FileParser::confirmOpenFile( output,  parameters.stringParameter("Site-specific substitution rates"), true );
        vector< vector<double> > probReconstruction;
        vector<unsigned int> outputRate;
        getRateCategory( probReconstruction, outputRate );
        output << "#Maximum a posteriori of the discrete rate categories (and/or invariant) at each site" << endl;
        for (unsigned int i = 0; i < outputRate.size(); ++i){
            output << outputRate[i] + 1 << ' ';
        }
        output << endl << endl;
        /* output BPP for each rate */
        output << "#Bayesian posterior probabilities" << endl;
        output.setf(ios::fixed);
        output << setprecision(5);
        for ( unsigned int i = 0; i < probReconstruction.size(); ++i ){
            output << setw(10) << i+1 << "   ";
            for ( unsigned int j = 0; j < probReconstruction[i].size(); ++j ){
                output << setw(10) << probReconstruction[i][j] << ' ';
            }
            output << endl;
        }        
    }

    /*now compute composition for the given rate categories for the chosen set of species*/
    unsigned int count = 0;
    if ( parameters.findParameter("Rate vs. composition") &&
        parameters.boolParameter("Rate vs. composition") ){
        vector< bool > specMask;
        specMask.resize(table->getNumberSpecies());
        if ( !parameters.findParameter("Selected species") ||
              parameters.stringParameter("Selected species") == "all" ){
            fill( specMask.begin(), specMask.end(), true);
            count = specMask.size();
        }
        else{
            string name = parameters.stringParameter("Selected species");
            /* if this is a species name set the mask accordingly */
            int speciesId = table->findSequence(name);
            if (speciesId!=-1){
                specMask[speciesId] = true;
                count = 1;
            }
            /* otherwise it must be a filename */
            else{
                ifstream speciesFile(name.c_str());
                speciesFile >> ws;
                if ( !speciesFile.good() || speciesFile.eof() ){
                    cerr << "Sorry, invalid file " << name << ". Cannot read your set of species." << endl;
                    exit(1);
                }
                while ( !speciesFile.eof() ){
                    string species;
                    speciesFile >> species;
                    int speciesId = table->findSequence(species);
                    if ( speciesId == -1 ){
                        cerr << "Sorry, unknown species name (" << species << ") in your species file." << endl;
                        exit(1);
                    }
                    specMask[speciesId] = true;
                    speciesFile >> ws;
                }
            }
        }
        
        vector< bool > rateMask;
        vector<double> stateFreq;
        for (unsigned int cat = 0; cat < model->getNumberSymbolCategory(); ++cat){
            rateMask.resize(numberRates[cat]);
            cout << "Estimated frequencies for each rate categories for the selected set of " << count << " species." << endl;
            cout << setprecision(5);
            double tot;
            for (unsigned int i = 0; i < model->getNumberStates( cat ); ++i){
                cout << setw(8) << model->getState( i, cat ) << setw(2) << ' ';
            }
            cout << setw(21) << "category BPP" << endl;
            for (unsigned int rate = 0; rate < rateMask.size(); ++rate){
                fill( rateMask.begin(), rateMask.end(), false);
                rateMask[rate] = true;
                tot = getStateFrequencies( specMask, cat, rateMask, stateFreq);
                for (unsigned int i = 0; i < stateFreq.size(); ++i){
                    cout << setw(10) << stateFreq[i];
                }
                cout << setw(20) << tot << endl;
            }
            cout << "Mean values" << endl;
            fill( rateMask.begin(), rateMask.end(), true );
            tot = getStateFrequencies( specMask, cat, rateMask, stateFreq);
            for (unsigned int i = 0; i < stateFreq.size(); ++i){
                cout << setw(10) << stateFreq[i];
            }
            cout << setw(20) << tot << endl;
        }
    }
    parameters.checkAllUsed();
    return (EXIT_SUCCESS);
}


Likelihood::~Likelihood(){
    delete table;
    delete model;
}



int main (int argc, char** argv){
    Likelihood likelihood;
    int res = likelihood.run( argc, argv );
    return res;
}
