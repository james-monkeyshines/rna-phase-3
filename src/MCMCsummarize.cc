#include "MCMCsummarize.h"

#include "Tree/ConsensusTree.h"

#include "Models/Model.h"
#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include <iostream>
#include <iomanip>

using namespace std;


MCMCsummarize::MCMCsummarize(){
    consensusTree = NULL;
    consensusModel1 = NULL;
    consensusModel2 = NULL;
}

int MCMCsummarize::run( int argc, char* argv[] ){
  
    Phase::run( argc, argv );

    if ( argc != 2 ) {
        cerr << "usage: " << argv[0] << " mcmc_control_file" << endl;
        exit(EXIT_FAILURE);
    }
    
    string mcmcControlFileName = string( argv[1] );
  
    FileParser fileParser( mcmcControlFileName );
    ParametersSet* ret = fileParser.retrieveParametersSet();
    if (!ret){
        cerr << "malformed control file " << mcmcControlFileName << endl;
        exit(EXIT_FAILURE);
    }
    ParametersSet controlMCMC = *ret;
    delete ret;

    expectedNumberSamples = controlMCMC.intParameter("Sampling iterations");
    expectedNumberSamples = expectedNumberSamples /
                    controlMCMC.intParameter("Sampling period");

    //create the factory
    Singleton < Factory<ConsensusTree> > & consensusTreeFactory =
          Singleton < Factory<ConsensusTree> >::instance();
    //retrieve the registration name of the corresponding consensus tree
    //(by convention replace "MCMC" with "consensus")
    string consensusName = controlMCMC("TREE").stringParameter( "Tree" );
    unsigned int pos = consensusName.find( "MCMC", 0 );
    if (pos==consensusName.length()){
        cerr << "the tree in your control file is not a MCMC tree" << endl;
        exit(EXIT_FAILURE);
    }
    consensusName.replace( pos, 4, "consensus" );
    //create the right consensus tree
    ConsensusTree* consensusTree =
            consensusTreeFactory.create( consensusName, controlMCMC );
    consensusTree->process();
    consensusTree->produceConsensus();
    consensusTree->writeConsensusFile();
    string basename = controlMCMC.stringParameter("Output file");
    string consensusTreeFileName = basename + string("-consensus.tre");
    ofstream consensusTreeFile( consensusTreeFileName.c_str() );

    Singleton < Factory<Model> > & modelFactory =
            Singleton < Factory<Model> >::instance();
    consensusModel1 = modelFactory.create( controlMCMC("MODEL").
                                                     stringParameter( "Model" ),
                                           controlMCMC( "MODEL" ) );
    consensusModel2 = modelFactory.create( controlMCMC("MODEL").
                                                     stringParameter( "Model" ),
                                           controlMCMC( "MODEL" ) );
    consensusModel1->initialisation( NULL );
    consensusModel2->initialisation( NULL );


    string modelParametersFileName = basename + string(".mp");
    ifstream modelParametersFile( modelParametersFileName.c_str() );

    string consensusModelFileName1 = basename + string("-consensus-1.mod");
    string consensusModelFileName2 = basename + string("-consensus-2.mod");

    ofstream consensusModelFile1;
    ofstream consensusModelFile2;
    FileParser::confirmOpenFile( consensusModelFile1, consensusModelFileName1, true );
    FileParser::confirmOpenFile( consensusModelFile2, consensusModelFileName2, true );

    modelParametersConsensus( modelParametersFile,
                              consensusModel1, consensusModel2 );

    ParametersSet consensusParameters( "consensusParameters" );

    cout << "Consensus 1 for the model" << endl;
    consensusModel1->printParameters( cout );
    consensusParameters = consensusModel1->getModelParameters();
    consensusParameters.saveToFile(consensusModelFile1);
    consensusModelFile1.close();
    cout << "Consensus 2 for the model" << endl;
    consensusModel2->printParameters( cout );
    consensusParameters = consensusModel2->getModelParameters();
    consensusParameters.saveToFile(consensusModelFile2);
    consensusModelFile1.close();

    cout << "Consensus for the tree" << endl;
    cout << consensusTree->toString( false );
    if (controlMCMC.stringParameter("Output format") != "bambe" ){
        cout << ';' << endl;
    }
    else{
        cout << endl;
    }
    consensusTreeFile << consensusTree->toString( false ) << endl;
    consensusTreeFile.close();

    return(0);
}

int MCMCsummarize::modelParametersConsensus( ifstream& modelParametersFile,
                              Model* modelConsensus1, Model* modelConsensus2 ){
    unsigned int numberParameters = modelConsensus1->getNumberLineParameters();
    vector<double> consensusParameters1( numberParameters );
    vector<double> consensusParameters2 = consensusParameters1;

    vector< map< double , int > > sortedValues(numberParameters);

    double value;
    int numberInstances = 0;

    modelParametersFile >> ws;
    while ( !modelParametersFile.eof() ){
        for ( unsigned int i = 0; i < numberParameters; ++i ) {
            modelParametersFile >> value;
            if ( modelParametersFile.eof() ) {
                cerr << "error while reading the \".mp\" file, sample no:"
                << numberInstances + 1 << ", param no :" << i + 1 << endl;
                exit(EXIT_FAILURE);
            }
            consensusParameters1[i] += value;
            ++((sortedValues[i])[value]);
        }
        ++numberInstances;
        if (numberInstances%10000==0){
            cout << numberInstances << endl;
        }
        modelParametersFile >> ws;
    }
    if (numberInstances != expectedNumberSamples){
        cerr << "Warning: the sample size (" << numberInstances
             << ") and the expected size (" << expectedNumberSamples
             << ") are different" << endl;
    }

    map< double, int >::iterator iter;

    // consensus for each parameter
    for ( unsigned int i = 0; i < numberParameters; ++i ) {
        //mean for the 1st consensus
        consensusParameters1[i] /= numberInstances;

        //median for the 2nd consensus
        double medianValue = ((double)numberInstances)/2.0;
        iter = sortedValues[i].begin();
        int sampleCounter = (*iter).second;
        while ( (double)sampleCounter < medianValue ){
            ++iter;
            sampleCounter += (*iter).second;
        }
        consensusParameters2[i] = (*iter).first;

        //max for the 2nd consensus
        //double maxFreq = 0.0;

/*
        //PDF of the parameter
        double minValue = (*(sortedValues[i].begin())).first;
        double maxValue = (*(--sortedValues[i].end())).first;
        double halfBinSize = (maxValue-minValue)*100.0/(double)numberInstances;
        double binMaxValue = minValue + (2 * halfBinSize);
        int counter = 0;
        for ( iter = sortedValues[i].begin();
              iter != sortedValues[i].end(); ++iter){
              //add the number of item with that value to the actual bin
              if( (*iter).first <= binMaxValue ){
                  counter += (*iter).second;
              }
              //or finish the bin and start the next one if we pass
              // the upper boundary
              else{
                  //output result for this bin
                  outputPDFFile[i] << setw(13) << binMaxValue-halfBinSize
                                   << ' ' << setw(13)
                                   << (double)counter/(double)numberInstances
                                   << endl;

                  if ( ((double)counter/(double)numberInstances) > maxFreq ){
                      consensusParameters2[i] = binMaxValue-halfBinSize;
                      maxFreq = (double)counter/(double)numberInstances;
                  }
                  //prepare the next bin
                  counter = (*iter).second;
                  binMaxValue += 2*halfBinSize;
              }
        }
        //last empty bin
        outputPDFFile[i] << setw(13) << binMaxValue-halfBinSize << ' '
                         << setw(13) << 0.0 << endl;
        outputPDFFile[i].close();
*/
    }

    modelConsensus1->fromLine( consensusParameters1 );
    modelConsensus1->validChange();
    modelConsensus2->fromLine( consensusParameters2 );
    modelConsensus2->validChange();

    return (EXIT_SUCCESS);
}



int main( int argc, char * * argv ) {
    MCMCsummarize summarize;
    return summarize.run( argc, argv );
}
