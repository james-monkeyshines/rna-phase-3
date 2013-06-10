#ifndef MCMCSUMMARIZE_H
#define MCMCSUMMARIZE_H

#include <fstream>
#include <iostream>
#include <string>

#include "Phase.h"

class ConsensusTree;
class Model;

using namespace std;

class MCMCsummarize : public Phase{
private:
    ConsensusTree* consensusTree;
    Model* consensusModel1;
    Model* consensusModel2;
    int expectedNumberSamples;

public:
    MCMCsummarize();

    int run( int argc, char* argv[] );
    
    int modelParametersConsensus( ifstream& modelParametersFile,
                            Model* consensusModel1, Model* consensusModel2 );
};

#endif //MCMCSUMMARIZE_H




