
#include "Sequence/SequenceTable.h"
#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include <vector>
#include <iostream>
#include <algorithm>


using namespace std;    
    
int main( int argc, char** argv ){
    //output stream
    ofstream *outputFileStreams;
    
    /** ************************************************************************
     *  Usage and explanations
     *
     ************************************************************************* */
    if (argc!=2){
        cerr << "Usage: " << argv[0] << " sequence_file" << endl << endl;
        cerr << "This program splits the nucleotides in sequence_file according to their" << endl
             << "associated model and writes back the sequences in the files sequence_file1," << endl
             << "sequence_file2, ... , sequence_fileX where X is the number of models found" << endl
             << "in your data file. cf DATAFILE section in the manual" << endl;
        exit(EXIT_FAILURE);
    }
    
    ParametersSet parameters("DATAFILE");
    string answer;
    parameters["Data file"] = argv[1];
    cout << "Interleaved data file?" << endl;
    cin >> answer;
    parameters["Interleaved data file"] = answer;
    if (parameters.boolParameter("Interleaved data file")){
        cout << "Interleaved structure?" << endl;
        cin >> answer;
        parameters["Interleaved structure"] = answer;
    }
    cout << "Heterogeneous data models?" << endl;
    cin >> answer;
    parameters["Heterogeneous data models"] = answer;
    SequenceTable sequenceTable( parameters );
    
    char fileName[150];
    sprintf( fileName, "%s.", argv[1] );
    char* numberPos = fileName;
    while(*numberPos) ++numberPos;
    outputFileStreams = new ofstream[sequenceTable.getNumberCategories()];
    
    for( unsigned int i = 0; i < sequenceTable.getNumberCategories(); ++i ){
        sprintf( numberPos, "%d", i + 1 );
        FileParser::confirmOpenFile(outputFileStreams[i], fileName, true);
    }
    const vector<int> & category = sequenceTable.getInitialCategoryAttrib();
    
    for(int i = 0; i < (int)sequenceTable.getNumberCategories(); ++i){
        outputFileStreams[i] << sequenceTable.getNumberSpecies() << "  ";
        outputFileStreams[i] << count(category.begin(), category.end(), i+1) << "   ";
        outputFileStreams[i] << "TYPE_UNKNOWN_PLEASE_REPLACE" << endl << endl;
    }
    const vector<char> & structure = sequenceTable.getStructure();
    for(unsigned int k = 0; k < structure.size(); ++k){
        outputFileStreams[category[k]-1] << structure[k];
    }
    
    for(unsigned int j = 0; j < sequenceTable.getNumberSpecies(); ++j){
        for (unsigned int k = 0; k < sequenceTable.getNumberCategories(); ++k ){
            outputFileStreams[k] << endl << endl << sequenceTable.species[j] << endl;
        }
        for(unsigned int k = 0; k < sequenceTable.getInitialLength(); ++k){
            outputFileStreams[category[k]-1] << sequenceTable.getInitialSequences( j, k);
        }
    }
    for(unsigned int i = 0; i < sequenceTable.getNumberCategories(); ++i){
        outputFileStreams[i].close();
    }
}
