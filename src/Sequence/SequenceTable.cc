#include "Sequence/SequenceTable.h"

#include "Models/Model.h"
#include "Util/ParametersSet.h"

#include "Util/FileParser.h"

#include <algorithm>
#include <functional>
#include <stack>
#include <queue>
#include <set>
#include <map>
#include <iostream>
#include <cstdlib>

#include <assert.h>
#include <ctype.h>

using namespace std;

void SequenceTable::processStructure() {

    //create a stack to store the positions of '(', '[', '{', '<' in the pairs list
    stack<unsigned int> lifo[4];

    //create stacks to store the positions of uppercase letter in the pairs list
    // Starting with a upper case letter means conventional nested pairing (12345.....54321).
    map< char, stack<unsigned int> > lifoMap;

    //create queues to store the positions of lowercase letter in the pairs list
    // The other pairing rule is used when starting with a lower case letter (123.4567......12..34567)
    map< char, queue<unsigned int> > fifoMap;

    bool error = false;

    //read the mask 'structure' and build the 'pairs' (RNA) and 'triplets' (codons) list.
    for (unsigned int position = 0; (error==false) && (position < structure.size()); ++position ) {
        char structElement = structure[position];

        // if '1' is found, then there must be a codon ...
        if (structElement == '1') {
            // ... and the next elements must be '2' and '3'
            if ( (position+2 >= structure.size()) || (structure[position+1]!='2') || (structure[position+2]!='3') ){
                cerr << "illogical structure :" << endl;
                cerr << "To describe a codon, a '1' must be followed by a '2' and a '3'" << endl;
                exit(EXIT_FAILURE);
            }
            triplets.push_back(position);
            position+=2;
        }
        else{
            //if lowercase letter is found ....
            if ( ( structElement >= 'a' ) && ( structElement <= 'z' ) ){
                char upper = toupper(structElement);
                //without a previous uppercase push structElement in the fifo queue,
                //this is the non-standard pairing first 'a' with first 'A'
                if ( lifoMap.find( upper ) == lifoMap.end() ){
                    fifoMap[structElement].push( pairs.size() );
                    pairs.push_back( pair < int, int > ( position, 0 ) );
                }
                //otherwise pop the lifo stack of the corresponding uppercase letter
                //this is the standard case last 'A' with first 'a' and we have the pair
                else{
                    if ( lifoMap[upper].empty() ) {
                        error = true;
                    }
                    else{
                        pairs[lifoMap[upper].top()].second = position;
                        lifoMap[upper].pop();
                    }
                }
            }
            else{
                if ( ( structElement >= 'A' ) && ( structElement <= 'Z' ) ){
                    char lower = tolower(structElement);
                    //without a previous lowercase push structElement in the lifo stack
                    //this is the standard case last 'A' with first 'a' and we have the first part
                    if ( fifoMap.find( lower ) == fifoMap.end() ){
                        lifoMap[structElement].push( pairs.size() );
                        pairs.push_back( pair < int, int > ( position, 0 ) );
                    }
                    //otherwise pop the fifo queue of the corresponding lowercase letter
                    //this is the non-standard pairing first 'a' with first 'A' and we have the pair
                    else{
                        if ( fifoMap[lower].empty() ) {
                            error = true;
                        }
                        else{
                            pairs[fifoMap[lower].front()].second = position;
                            fifoMap[lower].pop();
                            //issue a strong warning, this pairing should not be used
                            cerr << "WARNING: you seem to be using the non-standard pairing lowercase then uppercase" << endl
                                 << "Are you sure you know what you are doing? (consult section structure mask in the manual)" << endl;
                        }
                    }
                }
                else{
                    bool opening = false;
                    int lifoId = -1;
                    switch (structElement){
                        case '(': case ')':
                            lifoId = 0;
                            opening = (structElement == '(');
                        break;
                        case '[': case ']':
                            lifoId = 1;
                            opening = (structElement == '[');
                        break;
                        case '{': case '}':
                            lifoId = 2;
                            opening = (structElement == '{');
                        break;
                        case '<': case '>':
                            lifoId = 3;
                            opening = (structElement == '<');
                        break;
                        case '.': case '-':
                        break;
                        default:
                            error = true;
                    }
                    if (lifoId!=-1){
                        if ( opening ){
                            lifo[lifoId].push( pairs.size() );
                            pairs.push_back( pair < int, int > ( position, 0 )     );
                        }
                        else{
                            if (lifo[lifoId].empty()){
                                error = true;
                            }
                            else{
                                pairs[lifo[lifoId].top()].second = position;
                                lifo[lifoId].pop();
                            }
                        }
                    }
                }
            }
            if (error){
                cerr << "structure error at position " << ( position + 1 ) << endl;
                cerr << "character " << structElement << " unexpected, no corresponding opening bracket." << endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    for (unsigned int i = 0; i < 4; ++i){
        if ( !lifo[i].empty() ) {
            error = true;
            cerr << "Unbalanced pairing mask" << endl;
            while ( !lifo[i].empty() ) {
                cerr << "No matching pair for the open bracket at position " <<
                pairs[lifo[i].top()].first << endl;
                lifo[i].pop();
            }
        }
    }
    for ( map< char, stack<unsigned int> >::iterator iter = lifoMap.begin();
          iter != lifoMap.end(); ++ iter){
        if ( !(*iter).second.empty() ) {
            error = true;
            cerr << "Unbalanced pairing mask" << endl;
            while ( !(*iter).second.empty() ) {
                cerr << "No matching pair for the open bracket " << (*iter).first
                     << " at position " << pairs[(*iter).second.top()].first << endl;
                (*iter).second.pop();
            }
        }
    }
    for ( map< char, queue<unsigned int> >::iterator iter = fifoMap.begin();
          iter != fifoMap.end(); ++ iter){
        if ( !(*iter).second.empty() ) {
            error = true;
            cerr << "Unbalanced pairing mask" << endl;
            while ( !(*iter).second.empty() ) {
                cerr << "No matching pair for the open bracket " << (*iter).first
                     << " at position " << pairs[(*iter).second.front()].first << endl;
                (*iter).second.pop();
            }
        }
    }
    if (error){
        cerr << "Aborting..." << endl;
        exit(EXIT_FAILURE);
    }
}

SequenceTable::SequenceTable() {
    assert(0);
}

SequenceTable::SequenceTable( ParametersSet& sequenceParameters ) {
    ifstream ff;
    string filename = sequenceParameters.stringParameter( "Data file" );

    // Attempt to open the file
    ff.open( filename.c_str() );
    if ( !ff ) {
        cerr << "Unable to open file " << filename << " for reading" << endl;
        exit(EXIT_FAILURE);
    }
    init( ff, sequenceParameters );
}

void SequenceTable::categoryAttribution( UserCategoryAttribution userCategoryAttribution,
                                      stringstream& input ){
    //fill the user's category attribution vector if provided
    if ( userCategoryAttribution == CATEGORY_LINE ){
        for ( unsigned int i = 0; i < initialCategoryAttrib.size(); ++i ) {
            input >> ws;
            input >> initialCategoryAttrib[i];
            if (input.fail()){
                if (i==0) cerr << "Error: a class line was expected after your sequences" << endl;
                else cerr << "Error: incomplete class line after your sequences" << endl;
                exit(EXIT_FAILURE);
            }
            if (initialCategoryAttrib[i] < 1) {
                cerr << "Error for the class of the element "
                     << i+1 << ": " << initialCategoryAttrib[i] << endl;
                exit(EXIT_FAILURE);
            }
        }
        //check against the structure if provided
        if ( structure.size() != 0 ){
            checkCategoryAgainstStructure();
        }
    }
    //no attribution provided or automatic selection
    else{
        // 'only one category'
        for ( unsigned int i = 0; i < initialCategoryAttrib.size(); ++i){
                initialCategoryAttrib[i] = 1;
        }
        if ( userCategoryAttribution == AUTO ){
            //paired site for the second category
            for ( vector< pair<int,int>  >::iterator i = pairs.begin();
                  i != pairs.end(); ++i ){
                initialCategoryAttrib[(*i).first] = 2;
                initialCategoryAttrib[(*i).second] = 2;
            }
            // triplet site for the third category (or second if no doublet)
            unsigned int labelCategory = (pairs.size()==0) ? 2 : 3;
            for (unsigned int i = 0; i<triplets.size() ; ++i){
                initialCategoryAttrib[triplets[i]] = labelCategory;
                initialCategoryAttrib[triplets[i]+1] = labelCategory;
                initialCategoryAttrib[triplets[i]+2] = labelCategory;
            }
        }
    }
}

void SequenceTable::checkCategoryAgainstStructure(){
    // check for the pairs
    for ( vector< pair<int,int>  >::iterator iter = pairs.begin(); iter != pairs.end(); ++iter ){
        if ( initialCategoryAttrib[(*iter).first] != initialCategoryAttrib[(*iter).second] ){
            cerr << "class error for the pair (" << (*iter).first
                 << "," << (*iter).second << ")" << endl;
            cerr << "category 1 = " << initialCategoryAttrib[(*iter).first]
                 << ", category 2 = "
                 << initialCategoryAttrib[(*iter).second] << endl;
            exit(EXIT_FAILURE);
        }
    }
    // check for the triplets
    for (unsigned int i=0; i < triplets.size(); ++i ){
        if ( (initialCategoryAttrib[triplets[i]] != initialCategoryAttrib[triplets[i]+1]) ||
             (initialCategoryAttrib[triplets[i]] != initialCategoryAttrib[triplets[i]+2]) ){
            cerr << "class error for the codon ("
                 << triplets[i] << ", " << triplets[i]+1 << ", " << triplets[i]+2
                 << ")" << endl;
            cerr << "category 1 = " << initialCategoryAttrib[triplets[i]] << ", "
                 << "category 2 = " << initialCategoryAttrib[triplets[i]+1] << ", "
                 << "category 3 = " << initialCategoryAttrib[triplets[i]+2]
                 << endl;
            exit(EXIT_FAILURE);
        }
    }
}

void SequenceTable::init( ifstream & inputFile, ParametersSet& sequenceParameters ) {
    stringstream input;

    unsigned int numberSpecies = 0;
    unsigned int initialLength = 0;

    bool interleaved = sequenceParameters.boolParameter("Interleaved data file");

    UserCategoryAttribution userCategoryAttribution;
    if ( sequenceParameters.findParameter("Heterogeneous data models") ){
        string userAtt = sequenceParameters.stringParameter("Heterogeneous data models");
        if ( (userAtt=="auto") || (userAtt=="Auto") ){
            userCategoryAttribution = AUTO;
        }
        else{
            if ( (userAtt=="no") || (userAtt=="No") || (userAtt=="n") || (userAtt=="N") ){
                userCategoryAttribution = NONE;
            }
            else{
                if ( (userAtt=="yes") || (userAtt=="Yes") || (userAtt=="y") || (userAtt=="Y") ){
                     userCategoryAttribution = CATEGORY_LINE;
                }
                else{
                    cerr << "Error, yes, no or auto is expected for the field \"Heterogeneous data models\"" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    else{
        userCategoryAttribution = NONE;
    }
    
    bool userInterleavedStructure =
       ( sequenceParameters.findParameter("Interleaved structure") &&
         sequenceParameters.boolParameter("Interleaved structure") );

    if (userInterleavedStructure && !interleaved){
        cerr << "Check your control file" << endl;
        cerr << "If the the structure is interleaved, the data must be interleaved" << endl;
        exit(EXIT_FAILURE);
    }

    FileParser::skipComments( inputFile, input );

    // Read in header information
    // number of species, length of sequences, type ("DNA", "PROT", "CODON", "STRUCT")
    input >> ws;
    input >> numberSpecies;
    input >> ws;
    input >> initialLength;
    input >> ws;
    input >> type;

    initialSequences.resize( numberSpecies, initialLength );
    initialCategoryAttrib.resize( initialLength );

    if (type == "CODON"){
        // The initial length of sequences has to be divided by 3 for the codons data
        if (initialLength %3 != 0){
            cerr << "The length of sequences can not be divided by 3" << endl;
            exit(EXIT_FAILURE);
        }
        unsigned int nbCodons = initialLength/3;
        // if type is "CODON", there is no structure line, but the structure has to be 123 123 ...
        structure.resize(initialLength);
        for (unsigned int i = 0; i<nbCodons; ++i){
            structure[i*3]='1';
            structure[i*3+1]='2';
            structure[i*3+2]='3';
        }
        if (userCategoryAttribution==AUTO){
            cerr << "The token \"CODON\" in the datafile is incompatible with \"Heterogeneous data models = auto\"" << endl;
            cerr << "\"Heterogeneous data models = auto\" is used with the token STRUCT and a structure line is necessary" << endl;
            exit(EXIT_FAILURE);
        }
        if ( (userInterleavedStructure) && (userCategoryAttribution==NONE) ){
            cerr << "Error, the field \"Interleaved structure\" is set to yes but there should be no structure line nor class line" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else if ( type == "STRUCT" ){
        if (!userInterleavedStructure){
            // read structure
            for ( unsigned int i = 0; i < initialLength; ++i ) {
                char ch;
                input >> ws;
                input >> ch;
                putCharStructure(ch);
            }
        }
    }
    else if ( (type == "DNA") || (type == "PROT") ){
        if (userCategoryAttribution==AUTO){
            cerr << "The token in the datafile (\"DNA\" or \"PROT\") is incompatible with \"Heterogeneous data models = auto\"" << endl;
            cerr << "\"Heterogeneous data models = auto\" is used with the token STRUCT and a structure line is necessary" << endl;
            exit(EXIT_FAILURE);
        }
        if ( (userInterleavedStructure) && (userCategoryAttribution==NONE) ){
            cerr << "Error, the field \"Interleaved structure\" is set to yes but there is no structure line nor class line." << endl;
            cerr << "With the tokens DNA or PROT, the model line can be interleaved if you use \"Heterogeneous data models = yes\"" << endl;
            exit(EXIT_FAILURE);
        }
    }
    else{
        if ( type == "MIXED" ){
            cerr << "The old token MIXED must be replaced by the token STRUCT." << endl;
            cerr << "You also have to use \"Heterogeneous data models = auto\" in the DATAFILE section of your control-file to avoid "
                 << "having to write a model line. \"Heterogeneous data models = auto\" assigns automatically single nucleotide to the first model, "
                 << "pairs to the second model and codon to the third model." << endl;
            exit(EXIT_FAILURE);
        }
        cerr << "The third token in your data file should be the type of the file" << endl
             << "(either DNA, PROT, CODON or STRUCT)." << endl
             << "Please check the name of the data file you provided" << endl;
        exit(EXIT_FAILURE);
    }
    

    //fill the initialSequences variable with the sequences in the datafile
    if ( interleaved ){
        constructInterleaved(input, userCategoryAttribution, userInterleavedStructure);
    }
    else{
        constructNonInterleaved(input);
    }

    processStructure();

    // Category Attribution line is interleaved and so already read
    if( (userCategoryAttribution==CATEGORY_LINE) && userInterleavedStructure){
        checkCategoryAgainstStructure();
    }
    // No Category Attribution line or Category Attribution line at the end of the file
    else{
        categoryAttribution(userCategoryAttribution, input );
    }
    
    numberCategories = *(max_element(initialCategoryAttrib.begin(),initialCategoryAttrib.end()));

    constructRawSequences();
    constructSequences();
}

void SequenceTable::putCharStructure(char ch){
    if ( (ch!='(') && (ch!=')') && (ch!='.') && (ch!='-') && (ch!='1') && (ch!='2') && (ch!='3')
                   && (ch!='[') && (ch!=']') && (ch!='<') && (ch!='>') && (ch!='{') && (ch!='}') ){
        if ( (ch >= 'a') && (ch <= 'z') || (ch >= 'A') && (ch <= 'Z') ){
            cout << "WARNING: letters found in your consensus secondary structure. Make sure you read" << endl
                 << "the manual properly regarding their use. Starting with a upper case letter means" << endl
                 << "conventional nested pairing (12345.....54321). The other pairing rule is used" << endl
                 << "when starting with a lower case letter (123.4567......12..34567)" << endl;
        }
        else{
            cerr << "Unrecognized symbol in the structure (" << ch << ")." << endl;
            exit(EXIT_FAILURE);
        }
    }
    structure.push_back(ch);
}



// space and alpha, some unary function to help STL algorithms below
struct space : public unary_function<char, bool> {
    double operator()(char x) const { return isspace(x); }
};
struct alpha : public unary_function<char, bool> {
    double operator()(char x) const { return isalpha(x); }
};
struct num : public unary_function<char, bool> {
    double operator()(char x) const { return isdigit(x); }
};


// Constructs a tree given the data in the file is interleaved
void SequenceTable::constructInterleaved(stringstream & input, UserCategoryAttribution userCategoryAttribution, bool userInterleavedStructure) {
    int block=0;
    unsigned int numberSpecies = initialSequences.numberRows();
    unsigned int sequencesLength = initialSequences.numberColumns();
    unsigned int totalLineLength = 0;

    string readLine;
    string structLabel;
    string modelLabel;

    input >> ws;

    while (!input.eof() && totalLineLength<sequencesLength ){
        unsigned int blockLineLength = 0;

        // Read structure line
        if (userInterleavedStructure){
            getline(input,readLine);
            //declare an iterator on the string and look for first non-space character
            string::iterator iter = find_if(readLine.begin(),readLine.end(), not1(space()) );
            //first block? do we have a label on the structure?
            if (block == 0){
                //Is there a white space close in the line that could mark the end of a name for the structure?
                //ie, less than 20 characters, most of them letter
                string::iterator iter2 = find_if( iter, iter+21, space() );
                if ( (distance(iter,iter2)<20) && (distance(iter,iter2)>3) ){
                    string forbid=".-()[]{}<>123";
                    //if we do not find a forbidden character in the potential label then use it
                    if ( find_first_of( iter, iter2, forbid.begin(), forbid.end() ) == iter2 ){
                         structLabel.assign(iter,iter2);
                         iter = find_if( iter2, readLine.end(), not1(space()) );
                    }
                }
            }
            else{
                //if a label is defined at block 0
                if ( structLabel.length() ){
                    string::iterator iter2 = iter + structLabel.length();
                    if (isspace(*iter2)){
                        //and the "word" beginning at iter is the same then jump it
                        if ( equal( structLabel.begin(), structLabel.end(), iter ) ){
                            iter = find_if( iter2, readLine.end(), not1(space()) );
                        }
                    }
                }
            }
            // Read structure
            while ( iter!=readLine.end() ){
                if ( !isspace( *iter ) ) {
                    putCharStructure(*iter);
                    ++ blockLineLength;
                }
                ++iter;
            }
        }

        // Read data for each specie
        for(unsigned int i=0; i<numberSpecies; ++i){
            unsigned int lineLength = 0;
            input >> ws;
            getline(input,readLine);
            //declare an iterator on the string and look for first non-space character
            string::iterator iter = find_if( readLine.begin(),readLine.end(), not1(space()) );
            // if first block
            if (block == 0){
                string::iterator iter2 = find_if( iter, iter+51, space() );
                if ( distance(iter,iter2) == 51 ){
                    cerr << "Error, expecting a species name (less than 50 characters) in the first block, for species " << i+1 << endl;
                    cerr << "Got \"" << string(iter,iter2) << "\" instead" << endl;
                    exit(EXIT_FAILURE);
                }
                string label( iter, iter2 );
                //case-insensitive name checking (to avoid dangerous duplicate)
                int prev = findSequence(label);
                if (prev != -1){
                    cerr << "Error: duplicate taxon name: " << label << "; taxons " << prev+1 << " and " << i+1 << endl;
                    exit(EXIT_FAILURE);
                }
                species.push_back(label);
                iter = find_if( iter2, readLine.end(), not1(space()) );
            }
            // if not first block
            else{
                string::iterator iter2 = iter + species[i].length();
                if (isspace(*iter2)){
                    //and the "word" beginnind at iter is the same then jump it
                    if ( equal( species[i].begin(), species[i].end(), iter ) ){
                        iter = find_if( iter2, readLine.end(), not1(space()) );
                    }
                }
            }
            // Read data
            while (iter!=readLine.end()) {
                char next = *iter;
                if ( !isspace( next ) ) {
                    initialSequences(i, totalLineLength+lineLength) = next;
                    ++lineLength;
                }
                ++iter;
            }

            // No structure line in the block and first species
            if ( (blockLineLength==0) && (i==0) ){
                blockLineLength = lineLength;
            }
            else if(lineLength!=blockLineLength){
                cerr << "Error in block "<< block+1 << ", species " << i+1
                     << ". The lengths do not match." << endl;
                cerr << "Please check this line and also check for possible mismatches of species label in this block." << endl;
                cerr << "Expected length for this block: " << blockLineLength << ";   Length for the species: " << lineLength << endl;
                exit(EXIT_FAILURE);
            }
        }

        // Read category line
        if ( (userCategoryAttribution==CATEGORY_LINE) && userInterleavedStructure){
            input >> ws;
            getline(input,readLine);
            //declare an iterator on the string and look for first non-space character
            string::iterator iter = find_if(readLine.begin(),readLine.end(), not1(space()) );
            //first block? do we have a label on the model line?
            if (block == 0){
                //Is there a white space close in the line that could mark the end of a name for the structure?
                //ie, less than 20 characters, most of them letter
                string::iterator iter2 = find_if( iter, iter+21, space() );
                if ( (distance(iter,iter2)<20) && (distance(iter,iter2)>3) ){
                    //if we do not find a forbidden digit in the potential label then consider it label
                    if( count_if( iter, iter2, not1(num()) ) == ( distance(iter,iter2) ) ){
                        modelLabel.assign(iter,iter2);
                        iter = find_if( iter2, readLine.end(), not1(space()) );
                    }
                }
            }
            else{
                //if a label is defined at block 0
                if ( modelLabel.length() ){
                    string::iterator iter2 = iter + modelLabel.length();
                    if ( isspace(*iter2) ){
                        //and the "word" beginnind at iter is the same then jump it
                        if ( equal( modelLabel.begin(), modelLabel.end(), iter ) ){
                            iter = find_if( iter2, readLine.end(), not1(space()) );
                        }
                    }
                }
            }

            // Read category
            unsigned int lineLength = 0;
            char* pchar = &(*iter);
            while ( *pchar ){
                //read the number (not pretty but quick...)
                int nb=0;
                sscanf( pchar, "%d%n", &(initialCategoryAttrib[totalLineLength+lineLength]), &nb );
                pchar += nb;
                if ( (initialCategoryAttrib[totalLineLength+lineLength]<1) || !nb ){
                    cerr << "Error in the model line of block " << block+1 << "(character " <<  lineLength+1 << ")." << endl;
                    exit(EXIT_FAILURE);
                }
                ++lineLength;
                while (isspace(*pchar)) ++pchar;
            }
            if(lineLength!=blockLineLength){
                cerr << "Error in block "<< block+1 << endl;
                cerr << "The length of the model line does not match" << endl;
                exit(EXIT_FAILURE);
            }
        }
        totalLineLength += blockLineLength;
        ++block;
        input >> ws;
    }

    if ( totalLineLength != sequencesLength ){
        cerr << "Error: " << totalLineLength
             << " nucleotides were found in each sequence but it does not match with the declared length: "
             << sequencesLength << endl;
        exit(EXIT_FAILURE);
    }
}

// Constructs a tree given the data in the file is non-interleaved
void SequenceTable::constructNonInterleaved( stringstream & input ) {

    unsigned int numberSpecies = initialSequences.numberRows();
    unsigned int sequencesLength = initialSequences.numberColumns();

    char endline;

    // Read in species and sequences
    for ( unsigned int i = 0; i < numberSpecies; ++i ) {
        // Skip white spaces
        input >> ws;
        string speciesName;
        // Read in species "name"
        if( input.eof() ){
                cerr << "Unexpected end of sequence file when parsing species "
                     << i + 1 << "name. Exit..." << endl;
                exit(EXIT_FAILURE);
        }
        input >> speciesName;
        //case-insensitive name checking (to avoid dangerous duplicate)
        int prev = findSequence(speciesName);
        if (prev != -1){
            cerr << "Error: duplicate taxon name: " << speciesName << "; taxons " << prev+1 << " and " << i+1 << endl;
            exit(EXIT_FAILURE);
        }
        species.push_back(speciesName);
        // Read in sequences
        for ( unsigned int j = 0; j < sequencesLength; ++j ) {
            input >> ws; // Skip white spaces
            if( input.eof() ){
                cerr << "Unexpected end of sequence file when parsing species "
                     << speciesName << ". Exit..." << endl;
                exit(EXIT_FAILURE);
            }
            input >> initialSequences(i, j);
        }
        input.get( endline );
        if ( !isspace(endline) ){
            cerr << "End of line after " << speciesName
                 << " sequence not found. Exit..." << endl;
            exit(EXIT_FAILURE);
        }
    }
}

void SequenceTable::constructRawSequences(){
    unsigned int nbColumnsInRawSequences;

    nbColumnsInRawSequences = initialSequences.numberColumns() - pairs.size() - 2*triplets.size() ;

    rawSequences.resize( initialSequences.numberRows(),nbColumnsInRawSequences);
    //we will apply the RNA mask to initialCategoryAttrib too
    categoryAttrib.resize(nbColumnsInRawSequences);

    initial2Raw.resize( initialSequences.numberColumns() );
    raw2Initial.resize(nbColumnsInRawSequences);

    vector< pair<int,int> >::iterator iterPairs = pairs.begin();
    vector< int >::iterator iterTriplets = triplets.begin();
    set<unsigned int> alreadyDone;
    int currentCol = 0;
    for (unsigned int j = 0; j < initialSequences.numberColumns(); ++j){
        set<unsigned int>::iterator jSearch = alreadyDone.find(j);
        //if the nucleotid has not been treated before
        if ( jSearch == alreadyDone.end() ){
            // The symbol is one nucleotide (if it is not a pair, neither a triplet)
            if ( ( (iterPairs==pairs.end()) || ((*iterPairs).first != (int)j) ) &&
                 ( (iterTriplets==triplets.end()) || ((*iterTriplets) != (int)j) ) ){
                initial2Raw[j] = currentCol;
                raw2Initial[currentCol].resize(1);
                raw2Initial[currentCol][0] = j;
                for (unsigned int i = 0; i < initialSequences.numberRows(); ++i){
                    rawSequences(i,currentCol) = initialSequences(i,j);
                }
                categoryAttrib[currentCol] = initialCategoryAttrib[j]-1;
            }
            // The symbol is two nucleotides (if it is not a triplet)
            else if ( (iterTriplets==triplets.end()) || ( (*iterTriplets) != (int)j) ){
                for (unsigned int i = 0; i < initialSequences.numberRows(); ++i){
                    rawSequences(i,currentCol).resize(2);
                    rawSequences(i,currentCol)[0] = initialSequences(i,j);
                    rawSequences(i,currentCol)[1] = initialSequences(i,(*iterPairs).second);
                }
                initial2Raw[j] = currentCol;
                initial2Raw[(*iterPairs).second] = currentCol;
                raw2Initial[currentCol].resize(2);
                raw2Initial[currentCol][0] = j;
                raw2Initial[currentCol][1] = (*iterPairs).second;
                categoryAttrib[currentCol] = initialCategoryAttrib[j]-1;
                //double-check structure vs. category attribution
                //already done before
                assert(categoryAttrib[currentCol] ==
                          initialCategoryAttrib[(*iterPairs).second]-1);
                alreadyDone.insert((*iterPairs).second);
                ++iterPairs;
            }
            // The symbol is three nucleotides
            else{
                for (unsigned int i = 0; i < initialSequences.numberRows(); ++i){
                    rawSequences(i,currentCol).resize(3);
                    rawSequences(i,currentCol)[0] = initialSequences(i,j);
                    rawSequences(i,currentCol)[1] = initialSequences(i,j+1);
                    rawSequences(i,currentCol)[2] = initialSequences(i,j+2);
                }
                initial2Raw[j] = currentCol;
                initial2Raw[j+1] = currentCol;
                initial2Raw[j+2] = currentCol;
                raw2Initial[currentCol].resize(3);
                raw2Initial[currentCol][0]=j;
                raw2Initial[currentCol][1]=j+1;
                raw2Initial[currentCol][2]=j+2;
                categoryAttrib[currentCol] = initialCategoryAttrib[j]-1;
                alreadyDone.insert(j+1);
                alreadyDone.insert(j+2);
                ++iterTriplets;
            }
            ++currentCol;
        }
        //if already treated, remove it from the set
        else{
            alreadyDone.erase( jSearch );
        }
    }
    assert( alreadyDone.empty() );
}

void SequenceTable::constructSequences(){

    sequences.resize( numberCategories );
    invariantBases.resize( numberCategories );

    bool invar;
    unsigned int speciesId;

    raw2Final.resize(rawSequences.numberColumns());
    final2Raw.resize(numberCategories);
    invariant2Raw.resize(numberCategories);
    for ( unsigned int i = 0; i < rawSequences.numberColumns() ; ++i ){
        invar = true;
        speciesId = 0;
        while( invar && ( speciesId<species.size() ) ){
            invar = ( rawSequences(speciesId, i) == rawSequences(0, i) );
            ++speciesId;
        }
        if (!invar){
            raw2Final[i] = pair<unsigned int, int>( categoryAttrib[i],
                                        (int)final2Raw[categoryAttrib[i]].size() );
            final2Raw[categoryAttrib[i]].push_back(i);
        }
        //new invariant symbol
        else{
            unsigned int j = 0;
            bool found = false;
            while ( ( j < invariantBases[categoryAttrib[i]].size() ) && !found ){
                if ( invariantBases[categoryAttrib[i]][j].first == rawSequences( 0, i ) ){
                    ++invariantBases[categoryAttrib[i]][j].second;
                    invariant2Raw[categoryAttrib[i]][j].push_back(i);
                    raw2Final[i] = pair<unsigned int, int>(categoryAttrib[i],-(j+1));
                    found = true;
                }
                ++j;
            }
            if (!found){
                pair<string,int> newInvariant( rawSequences( 0, i ), 1);
                invariantBases[categoryAttrib[i]].push_back(newInvariant);
                vector< unsigned int > vec(1);
                vec[0] = i;
                invariant2Raw[categoryAttrib[i]].push_back( vec );
                raw2Final[i] = pair<unsigned int, int>(categoryAttrib[i],
                                  -(invariantBases[categoryAttrib[i]].size()));
            }
        }
    }
    for ( unsigned int categoryId = 0; categoryId < numberCategories; ++categoryId ){
        sequences[categoryId].resize(species.size(),final2Raw[categoryId].size());
        for ( unsigned int i = 0; i < final2Raw[categoryId].size() ; ++i ){
            for ( unsigned int j = 0; j < species.size(); ++j ){
                sequences[categoryId](j, i) = rawSequences( j, final2Raw[categoryId][i] );
            }
        }
    }
}

SequenceTable::~SequenceTable() {
}

int SequenceTable::findSequence( string label ) {

    string temp;

    transform( label.begin(), label.end(), label.begin(), (int(*)(int))tolower );

    for ( unsigned int i = 0; i < species.size(); ++i ) {
        temp = species[i];
        transform( temp.begin(), temp.end(), temp.begin(), (int(*)(int))tolower );
        if ( temp == label ) return i;
    }
    return ( -1 );
}

void SequenceTable::print( ostream & output ) {

    output << "sequence type : " << type << endl;
    output << "number of species : " << species.size() << endl;
    output << "number of categories used : " << sequences.size() << endl;
    for ( unsigned int i = 0; i < species.size(); ++i ) {
        output << "species "<< i + 1 << " : " << species[i] << endl;
        for ( unsigned int j = 0; j < numberCategories; ++j ) {
            output << "category "<< j+1 << ", " << sequences[j].numberColumns()
                   << " non-invariant symbols" << endl;
            for ( unsigned int k = 0; k < sequences[j].numberColumns(); ++k ) {
                output << sequences[j](i,k);
            }
            output << endl;
        }
        output << endl;
    }
    output << "invariants " << endl;
    for ( unsigned int j = 0; j < numberCategories; ++j ) {
        output << "category "<< j+1 << " : ";
        for ( vector< pair<string, unsigned int> > ::iterator iter =
              invariantBases[j].begin(); iter != invariantBases[j].end();
              ++iter) {
            output << (*iter).first << "-" << (*iter).second << "  ";
        }
        output << endl;
    }
    output << endl << "categories" << endl;
    for (unsigned int i=0; i<categoryAttrib.size(); ++i){
        output << categoryAttrib[i] ;
    }
    output << endl;
}

vector < vector< unsigned int > > SequenceTable::retrieveInitialNucleotides( int categoryId,
               int symbolIndex ){
    vector < vector< unsigned int > > ret;
    if (symbolIndex >= 0){
        assert((unsigned int)symbolIndex < final2Raw[categoryId].size());
        unsigned int rawIndex = final2Raw[categoryId][symbolIndex];
        ret.push_back( raw2Initial[rawIndex] );
    }
    else {
        symbolIndex = -1-symbolIndex;
        assert((unsigned int)symbolIndex < invariantBases[categoryId].size());
        for ( vector<unsigned int>::iterator invSymbolIterator =
                invariant2Raw[categoryId][symbolIndex].begin();
              invSymbolIterator != invariant2Raw[categoryId][symbolIndex].end();
              ++invSymbolIterator ){
            ret.push_back( raw2Initial[*invSymbolIterator] );
        }
    }
    return ret;
}

vector< unsigned int > SequenceTable::retrieveRawNucleotides( int categoryId,
               int symbolIndex ){
    vector< unsigned int > ret;
    if (symbolIndex >= 0){
        ret.push_back(final2Raw[categoryId][symbolIndex]);
    }
    else{
        symbolIndex = -1-symbolIndex;
        assert((unsigned int)symbolIndex < invariantBases[categoryId].size());
        for ( vector<unsigned int>::iterator invSymbolIterator =
                invariant2Raw[categoryId][symbolIndex].begin();
              invSymbolIterator != invariant2Raw[categoryId][symbolIndex].end();
              ++invSymbolIterator ){
            ret.push_back( *invSymbolIterator );
        }
    }
    return ret;
}

void SequenceTable::printStats( ostream & output, vector<int>, vector<int> ){

    output << "initialSequences" << endl;
    output << "... number of Rows: " << initialSequences.numberRows() << endl;
    output << "... number of Columns: " << initialSequences.numberColumns() << endl;

    output << "rawSequences" << endl;
    output << "... number of Rows: " << rawSequences.numberRows() << endl;
    output << "... number of Columns: " << rawSequences.numberColumns() << endl;

    output << "sequences" << endl;
    output << "... number of models" << sequences.size() << endl;
    for (unsigned int i=0; i<sequences.size(); ++i){
        output << "... model " << i << ": " << endl;
        output << "...... number of Raws: " << sequences[i].numberRows() <<endl;
        output << "...... number of Columns: " << sequences[i].numberColumns() <<endl;
    }

    output << "structure size: " << structure.size() <<endl;
    output << "pairs size: " << pairs.size() <<endl;
    output << "triplets size: " << triplets.size() <<endl;
    output << "initialCategoryAttrib - size: " << initialCategoryAttrib.size() <<endl;
    output << "categoryAttrib size: " << categoryAttrib.size() <<endl;
}
