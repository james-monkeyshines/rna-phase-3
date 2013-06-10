#include "Util/FileParser.h"

#include "Util/ParametersSet.h"


using namespace std;

#include <iostream>
#include <sstream>
#include <vector>
#include <stack>
#include <cstdlib>

FileParser::FileParser( const string & fileName ) :
inputFile( fileName.c_str() ) {
    if ( !inputFile ) {
        cerr << fileName << " does not exist ! exit..." << endl;
        exit(EXIT_FAILURE);
    }
    this->fileName = fileName;
}

FileParser::~FileParser() {
    inputFile.close();
}

ParametersSet* FileParser::retrieveParametersSet() {

    char parameterLabel[160];
    char parameterValue[500];

    ParametersSet* parameters = new ParametersSet( fileName );
    stack <  ParametersSet * > categories;
    ParametersSet * currentCategory = parameters;

    stringstream input;
    skipComments( inputFile, input );

    int line = 0;

    int index;

    while ( !input.eof() ) {
        index = 0;
        do {
            input.get( parameterLabel[0] );
            if ( parameterLabel[0] == '\n' ) {
                ++line;
            }
        }
        while ( isspace( parameterLabel[0] ) );
        // read the label
        if ( parameterLabel[0] != 0 ) {
            do {
                input.get( parameterLabel[++index] );
            }
            while ( ( parameterLabel[index] != '=' ) &&
            ( parameterLabel[index] != '\n' ) &&
            parameterLabel[index] != 0 );

            if ( parameterLabel[index] == 0 ) {
                parameterLabel[index] = '\n';
            }
            switch ( parameterLabel[index] ) {
                case '\n' : {
                        ++line;
                        // category command
                        if ( ( parameterLabel[0] ) != '{' ) {
                            cerr << "Malformed file " << fileName << " line : " <<
                            line << endl ;
                            delete parameters;
                            return NULL;
                        }
                        while ( parameterLabel[--index] != '}' ) {
                            if ( index == 0 ) {
                                cerr << "Malformed file " << fileName <<
                                " line : " << line << endl ;
                                delete parameters;
                                return NULL;
                            }
                        }
                        //null terminate the character string, erase the '}'
                        parameterLabel[index] = 0;
                        //end of category
                        if ( parameterLabel[1] == '\\' ) {
                            if ( ( ( int )( currentCategory->getName() ).length()
                            != index - 2 ) ||
                            ( !( currentCategory->getName()
                                 ).compare( 2, index - 2, parameterLabel ) ) ) {
                                     cerr << "Malformed file " << fileName <<
                                     " line : " << line << endl ;
                                     cerr << "Only the " <<
                                     currentCategory->getName() <<
                                     " category can end there" << endl;
                                     delete parameters;
                                     return NULL;
                            }
                            //ok, end of category
                            else {
                                currentCategory = categories.top();
                                categories.pop();
                            }
                        }
                        //begin a category
                        else {
//                            cout << "Add category : " << string( parameterLabel +
//                            1, index - 1 ) << "." << endl;
                            categories.push( currentCategory );
                            currentCategory = &( 
                                  (*currentCategory)(string(parameterLabel+1,index-1)) );
                        }
                    }
                    break;
                case '=' : {
                        ++line;
                        //cut the parameter label character string at its end (remove
                        // the '=' and the space at the end of the name
                        parameterLabel[index] = 0;
                        while ( isspace( parameterLabel[--index] ) ) {
                            parameterLabel[index] = 0;
                        }
                        index = 0;
                        do {
                            input.get( parameterValue[0] );
                        }
                        while ( isspace( parameterValue[0] ) );
                        //read the parameter value
                        do {
                            input.get( parameterValue[++index] );
                        }
                        while ( ( parameterValue[index] != '\n' ) &&
                        ( parameterValue[index] != 0 ) );
                        //cut the parameterValue character string at its end (remove
                        // the '\n' and the space at the end of the name
                        if ( parameterValue[index] == '\n' ) {
                            parameterValue[index] = 0;
                        }
                        while ( isspace( parameterValue[--index] ) ) {
                            parameterValue[index] = 0;
                        }
//                        cout << "line " << line << ":" << parameterLabel <<
//                        " = " << parameterValue << "." << endl ;
                        //add the parameter to the map
                        ( ( * currentCategory )
                        [parameterLabel] ).assign( parameterValue );
                    }
                    break;
                default : {
                        cerr << "Malformed file " << fileName << " line : " <<
                        line << endl ;
                        delete parameters;
                        return NULL;
                    }
            }
        }
    }
    if( categories.size() != 0 ){
        cerr << "Malformed file " << fileName << ':'
             << " unclosed block(s)" << endl;
        delete parameters;
        return NULL;
    }
    return parameters;
}


void FileParser::skipComments( ifstream & inputFile,
stringstream & inputStream ) {
    char junk;

    while ( !inputFile.eof() ) {
        junk = 0;
        while ( ( junk != '\n' ) && ( !inputFile.eof() ) ) {
            inputFile.get( junk ) ;
            //skip commentary
            if ( junk == '#' ) {
                while ( ( junk != '\n' ) && ( !inputFile.eof() ) ) {
                    inputFile.get( junk );
                }
                inputStream.put('\n');
            }
            else {
                inputStream.put(junk);
            }
        }
    }
    inputStream.seekg( 0 );
    inputFile.seekg( 0 );
}

void FileParser::readTree( istream & inputFile, string& stringTree ) {
  
    stringTree.clear();
    char next = ' ';
    int parenthesisCount = 0;

    if (!inputFile.good()){
        cerr << "WARNING : error while reading a tree" << endl;
        stringTree.clear();
        return;
    }
    inputFile >> ws;
    while( parenthesisCount ||
           ( (next != '\n') && !inputFile.eof()) ){
        inputFile.get( next );
        if ( inputFile.eof() ){
            cerr << "WARNING : error while reading a tree" << endl;
            stringTree.clear();
            return;
        }
        switch(next){
            case '(' :
                ++parenthesisCount;
                break;
            case ')' :
                if ( !parenthesisCount ){
                    cerr << "WARNING : error while reading a tree, unmatched ')'" << endl;
                    stringTree.clear();
                    return;
                }
                --parenthesisCount;
                break;
            //exit without adding the ';' at the tree
            case ';' :
                if ( parenthesisCount ){
                    cerr << "WARNING : error while reading a tree, unmatched '('" << endl;
                    stringTree.clear();
                }
                return;
        }
        stringTree += next;
        if (parenthesisCount){
            inputFile >> ws;
        }
   }
   
   if(parenthesisCount){
        cerr << "WARNING : error while reading a tree" << endl;
        stringTree.clear();
   }
}


void FileParser::confirmOpenFile( ofstream & file, const string & filename,
                                  bool overwrite ) {
    
    string answer ;

    ifstream temp( filename.c_str() );
    // Find out if the file with name "filename" already exists and if it does
    // confirm whether user wants to overwrite the file; unless we are not in
	// 'cautious' mode, in which case files are automatically overwritten.
    if ( temp.is_open() ) {
        if (overwrite){
			if (cautious) {
            	cout << "The file \"" << filename << "\" already exists" << endl;
            	cout << "Overwrite (yes/no) ? " ;
            	cin >> answer  ;
            	if ( ( answer[0] == 'y' ) || ( answer[0] == 'Y' ) ) {
                	temp.close() ;
                	file.open( filename.c_str(), ifstream::out );
            	}
            	else {
            	    cout << "Exiting" << endl ;
            	    exit ( 0 ) ;
            	}
			} else {
				temp.close() ;
                file.open( filename.c_str(), ifstream::out );
			}
        }
        else{ //append
            file.open( filename.c_str(), ofstream::in | ofstream::out | ofstream::ate );
        }
    }
    else{
        file.open( filename.c_str(), ofstream::out ) ;
    }

    if ( !file.is_open() ) {
        cerr << "Unable to create file \"" << filename << "\"" << endl ;
        exit ( 1 ) ;
    }
}




