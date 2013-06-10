#ifndef FILEPARSER_H
#define FILEPARSER_H

#include <string>
#include <fstream>
#include <sstream>
using namespace std;

const bool cautious = 1;

class ParametersSet;

class FileParser {
public:

    /** .********************************************************************
     * FileParser
     * @input      the input stream to parse
     * @semantics  construct a FileParser object to read the content of the
     *             input stream
     ************************************************************************ */
    FileParser( const string & fileName );

    /** ************************************************************************
     * ~FileParser
     * @semantics  destructor, close the file
     ************************************************************************ */
    ~FileParser();

    /** ************************************************************************
     * retrieveParametersSet
     * @semantics  this routine parses the content of the input stream and
     *             return a ParametersSet instance
     ************************************************************************ */
    ParametersSet* retrieveParametersSet();

    /** ************************************************************************
     * skipComments
     * @input      ifstream, the file stream to be filtered
     * @return     a stringstream, a copy of the content of the file
     *             without comment
     * @semantics  skip all whitespace and skip all lines
     *             beginning with a hash (#) symbol
     ************************************************************************ */
    static void skipComments( ifstream & inputFile, stringstream & inputStream );

    /** ************************************************************************
     * readTree
     * @input      istream, the stream from where the tree should be read
     * @return     the string representation of the tree (without space and \n)
     * @semantics  read (and check) the string representation of a tree from
     *             an input stream
     ************************************************************************ */
    static void readTree( istream & inputFile, string& stringTree );
    
    /** ************************************************************************
     * confirmOpenFile
     * @input      fileName, the name of the file to open in writing mode
     * @input      overwrite, (false = append,
     *                         true = overwrite with confirmation)
     * @output     file, an output file stream
     * @semantics  open a file for writing, if the file already exists ask for
     *             a confirmation from the user
     ************************************************************************ */
    static void confirmOpenFile( ofstream & file , const string & fileName,
                                 bool overwrite );

private:
    ifstream inputFile;
    string fileName;
};

#endif //FILEPARSER_H




