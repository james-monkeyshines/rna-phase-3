#include "Models/CodonParent.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

using namespace std;

CodonParent::CodonParent(){
}

CodonParent::CodonParent( const string & registrationName ) :
        MatrixModel( registrationName ){
}

CodonParent::CodonParent( ParametersSet & parameters ) : MatrixModel( parameters ){
    aminoAcid.reserve(20);
    aminoAcid.push_back('A');
    aminoAcid.push_back('R');
    aminoAcid.push_back('N');
    aminoAcid.push_back('D');
    aminoAcid.push_back('C');
    aminoAcid.push_back('Q');
    aminoAcid.push_back('E');
    aminoAcid.push_back('G');
    aminoAcid.push_back('H');
    aminoAcid.push_back('I');
    aminoAcid.push_back('L');
    aminoAcid.push_back('K');
    aminoAcid.push_back('M');
    aminoAcid.push_back('F');
    aminoAcid.push_back('P');
    aminoAcid.push_back('S');
    aminoAcid.push_back('T');
    aminoAcid.push_back('W');
    aminoAcid.push_back('Y');
    aminoAcid.push_back('V');
}

void CodonParent::readGeneticCode(string fileName){
    ifstream fileStream;
    
    // Attempt to open the file
    fileStream.open( fileName.c_str() );
    if ( !fileStream ) {
        cerr << "Unable to open file " << fileName << " for reading" << endl;
        exit(EXIT_FAILURE);
    }
    
    while(!fileStream.eof()){
        string codon;
        char aa;
        fileStream >> codon;
        if (codon.length()!=3){
            cerr << "Invalid codon \"" << codon << "\" in your genetic code" << endl;
        }
        for ( unsigned int n = 0; n < 3; ++n ){
            unsigned int i = getSingleSymbolNumber(codon[n]);
            if ( (i!=0x01) && (i!=0x02) && (i!=0x04) && (i!=0x08) ){
                cerr << "Invalid codon \"" << codon << "\" in your genetic code" << endl;
                exit(EXIT_FAILURE);
            }
        }
        fileStream >> ws;
        fileStream >> aa;
        if ( find(aminoAcid.begin(), aminoAcid.end(), aa) == aminoAcid.end() ){
            cerr << "Unrecognized amino-acid \"" << aa << "\" in your genetic code" << endl;
            exit(EXIT_FAILURE);
        }
        geneticCode.push_back( pair <string, char> (codon, aa) );
        fileStream >> ws;
    }
}

CodonParent::~CodonParent() {
}


int CodonParent::getSingleSymbolNumber(const char & base) const {
    //for efficiency reason we use binary bits and mask there
    switch (base) {
        case 'a': case 'A':
            return 0x01;
        case 'c': case 'C':
            return 0x02;
        case 'g': case 'G':
            return 0x04;
        case 'u': case 't': case 'U': case 'T':
            return 0x08;
        case 'r': case 'R':
            return 0x05; // Unknown purine       (A or G)
        case 'y': case 'Y':
            return 0x0A; // Unknown pyramidine   (C or U)
        case 'm': case 'M':
            return 0x03; //                      (A or C)
        case 'k': case 'K':
            return 0x0C; //                      (U or G)
        case 'w': case 'W':
            return 0X09; //                      (U or A)
        case 's': case 'S':
            return 0X06; //                      (C or G)
        case 'b': case 'B':
            return 0X0E; //                     (U, C or G)
        case 'd': case 'D':
            return 0X0D; //                     (U, A or G)
        case 'h': case 'H':
            return 0x0B; //                     (U, A or C)
        case 'v': case 'V':
            return 0X07; //                     (A, C or G)
        case 'n': case 'x': case '?': case 'N': case 'X': case '-':
            return 0x0F; // Unknown nucleotide  (A, C, G or T)
        default :
            return -1; // Out of range (unrecognized symbol)
    }
}

string CodonParent::getSingleSymbol( unsigned int symbolNumber) const{

    switch ( symbolNumber ) {
        case 0x01 :
            return ( string( "A" ) );
        case 0x02 :
            return ( string( "C" ) );
        case 0x04 :
            return ( string( "G" ) );
        case 0x08 :
            return ( string( "U" ) );
        case 0x05 :
            return ( string( "R" ) ); // Unknown Purine     (A or G)
        case 0x0A :
            return ( string( "Y" ) ); // Unknown Pyramidine (C or T)
        case 0x03 :
            return ( string( "M" ) ); //                    (A or C)
        case 0x0C :
            return ( string( "K" ) ); //                    (U or G)
        case 0x09 :
            return ( string( "W" ) ); //                    (U or A)
        case 0x06 :
            return ( string( "S" ) ); //                    (C or G)
        case 0x0E :
            return ( string( "B" ) ); //                 (C, G or T)
        case 0x0D :
            return ( string( "D" ) ); //                 (A, G or T)
        case 0x0B :
            return ( string( "H" ) ); //                 (A, C or T)
        case 0x07 :
            return ( string( "V" ) ); //                 (A, C or G)
        case 0x0F :
            return ( string( "N" ) ); // Unknown    (A , C , G or T)
        default :
            return ( string( "Z" ) ); // Unrecognized symbol
    }
}



int CodonParent::getSymbolNumber( const string & codonSymbol, unsigned int ) const {
    if (codonSymbol.length() != 3) return getNumberSymbols();

    int nuc1 = getSingleSymbolNumber(codonSymbol[0]) - 1 ;
    if (nuc1==-2) return getNumberSymbols();
    int nuc2 = getSingleSymbolNumber(codonSymbol[1]) - 1 ;
    if (nuc2==-2) return getNumberSymbols();
    int nuc3 = getSingleSymbolNumber(codonSymbol[2]) - 1 ;
    if (nuc3==-2) return getNumberSymbols();

    return (nuc1*225+nuc2*15+nuc3);
}

unsigned int CodonParent::getMask( unsigned int symbolNumber ) const {
    unsigned int n1 = symbolNumber/225;
    unsigned int r = symbolNumber%225;    
    unsigned int n2 = r/15;
    unsigned int n3 = r%15;
    return (((n1+1)<<8) | ((n2+1)<<4) | (n3+1));
}

string CodonParent::getSymbol( unsigned int symbolNumber, unsigned int ) const {

    if (symbolNumber>=getNumberSymbols()) return "ZZZ";
    
    unsigned int q = symbolNumber/225;
    unsigned int r = symbolNumber%225;
    string s = getSingleSymbol(q+1);
    s += getSingleSymbol(r/15+1);
    s += getSingleSymbol(r%15+1);
    
    return s;
}

