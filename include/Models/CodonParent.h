#ifndef CODONPARENT_H
#define CODONPARENT_H

#include "Models/MatrixModel.h"


using namespace std;

/* Ancestral class for codon and amino-acid substitution models
 * deals with procedure related to the genetic code */

class CodonParent: public MatrixModel{

protected:
    
    //geneticCode is a vector of string+char
    //it contains the mapping codon->amino-acid
    vector< pair<string,char> > geneticCode;
    
    // list of the char coding for the amino acid in the usual order
    vector<char> aminoAcid;


    CodonParent( const string & registrationName );
    CodonParent( ParametersSet & parameter );
    CodonParent();
    
    ~CodonParent();

    /** ************************************************************************
     * readGeneticCode()
     * @input     fileName, name of the genetic code file
     * @semantics read the genetic code that the user wants to use in a special
     *            file
     ************************************************************************ */
    void readGeneticCode(string fileName);
    
    /** ************************************************************************
     * getNumberSymbols
     * @return  The number of distinct symbols recognized by the model
     ************************************************************************ */
    virtual unsigned int getNumberSymbols( unsigned int = 0 ) const{
        return 3375;
    }
    
    /** ************************************************************************
     * getSymbolNumber
     * @input   base, a string representing a symbol of the model
     * @return  The basenumber of the symbol
     *           e.g. "AAA"=0 , "AAC"=1 , "AAG"=2, "AAT"=3, ... "ACA"=16
    ************************************************************************ */
    virtual int getSymbolNumber( const string & codonSymbol, unsigned int symbolCategory = 0 ) const;
    int getSingleSymbolNumber(const char & base) const;

    /** ************************************************************************
     * getSymbol
     * @input   symbolNumber, the index of a symbol of the model
     * @return  The string which represents the state
     ************************************************************************ */
    virtual string getSymbol( unsigned int symbolNumber, unsigned int symbolCategory = 0 ) const;
    string getSingleSymbol( unsigned int symbolNumber) const;
    unsigned int getMask( unsigned int symbolNumber ) const;
};

#endif //CODONPARENT_H
