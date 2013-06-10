#ifndef SEQUENCETABLE_H
#define SEQUENCETABLE_H


#include <stdlib.h>
#include <string.h>
#include <vector>
#include <fstream>
#include <sstream>
#include "Util/array2D.h"

using namespace std;

class Model;
class ParametersSet;


// Sequence table class
// This contains the molecular data for species indexed by a species name
class SequenceTable {
    typedef enum {
        NONE,
        AUTO,
        CATEGORY_LINE,
        ATTR_NUMBER
    } UserCategoryAttribution;

public:
    /** .***********************************************************************
     * SequenceTable
     * default constuctor, contruct an empty SequenceTable
     ************************************************************************ */
    SequenceTable();

    /** ************************************************************************
     * SequenceTable
     * constructs a sequence table from data in file with name "filename"
     ************************************************************************ */
    SequenceTable( ParametersSet& sequenceParameters );

    /** ************************************************************************
     * ~SequenceTable
     * The destructor, free the memory used by dynamic data
     ************************************************************************ */
    ~SequenceTable();

    /** ************************************************************************
     * normalise
     * Normalises the data in the table by removing columns which contain
     * characters not specified in the "valids" array
     * NB : if there is an invalid caracter, assert for the moment....
     ************************************************************************ */
    //void normalise( char * valids, int number );

    /** ************************************************************************
     * findsequence
     * @input   label, the name of the species
     * @return  The index of the sequence of a given species,
     *          -1 if name does not exist
     ************************************************************************ */
    int findSequence( string label );

    /** ************************************************************************
     * getRawSequences
     * @input   speciesIndex, the index of a species
     * @input   symbolIndex, the index of a symbol in rawSequences
     * @pre     (speciesIndex >= 0) && (speciesIndex <= numberSpecies)
     * @pre     (symbolIndex  >= 0) && (symbolIndex <= rawLength)
     * @return  the symbolIndex th string of rawSequences for the
     *          speciesIndex th species
     ************************************************************************ */
    const string & getRawSequences( int speciesIndex,
    int symbolIndex ) {
        return rawSequences(speciesIndex,symbolIndex);
    }

    /** ************************************************************************
     * getInvariantBases
     * @input   the category ID used to treat the desired sequences
     * @return  the invariantBases (vector of pair<string,int>) in the
     *          sequences treated zith this model
     ************************************************************************ */
    inline const vector< pair< string, unsigned int > > & getInvariantBases( int categoryId ){
        return invariantBases[categoryId];
    }

    /** ************************************************************************
     * getNumberSpecies
     * @return  the number of species (=number sequences) in the table
     ************************************************************************ */
    inline unsigned int getNumberSpecies() {
        return species.size();
    }

    /** ************************************************************************
     * getInitialLength
     * @return  the length of the sequences read in the data file
     ************************************************************************ */
    inline unsigned int getInitialLength() {
        return initialSequences.numberColumns();
    }

    /** ************************************************************************
     * getRawLength
     * @return  the length of the sequences of symbol in rawSequences
     ************************************************************************ */
    inline unsigned int getRawLength() {
        return rawSequences.numberColumns();
    }

    /** ************************************************************************
     * getSequencesLength
     * @input   the category ID used to treat to desired sequences
     * @return  the length of the sequences treated with this category
     ************************************************************************ */
    inline unsigned int getSequencesLength( int categoryId ) {
        return sequences[categoryId].numberColumns();
    }

    /** ************************************************************************
     * getInvariantsLength
     * @input   the category ID used to treat to desired sequences
     * @return  the number of different invariant symbols for the given ID
     ************************************************************************ */
    inline unsigned int getInvariantsLength( int categoryId ) {
        return invariantBases[categoryId].size();
    }

    /** ************************************************************************
     * getStructure
     * @return  the structure as it was read from the file
     ************************************************************************ */
    inline const vector<char> & getStructure(){
        return structure;
    }

    /** ************************************************************************
     * getInitialSequences
     * @input   speciesIndex, the index of a species
     * @input   baseIndex, the index of the base
     * @pre     (speciesIndex >= 0) && (speciesIndex <= numberSpecies)
     * @return  the baseIndex th nucleotids for the speciesIndex th specie
     *          as they were read from the file
     ************************************************************************ */
    inline char getInitialSequences( int speciesIndex, int baseIndex ) {
        return initialSequences(speciesIndex, baseIndex );
    }
    
    const vector<int> & getInitialCategoryAttrib(){
        return initialCategoryAttrib;
    }

    /** ************************************************************************
     * getSequences
     * @input   the category ID used to treat to desired sequences
     * @return  the sequences treated with this categrory
     ************************************************************************ */
    inline const array2D< string > & getSequences( int categoryId ) {
        return sequences[categoryId];
    }

    void print( ostream & output );

    /* The type of the sequences contained in the file     */
    /* ("DNA", "RNA" or "MIXED")                           */
    string type;

    /* The species labels */
    vector<string> species;

    /** ************************************************************************
     * retrieveInitialNucleotides
     * retrieve the nucleotides responsible for the apparition of the given
     * symbol with the given category
     * @return    a vector of vector of nucleotide position :
     *            a symbol can appear in the final sequence because of multiple
     *            symbols in the raw sequences, the elements of the first vector
     *            correspond anto these symbols. Each symbol is made with
     *            one or more nucleotides
     ************************************************************************ */
    vector < vector< unsigned int > > retrieveInitialNucleotides( int categoryId, int symbolIndex );
    
    /** ************************************************************************
     * retrieveRawNucleotides
     * retrieve the raw responsible for the apparition of the given
     * symbol with the given category
     * @return    a vector of nucleotide position :
     *            a symbol can appear in the final sequence because of multiple
     *            symbols in the raw sequences, the returned vector
     *            correspond to these symbols.
     ************************************************************************ */
    vector< unsigned int > retrieveRawNucleotides( int categoryId, int symbolIndex );
    
    /** ************************************************************************
     * getNumberCategories
     * @semantics return the number of categories according to the data file
     ************************************************************************ */
     inline unsigned int getNumberCategories(){
         return numberCategories;
     }
    
    void printStats( ostream & output, vector<int> species, vector<int> sites );

protected:

    /** ************************************************************************
     * the initialisation function
     ************************************************************************ */
    void init( ifstream & ff, ParametersSet& sequenceParameters );

    /** ************************************************************************
     * constructInterleaved
     * @input     input, a stringstream where the data are
     * @input     userCategoryAttribution, true if the user has defined the
     *            category attribution
     * @input     userInterleavedStructure, true if the user has defined the 
     *            an interleaved structure line
     * @semantics Constructor primitive for interleaved data, construct 
     *            initialSequences from the data in input
     ************************************************************************ */
    void constructInterleaved(stringstream & input, UserCategoryAttribution userCategoryAttribution, bool userInterleavedStructure);
    
    /** ************************************************************************
     * constructNonInterleaved
     * @input     input : a stringstream where the data are
     * @semantics Constructor primitive for non-interleaved data, construct 
     *            initialSequences from the data in input
     ************************************************************************ */
    void constructNonInterleaved( stringstream & input );

    /** ************************************************************************
     * processStructure
     * Compute the list of paired position of the base-paired nucleotide for a
     * "RNA"-type file
     * @semantics Verify that the pairing mask is complete (i.e.
     *            brackets match)
     *            and build the 'rnaPairs' list of paired position accordingly.
     *            ignore characters different from '(' and ')'.
     ************************************************************************ */
    void processStructure();

    /** ************************************************************************
     * categoryAttribution
     * @input      userCategoryAttribution, true if categories are assigned by the user
     * @input      input, stringstream used to read the user attribution
     * @semantics  Verify that the attribution is consistent with the structure.
     *             if userCategoryAttribution = false, single symbol are assigned
     *             to the first model and paired symbol to the second
     ************************************************************************ */
    void categoryAttribution( UserCategoryAttribution userCategoryAttribution , stringstream& input );

    /** ************************************************************************
     * checkCategoryAgainstStructure
     * @semantics Check the category definition against the structure defintion
     ************************************************************************ */    
    void checkCategoryAgainstStructure();

    /** ************************************************************************
     * constructRawSequences
     * Construct the raw sequences of symbols from the initialSequences of
     * character according to the structure provided (if any)
     * @semantics combine characters into symbol according to the structure
     ************************************************************************ */
    void constructRawSequences();

    
    /** ************************************************************************
     * constructSequences
     * Construct the sequences of symbols from the rawSequences of
     * symbol
     * @semantics split symbols according to their model and extract invariants
     ************************************************************************ */
    void constructSequences();
    
    /** ************************************************************************
     * putCharStructure
     * @input      ch, character to put in the structure
     * @semantics  Put the character in the the structure
     ************************************************************************ */    
    void putCharStructure(char ch);

    /* Protected data structures                                              */

    // The sequences contained in the data file
    array2D< char > initialSequences;

    // vector used to store the RNA secondary structure
    vector<char> structure;
    // list of paired position (RNA)
    vector < pair < int, int > > pairs;
    // list of the codon first nucleotide positions
    vector < int > triplets;

    // list of user model attribution for each base (must be consistent with
    // the structure)
    vector<int> initialCategoryAttrib;

    // Store the model used for each base (or paired base)
    // the size of modelAttrib is the size once the structure
    // has been "applied" (rawSequences size)
    vector<int> categoryAttrib;

    unsigned int numberCategories;

    // The raw sequences read from the file
    // rawSequences(i,j) stores the partial sequence of the ith species
    // treated with the jth model
    array2D< string > rawSequences;
    
    //link the initial sequences with the rawSequences
    vector<unsigned int> initial2Raw;
    vector< vector<unsigned int> > raw2Initial;

    // sequence after processing : sequences[modelId](speciesId,symbolId)
    // invariantBases are removed
    vector< array2D< string > > sequences;

    //link the raw sequences with the final sequences
    vector< pair<unsigned int,int> > raw2Final;
    //final2Raw[modelId][finalId]
    vector < vector < unsigned int > > final2Raw;
    //invariant2Raw[modelId][invariantId][eltId]
    vector < vector < vector < unsigned int > > > invariant2Raw;
    
    // store the invariant bases (shortcut for the likelihood computation)
    // for each model : invariantBases[modelId][pairId]
    vector < vector< pair<string,unsigned int> > > invariantBases;
    
};
#endif




