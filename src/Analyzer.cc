#include "Analyzer.h"

#include "PatternDesign/Singleton.h"
#include "Util/randombox.h"

#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include "Sequence/SequenceTable.h"

#include <iostream>
#include <iomanip>
#include <map>
#include <algorithm>
#include <functional>
#include <ctype.h>

using namespace std;

int Analyzer::run( int argc, char * argv[] ){
    Phase::run( argc, argv );
    
    // Retrieve the control-file parameters
    ParametersSet parameters("DATAFILE");
    
    if ( argc > 2 ){
        cerr << "usage : analyser data_file" << endl;
    }
    if ( argc == 1 ){
        cout << "Data file ? ";
        cin >> parameters["Data file"];
    }
    else{
        parameters["Data file"] = argv[1];
    }
    cout << "Interleaved data file? ";
    cin >> parameters["Interleaved data file"];
    
    cout << "Heterogeneous data models? ";
    cin >> parameters["Heterogeneous data models"];
    
    if (parameters.boolParameter("Interleaved data file")){
        cout << "Interleaved structure? ";
        cin >> parameters["Interleaved structure"];
    }
    
    /** ************************************************************************
     * Data File Initialization
     ************************************************************************ */
    // Read in molecular data and remove unknown characters from the data matrix
    SequenceTable * table = new SequenceTable( parameters );
        
    unsigned int cat = table->getNumberCategories();
    //if more than one model
    if (cat>1){
        cout << "class to check (1-" << table->getNumberCategories() << ") ? ";
        cin >> cat;
    }
    --cat;
    
    string lumpFile;
    cout << "name of your .lmp file (e.g., input-data/dna.lmp or input-data/rna.lmp) ? ";
    cin >> lumpFile;
    
    // Create a FileParser from the control-file
    FileParser * fileParser = new FileParser( lumpFile.c_str() );
    // Retrieve the control-file parameters
    ParametersSet* ret = fileParser->retrieveParametersSet();
    delete fileParser;
    if (!ret){
        cerr << "malformed lump file" << endl;
        exit(EXIT_FAILURE);
    }
    ParametersSet lump = *ret;
    delete ret;

    unsigned int nbSpecies = table->getNumberSpecies();
    
    //assign each symbol an id
    map< string, int > symbol;
    //store the number of each symbol (for this model) found over all the sequences
    vector< int > nbSymb; 
    //the same but keep a separate counter for each species
    vector< int > * nbSymbSpecies = new vector< int >[nbSpecies];
    
    unsigned int sequenceLength = table->getSequencesLength( cat ); //TEMPORARILY, use the number of variables sites for sequenceLength
    
    //parse variable sites first
    int idSymb = 0;
    for ( unsigned int j = 0; j < nbSpecies; ++j ){
        string strSymb;
        for ( unsigned int index = 0; index < sequenceLength; ++index ){
            if  (!lump.findParameter(table->getSequences( cat )(j,index))){
                cerr << "Warning, the symbol \"" << table->getSequences( cat )(j,index)
                     << "\" in your dataset is not recognized in your lump file" << endl;
                strSymb = "";
            }
            else{
                strSymb = lump.stringParameter(table->getSequences( cat )(j,index));
            }
            //find the iterator of the symbol (to obtain its id)
            map< string, int >::iterator iterSymb = symbol.find(strSymb);
            //new symbol ?
            if (iterSymb == symbol.end()){
                // reference the new symbol
                symbol[strSymb] = idSymb;
                //and initialise its counter
                nbSymb.push_back(1);
                // increase the size of the vector for each species (and init at 0)
                for (unsigned int k = 0; k < nbSpecies; ++k){
                    nbSymbSpecies[k].push_back(0);
                }
                //for the species j, 1 element is found
                nbSymbSpecies[j][idSymb] = 1;
                //increase the idSymb for the next new symbol
                ++idSymb;
            }
            //otherwise add 1 to the total number of element and add 1 to
            //the number for this species
            else{
                ++nbSymb[(*iterSymb).second];
                ++nbSymbSpecies[j][(*iterSymb).second];
            }
        }
    }

    //invariant sites now, update sequenceLength on the way
    for ( unsigned int index = 0; index < table->getInvariantsLength( cat );
          ++index ){
        string strSymb;
        if  (!lump.findParameter((table->getInvariantBases( cat )[index]).first)){
           cerr << "Warning, the symbol \"" << (table->getInvariantBases( cat )[index]).first
                << "\" in your dataset is not recognized in your lump file" << endl;
            strSymb = "";
        }
        else{
            strSymb = lump.stringParameter((table->getInvariantBases( cat )[index]).first);
        }
        unsigned int nbSites = (table->getInvariantBases( cat )[index]).second;
        sequenceLength += nbSites;
        //find the iterator of the symbol (to obtain its id)
        map< string, int >::iterator iterSymb = symbol.find(strSymb);
        //new symbol ?
        if (iterSymb == symbol.end()){
            // reference the new symbol
            symbol[strSymb] = idSymb;
            //and initialise its counter (the symbol is invariant nbSites times so it is
            //present nbSites*nbSpecies times in the alignment)
            nbSymb.push_back(nbSites*nbSpecies);
            // increase the size of the vector for each species (initialise at nbSites for each)
            for (unsigned int k = 0; k < nbSpecies; ++k){
                nbSymbSpecies[k].push_back(nbSites);
            }
            //increase the idSymb for the next new symbol
            ++idSymb;
        }
        //otherwise add 1 to the total number of element and add 1 to
        //the number per species
        else{
            nbSymb[(*iterSymb).second] += nbSites*nbSpecies;
            for (unsigned int k = 0; k < nbSpecies; ++k){
                nbSymbSpecies[k][(*iterSymb).second] += nbSites;
            }
        }
    }
    
    // for each symbol (there is idSymb symbol) :
    // 1) compute the mean percentage
    double mean[idSymb];
    for (int i = 0; i < idSymb; ++i){
        mean[i] = (double)nbSymb[i]/(double)( nbSpecies * sequenceLength );
    }

    // for each symbol, browse all the species in order to
    // 2)find the id of the species with the minimum and maximum content
    // 3)store the value of this min/max content
    // 4)compute a standard deviation of the content
    int minSpecies[idSymb];
    double min[idSymb];
    int maxSpecies[idSymb];
    double max[idSymb];
    double sigma[idSymb];
        
    //for all the symbols found,
    for (int i = 0; i < idSymb; ++i ){
        //initialise minmax values
        min[i] = mean[i];
        max[i] = mean[i];
        //initialise the species with the min and max content for the symbol too
        minSpecies[i] = 0;
        maxSpecies[i] = 0;
        //and the SD
        sigma[i] = 0;
        for ( unsigned int j = 0; j < nbSpecies; ++j ){
            double percent = ( (double)nbSymbSpecies[j][i] )/sequenceLength;
            if ( min[i] > percent ){
                min[i] = percent;
                minSpecies[i] = j;
            }
            if ( max[i] < percent ){
                max[i] = percent;
                maxSpecies[i] = j;
            }
            sigma[i] += ( (percent-mean[i])*(percent-mean[i]) );
        }
        sigma[i] = sqrt(sigma[i]/(double)nbSpecies);
    }
    //print the results
    cout.setf(ios::fixed);
    cout << setprecision(2);
    for ( map< string, int >::iterator iter = symbol.begin();
          iter != symbol.end(); ++iter ){
        int id = (*iter).second;
        cout << "Symbol " << (*iter).first << " :" << endl
             << "mean = " << mean[id] * 100.0 << '%' << endl
             << "min = " << min[id]*100.0 << "%, species : "
             << table->species[minSpecies[id]] << endl
             << "max = " << max[id]*100.0 << "%, species : "
             << table->species[maxSpecies[id]] << endl
             << "standard deviation = " << sigma[id]*100.0 << '%' << endl;
    }
    
    double cutoff = -1.0;
    while ( (cutoff<0.0) || (cutoff>1.0) ){
        cout << "cut-off ([0, 1.0])? ";
        cin >> cutoff;
    }
    
    cout << "symbols for cut-off : ";
    for ( ParametersSet::iterator iter = lump("MISMATCH").begin();
          iter != lump("MISMATCH").end(); ++iter ){
        if (iter!=lump("MISMATCH").begin()){
            cout << ", ";
        }
        cout << (*iter).first;
    }
    cout << endl;
    vector< pair<unsigned int, string> > unrecognized;
    for ( int index = -table->getInvariantsLength(cat);
          index < (int)table->getSequencesLength( cat ); ++index ){
        int mismatch = 0;
        unrecognized.clear();
        for ( unsigned int j = 0; j < nbSpecies; ++j ){
            string symbol;
            if (index < 0){
                symbol = table->getInvariantBases( cat )[-(1+index)].first;
            }
            else{
                symbol = table->getSequences( cat )( j, index);
            }
            if ( !lump.findParameter(symbol) ){
                unrecognized.push_back(pair<unsigned int, string>(j,symbol));
            }
            else{
                if ( lump("MISMATCH").findParameter(lump.stringParameter(symbol)) ){
                    ++mismatch;
                }
            }
        }
        double freqMismatch = ((double)mismatch)/((double)nbSpecies);
        if(unrecognized.size()){
            cout << "unrecognized symbol(s): ";
            vector < vector< unsigned int > > pos =
                table->retrieveInitialNucleotides( cat, index );
            for ( vector < vector< unsigned int > >::iterator iter = pos.begin();
                  iter != pos.end(); ++iter ){
                for ( vector< unsigned int >::iterator iter2 = (*iter).begin();
                      iter2 != (*iter).end(); ++iter2 ){
                     if (iter2 != (*iter).begin() ){
                         cout << ", ";
                     }
                     cout << *iter2 + 1;
                }
                cout << "; ";
            }
            cout << "#species: ";
            for ( vector< pair<unsigned int, string> >::iterator iter = unrecognized.begin();
                  iter != unrecognized.end(); ++iter ){
                if (iter != unrecognized.begin() ){
                    cout << ", ";
                }
                cout << '(' << iter->second << ','
                     << table->species[iter->first] << ')';
            }
            cout << endl;
        }
        if ( freqMismatch > cutoff ){
            cout << "mismatch = " << freqMismatch*100.0 << "% : ";
            vector < vector< unsigned int > > pos =
                table->retrieveInitialNucleotides( cat, index );
            for ( vector < vector< unsigned int > >::iterator iter = pos.begin();
                  iter != pos.end(); ++iter ){
                for ( vector< unsigned int >::iterator iter2 = (*iter).begin();
                      iter2 != (*iter).end(); ++iter2 ){
                     if (iter2 != (*iter).begin() ){
                         cout << ", ";
                     }
                     cout << *iter2 + 1;
                }
                cout << "; ";
            }
            cout << endl;
        }
    }
    delete [] nbSymbSpecies;
    return (EXIT_SUCCESS);
}

int main( int argc, char* argv[]){
    Analyzer analyzer;
    int res = analyzer.run( argc, argv );
    return res;
}

