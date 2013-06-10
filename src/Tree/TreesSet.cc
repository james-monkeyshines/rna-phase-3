#include <cstdlib>

#include "Tree/TreesSet.h"

#include "Tree/TreeModelValue.h"

TreesSet::TreesSet(){
    numberEntries = sofar = 0 ;
    entries = NULL ;
    smallestIndex = -1 ;

}

TreesSet::TreesSet(int numberEntries){

    this->numberEntries = numberEntries ;
    entries = new TreeModelValue*[numberEntries] ;

    for ( int i = 0; i < numberEntries ; ++i){
        entries[i] = NULL;
    }
    smallestIndex = -1;
    sofar = 0;
}

TreesSet::~TreesSet() {
    for ( int i=0 ; i < size(); ++i){
        delete entries[i] ;
    }
    delete[] entries ;
}

int TreesSet::size() {
  if (numberEntries <= sofar){
      return numberEntries;
  }
  return (sofar) ;
}

void TreesSet::operator=(const TreesSet& src) {

  for ( int i=0 ; i < size() ; ++i) {
    if (entries[i] != NULL)
      delete entries[i] ;
  }
  
    if (entries != NULL){
        delete[] entries ;
    }
  

    numberEntries = src.numberEntries ;
    entries = new TreeModelValue*[numberEntries] ;
    smallestIndex = src.smallestIndex ;
    sofar = src.sofar ;

  
    for ( int i = 0; i < numberEntries; ++i) {
        if ( src.entries[i] != NULL ){
          entries[i] = new TreeModelValue( src.entries[i]->stringTree,
                                           src.entries[i]->branchLengths,
                                           src.entries[i]->modelParameters,
                                           src.entries[i]->logLikelihood );
        }
    }
}

bool TreesSet::insert(const string& stringTree,
                      const vector <double>& branchLengths,
                      const vector <double>& modelParameters,
                      double logLikelihood) {
    TreeModelValue* newEntry = NULL ; 
  
    if (sofar < numberEntries) {
        newEntry = new TreeModelValue(stringTree, branchLengths,
                                      modelParameters, logLikelihood);
        entries[sofar] = newEntry;
        ++sofar;
    }
    else {
        if (entries[smallestIndex]->logLikelihood < logLikelihood) {
            newEntry = new TreeModelValue(stringTree, branchLengths,
                                          modelParameters, logLikelihood);
            if (entries[smallestIndex] != NULL){
                delete entries[smallestIndex];
            }
            entries[smallestIndex] = newEntry ;
        }
        else{
            return false;
        }
    }
  
    smallestIndex = numberEntries;
    if (sofar<smallestIndex){
        smallestIndex = sofar;
    }
    qsort(entries, smallestIndex, sizeof(newEntry),
                               TreeModelValue::compareEntries);
    --smallestIndex;
    return true ;
}


TreeModelValue * TreesSet::get(int index) {
    return (entries[index-1]) ;
}
