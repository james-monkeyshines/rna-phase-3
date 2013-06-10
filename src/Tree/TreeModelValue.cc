#include "Tree/TreeModelValue.h"

#include <assert.h>

int TreeModelValue::compareEntries(const void *v1 , const void * v2){
    assert(v1 != NULL);
    assert(v2 != NULL);
    TreeModelValue * e1 = *((TreeModelValue **) v1) ;
    TreeModelValue * e2 = *((TreeModelValue **) v2) ;


    if (e1->logLikelihood == e2->logLikelihood)
        return 0 ;
    if (e1->logLikelihood > e2->logLikelihood)
        return -1 ;
    return 1 ;
}


// String rep. of tree , Model parameter vector and ML value
TreeModelValue::TreeModelValue() : stringTree(""){
  logLikelihood = 0.0 ;
}

TreeModelValue::TreeModelValue(const string& stringTree,
                               const vector <double> & branchLengths,
                               const vector <double> & modelParameters,
                               double logLikelihood) {
  this->stringTree = stringTree ;
  this->branchLengths = branchLengths ;
  this->modelParameters = modelParameters ;
  this->logLikelihood = logLikelihood ;
}





