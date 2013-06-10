#include "Tree/MLOptimizerTree.h"

#include "Util/ParametersSet.h"
#include "Util/FileParser.h"
#include "Util/randombox.h"

#include "PatternDesign/Singleton.h"

#include <iostream>
#include <iomanip>
#include <assert.h>

//register the tree to the tree factory with the name ML optimizer tree
MLOptimizerTree MLOptimizerTree::prototype("ML optimizer tree");


MLOptimizerTree::MLOptimizerTree( const string & registrationName ) :
InferenceTree(), OptimizerTreeBasic(registrationName){
}

MLOptimizerTree::MLOptimizerTree(ParametersSet& parameters) :
InferenceTree(),OptimizerTreeBasic(parameters){
    optimizedTree = parameters.stringParameter( "Tree" );
}

MLOptimizerTree::~MLOptimizerTree(){
}


bool MLOptimizerTree::initialisation( SequenceTable * ptable, Model * pmodel ){
    loadDataAndModel( ptable, pmodel );
    constructFromString( optimizedTree );
    vector<double> branchLengths;
    Singleton<randombox>& randBox = Singleton<randombox>::instance();
    vector<BasicNode*>::iterator iter = nodeRefVector.begin();
    ++iter;
    while (iter != nodeRefVector.end()){
        if ((*iter)->getParentDistance()<=0.0){
            (*iter)->setParentDistance( 0.7 * randBox.ran() );
        }
        ++iter;
    }
    return true;
}


OptimizerTree* MLOptimizerTree::clone(ParametersSet& parameters) const{
    return new MLOptimizerTree(parameters);
}
