#include "Tree/SearchTreeBasic.h"

#include "Util/ParametersSet.h"
#include "Util/FileParser.h"
#include "Util/randombox.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Tree/OptimizerNode.h"
#include "Tree/ClustersTreeNode.h"

#include "Models/Model.h"

#include <iostream>
#include <iomanip>
#include <assert.h>


SearchTreeBasic::SearchTreeBasic(const string & registrationName) : OptimizerTreeBasic(){
    Singleton < Factory<SearchTree> > & treeFactory = Singleton < Factory<SearchTree> >::instance();
    treeFactory.subscribe( this, registrationName );
}

SearchTreeBasic::SearchTreeBasic(ParametersSet& parameters):OptimizerTreeBasic(parameters){
}

bool SearchTreeBasic::initialisation( SequenceTable * ptable, Model * pmodel ){
    OptimizerTreeBasic::loadDataAndModel( ptable, pmodel );
    //the optimizeModel flag should have been set up already
    optimizeModelPenalty = (pmodel->getNumberPenaltyParameters()!=0);
    if ( !optimizeModelPenalty && optimizeModel ){
        optimizeModelPenalty = (pmodel->getLnPenalty() != 0.0);
    }
    return false;
}

SearchTreeBasic::~SearchTreeBasic(){
}

void SearchTreeBasic::diffloglikelihood(vector < double > & gradientVector){
    OptimizerTreeBasic::diffloglikelihood( gradientVector );
}

