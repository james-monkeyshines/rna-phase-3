#include "Tree/SimulationNode.h"

#include "Tree/SimulationTree.h"

#include "Tree/TreeMap.h"

#include "Util/array2D.h"
#include "Util/array3D.h"

#include "Models/Model.h"

SimulationNode::SimulationNode( SimulationTree* ptree ):
BasicNode(){
    //initialise the inference node like the basic node
    this->ptree = ptree;
}

SimulationNode::SimulationNode( const string& label, SimulationTree* ptree ):
BasicNode( label ) {
    this->ptree = ptree;
}

SimulationNode::SimulationNode( const TreeMap& treeMap, SimulationTree* ptree ){
    this->ptree = ptree;
    //if this node is an internal node
    if ( treeMap.getNumberChildren() ){
        parent = NULL;
        for ( unsigned int i = 0; i < treeMap.getNumberChildren() ; ++i ) {
            addChild( new SimulationNode( treeMap.getChildMap(i), ptree ),
                      treeMap.getDistance(i) );
        }
        updateNodesCount();
    }
    //if this node is a leaf
    else{
        parent = NULL;
        nbInternalNodes = 0;
        nbLeaves = 1;
        this->label = treeMap.getLabel();
    }
}

void SimulationNode::simulate( const vector< vector<unsigned int> > & categories,
                          const vector< vector<unsigned int> > & parent ) {

    Singleton < randombox > & randBox = Singleton < randombox >::instance();

    unsigned int numberModels = categories.size();
            
    sequence.resize(numberModels);
    
    //if root node
    if (this==ptree->getRoot()){
        array2D<double> freqVector;
        for ( unsigned int modelId = 0; modelId < numberModels; ++modelId ) {
            sequence[modelId].reserve(categories[modelId].size());
            freqVector.resize( ptree->pmodel->getNumberRatesCategories(modelId),
                               ptree->pmodel->getNumberStates(modelId) );
            for ( unsigned int rate = 0; rate < ptree->pmodel->getNumberRatesCategories(modelId); ++rate ) {
                double cumul = 0.0;
                for ( unsigned int state = 0; state < ptree->pmodel->getNumberStates(modelId)-1; ++state ) {
                    cumul += ptree->pmodel->getFrequency( state, rate, modelId );
                    freqVector(rate,state) = cumul;
                }
                freqVector(rate, ptree->pmodel->getNumberStates(modelId)-1) = 1.0;
            }
            for ( unsigned int site = 0; site < categories[modelId].size(); ++site ) {
               /* assign a frequency according to the probabilities for the rate */
                unsigned int rate = categories[modelId][site];
                unsigned int state = 0;
                double number = randBox.ran();
                while ( number > freqVector(rate,state) ) ++state;
                assert ( state < ptree->pmodel->getNumberStates(modelId) );
                sequence[modelId].push_back(state);
            }
        }
    }
    //else, normal internal node or leaf
    else{
        array3D<double> lookupTable;
        for ( unsigned int modelId = 0; modelId < numberModels; ++modelId ) {
            sequence[modelId].reserve(categories[modelId].size());
            lookupTable.resize( ptree->pmodel->getNumberRatesCategories(modelId),
                                ptree->pmodel->getNumberStates(modelId),
                                ptree->pmodel->getNumberStates(modelId) );
            for ( unsigned int rate = 0; rate < ptree->pmodel->getNumberRatesCategories(modelId); ++rate ) {
                for ( unsigned int parentState = 0; parentState < ptree->pmodel->getNumberStates(modelId); ++parentState ) {
                    double cumul = 0.0;
                    for ( unsigned int state = 0; state < ptree->pmodel->getNumberStates(modelId)-1; ++state ) {
                        cumul += ptree->pmodel->probability( parentState, state, pdistance, rate, modelId );
                        lookupTable(rate, parentState, state) = cumul;
                    }
                    lookupTable(rate, parentState, ptree->pmodel->getNumberStates(modelId)-1) = 1.0;
                }
            }
            for ( unsigned int site = 0; site < categories[modelId].size(); ++site ) {
                unsigned int state = 0;
                unsigned int rate = categories[modelId][site];
                double number = randBox.ran();
                while ( number > lookupTable(rate, parent[modelId][site], state ) ) ++state;
                assert ( state < ptree->pmodel->getNumberStates(modelId) );
                assert ( state < ptree->pmodel->getNumberStates(modelId) );
                sequence[modelId].push_back(state);
            }
        }
    }
    if ( !isLeaf() ){
        for ( list<BasicNode*>::iterator iter = getChildrenList().begin();
              iter != getChildrenList().end(); ++iter ){
            ((SimulationNode*)(*iter))->simulate( categories, sequence );
        }
    }
}
