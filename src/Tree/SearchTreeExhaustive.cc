#include "Tree/SearchTreeExhaustive.h"

#include "Util/ParametersSet.h"
#include "Util/Optimise.h"
#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/ParametersSet.h"
#include "Util/Optimise.h"

#include "Models/Model.h"

#include "Tree/OptimizerNode.h"

#include <float.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the tree to the tree factory with the name Exhaustive
SearchTreeExhaustive SearchTreeExhaustive::prototype("Exhaustive");


SearchTreeExhaustive::SearchTreeExhaustive( const string & registrationName ) :
SearchTreeInsertionBasic(registrationName){
}

SearchTreeExhaustive::SearchTreeExhaustive(ParametersSet& parameters) :
InferenceTree(),OptimizerTree(),SearchTreeInsertionBasic(parameters){
}

SearchTreeExhaustive::~SearchTreeExhaustive(){
}

bool SearchTreeExhaustive::initialisation( SequenceTable * ptable, Model * pmodel ){
    SearchTreeInsertionBasic::initialisation( ptable, pmodel );

    if(optimizeModelPenalty){
        pmodel->getAllPenaltyParameters( initialModelPenalty );
    }
    if(optimizeModel){
        pmodel->getAllParameters( initialModelParameters );
    }
    
    //random order for the addition of species to the tree
    additionOrder.resize(ptable->getNumberSpecies());
    for ( unsigned int i = 0 ; i < additionOrder.size() ; ++i ) {
        additionOrder[i] = i ;
    }
    additionOrder = statlib::random_permutation( additionOrder );

    waitingClusters.clear();

    //first species...
    createInitialThreeLeavesTree( additionOrder[0], additionOrder[1], additionOrder[2] );

    //create all the nodes
    for (unsigned int i = 3; i < ptable->getNumberSpecies(); ++i){
       //create the node for the species
        OptimizerNode * node = new OptimizerNode( ptable->species[additionOrder[i]], this );
        //and an insertion node
        OptimizerNode * newParent = new OptimizerNode( this );
        insertNodes.push_back( pair< OptimizerNode*, OptimizerNode* > (newParent,node) );
    }

    bestLnLik = -DBL_MAX;
    bestLnPenalty = -DBL_MAX;
    //don't display this initial state
    return false;
}



OptimizerTree* SearchTreeExhaustive::clone(ParametersSet& parameters) const{
    return new SearchTreeExhaustive(parameters);
}

bool SearchTreeExhaustive::backTrack(){

    OptimizerNode* oldInsert = *(branchIter.top());
    //remove the last added node (removalValidate and removal are always called together in this class)
    removalValidate( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                           oldInsert );
    removal( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                           oldInsert );
    ++(branchIter.top());

    //and continue to backtrack until we reach a subtree that leads to a new unexplored portion of the topology space
    while(branchIter.top()==branchesStack.top().end()){
        //remove node, reassign clusters
        branchIter.pop();
        branchesStack.pop();
        //finished?
        if (branchIter.empty()){
            return false;
        }
        //the node is removed definitely (removalValidate and removal are always called together in this class)
        removalValidate( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                                *(branchIter.top()) );
        removal( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                                *(branchIter.top()) );
        ++(branchIter.top());
    }
    //try next, restore the state as it was
    //insert the next node (id=insertionPoints.size()) to the next possible branch
    //(insertion and insertionValidate are always called together in this class)
    insertion( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                        *(branchIter.top()) );
    insertionValidate( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                       *(branchIter.top()) );
    return true;
}

bool SearchTreeExhaustive::buildRec(){

     //candidate tree found? ie, no more species to add?
    if (branchesStack.size()+3 == ptable->getNumberSpecies()){
        if(optimizeModelPenalty){
            pmodel->setAllPenaltyParameters( initialModelPenalty );
        }
        if(optimizeModel){
            pmodel->setAllParameters( initialModelParameters );
        }
        if (optimizeModel || optimizeModelPenalty){
            pmodel->validChange();
        }
        //call the unique instance of randombox and randomize branch lenghts
        Singleton< randombox > & randBox = Singleton< randombox >::instance();
        vector< double > bv;
        getBranchVector(bv);
        for( vector<double>::iterator iter = bv.begin(); iter != bv.end(); ++iter ){
            *iter = randBox.ran()*.6;
        }
        Optimise::optimiseQuasiNewton( *this, *pmodel, optimizeModel, empiricalFreqs, optimizeModelPenalty ? 1e-9 : 1e-8, false, false );
//        cout << endl << toString() << ';' << endl;
        double lik = loglikelihood();
        double pen = pmodel->getLnPenalty();
//        if(optimizeModelPenalty){
//            cout << "optimized value=" << lik + pen << "     (";
//        }
//        cout << "maxLikelihood = " << lik;
//        if(optimizeModelPenalty){
//            cout << ", maxPenalty=" << pen << ')';
//        }
//        else{
//            if (pen){
//                cout << "      (penalty=" << bestLnPenalty << ')';
//            }
//        }
//        cout << endl << endl;

        if (lik+pen>bestLnLik+bestLnPenalty){
            bestLnLik = lik;
            bestLnPenalty = pen;
            if (optimizeModel){
                pmodel->getAllParameters(  bestModelParameters );
            }
            if (pmodel->getNumberPenaltyParameters()){
                pmodel->getAllPenaltyParameters( bestModelPenalty );
            }
            bestTree = toString();
        }
        
        //better or not, return the result (the search will restart by a call to processNext)
        return true;
    }
    //usual case, the tree is still a partial tree
    else{
        //add a layer of possible insertion point for the next species:
        //we find where it is possible to insert the next species and update clusters
        vector< InferenceNode* > allowedBranches;
        //find possible insertion points
        retrievePossibleBranches( allowedBranches, additionOrder[branchesStack.size()+3] );
        assert(allowedBranches.size());

        //push the possible insertion points in possibleBranches
        branchesStack.push( vector<OptimizerNode*> ());
        for ( vector< InferenceNode* >::const_iterator iter = allowedBranches.begin();
              iter != allowedBranches.end(); ++iter ){
            branchesStack.top().push_back( (OptimizerNode*)*iter );
        }
        branchIter.push( branchesStack.top().begin() );
        
        //insert the next node (id=insertionPoints.size()) to the next possible branch
        //(insertion and insertionValidate are always called together in this class)
        insertion( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                            *(branchIter.top()) );
        insertionValidate( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                            *(branchIter.top()) );
                            
        return buildRec();
    }
}

bool SearchTreeExhaustive::processNext(){
    //the algorithm:
    //The tree is constructed by successive addition of species and backtracking.
    //Each topology is visited.
    //processNext returns false once the topology space was covered


    bool contin = true;

    //each time processNext is called a backtrack is needed.
    //(each time but the first one)
    if ( !branchIter.empty() ){
        //backTrack destroys the current tree (it should be a complete one) and
        //return to the next state we have to try.
        contin = backTrack();
    }
    //buildRec continues the building procedure at the point designed by the current stacks if it is not finished
    if (contin==true){
        contin=buildRec();
    }
    return contin;
}

void SearchTreeExhaustive::printResults( ostream& outputStream ){
    //this one should be useless
    if ( pmodel->getNumberPenaltyParameters() ){
        pmodel->setAllPenaltyParameters( bestModelPenalty );
    }
    if (optimizeModel){
        outputStream << "Best model found:" << endl;
        pmodel->setAllParameters(  bestModelParameters );
        pmodel->validChange();
        pmodel->printParameters( outputStream );
    }
    outputStream << "Best Tree found:" << endl;
    outputStream << bestTree << ';' << endl;
    
    outputStream.setf(ios::fixed);
    outputStream << setprecision(4);
    if(optimizeModelPenalty){
        outputStream << "optimized value=" << bestLnLik + bestLnPenalty << "     (";
    }
    outputStream << "maxLikelihood=" << bestLnLik;
    if(optimizeModelPenalty){
        outputStream << ", maxPenalty=" << bestLnPenalty << ')';
    }
    else{
        if (bestLnPenalty){
            outputStream << "      (penalty=" << bestLnPenalty << ')';
        }
    }
    outputStream << endl;
}

void SearchTreeExhaustive::printTree( ostream& outputStream ){
    // Save a tree to a file on its own
    outputStream << bestTree << ';' << endl;
}
