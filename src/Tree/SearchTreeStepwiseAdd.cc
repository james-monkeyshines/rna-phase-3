#include "Tree/SearchTreeStepwiseAdd.h"

#include "Util/ParametersSet.h"
#include "Util/Optimise.h"

#include "Models/Model.h"

#include "Tree/OptimizerNode.h"

#include <float.h>
#include <iostream>
#include <iomanip>

using namespace std;

//register the tree to the tree factory with the name Stepwise addition
SearchTreeStepwiseAdd SearchTreeStepwiseAdd::prototype("Stepwise addition");


SearchTreeStepwiseAdd::SearchTreeStepwiseAdd( const string & registrationName ) :
SearchTreeInsertionBasic(registrationName){
}

SearchTreeStepwiseAdd::SearchTreeStepwiseAdd(ParametersSet& parameters) :
InferenceTree(),OptimizerTree(),SearchTreeInsertionBasic(parameters){
}

SearchTreeStepwiseAdd::~SearchTreeStepwiseAdd(){
}

bool SearchTreeStepwiseAdd::initialisation( SequenceTable * ptable, Model * pmodel ){
    SearchTreeInsertionBasic::initialisation( ptable, pmodel );

    //random order for the addition of species to the tree
    additionOrder.resize(ptable->getNumberSpecies());
    for ( unsigned int i = 0 ; i < additionOrder.size() ; ++i ) {
        additionOrder[i] = i ;
    }
    additionOrder = statlib::random_permutation( additionOrder );


    //first species
    createInitialThreeLeavesTree( additionOrder[0], additionOrder[1], additionOrder[2] );
    index = 3;
    Optimise::optimiseQuasiNewton( *this, *pmodel, optimizeModel, empiricalFreqs, optimizeModelPenalty ? 1e-9 : 1e-8, false, false );
    cout << endl;
    //return true since we have an initial state with three species
    return true;
}



OptimizerTree* SearchTreeStepwiseAdd::clone(ParametersSet& parameters) const{
    return new SearchTreeStepwiseAdd(parameters);
}

bool SearchTreeStepwiseAdd::processNext(){

    if (index == additionOrder.size()){
       return false;
    }

    OptimizerNode* node = new OptimizerNode(ptable->species[additionOrder[index]], this );
    OptimizerNode* newParent = new OptimizerNode( this );

    cout << "Adding: " << ptable->species[additionOrder[index]] << endl;

    //find possible insertion points
    vector< InferenceNode* > allowedBranches;
    //find possible insertion points
    retrievePossibleBranches( allowedBranches, additionOrder[index] );
    assert(allowedBranches.size());

    unsigned int chosenIndice = 0;
    if (allowedBranches.size()>1){
        cout << "Trying: " << allowedBranches.size() << " tree(s)" << endl;
        double bestLikelihood = -DBL_MAX/3.0;
        OptimizerNode* bestNode = NULL;
        vector< double > oldBranches;
        vector< double > bestBranches;
        vector< double > bestModelParameters;
        vector< double > oldModelParameters;
        vector< double > bestModelPenalty;
        vector< double > oldModelPenalty;
        double bestPenalty = 0.0;
        if (optimizeModel){
            pmodel->getAllParameters( oldModelParameters );
        }
        if (optimizeModelPenalty){
            bestPenalty = -DBL_MAX/3.0;
            pmodel->getAllPenaltyParameters( oldModelPenalty );
        }
        else{
            bestPenalty = pmodel->getLnPenalty();
        }
        getBranchVector( oldBranches );
        for (unsigned int i = 0; i< allowedBranches.size(); ++i){
            //cosmetic separation
            if (i!=0) cout << "--------" << endl;
            
            insertion( node, newParent, (OptimizerNode*)(allowedBranches[i]) );
            cout << i + 1 << '.' << toString(false) << endl;

            Optimise::optimiseQuasiNewton( *this, *pmodel, optimizeModel, empiricalFreqs, optimizeModelPenalty ? 1e-9 : 1e-8, false, false  );
            double lik = loglikelihood();
            double pen = pmodel->getLnPenalty();
            if(optimizeModelPenalty){
                cout << "optimized value=" << lik + pen << "     (";
            }
            cout << "maxLikelihood=" << lik;
            if(optimizeModelPenalty){
                cout << ", maxPenalty=" << pen << ')';
            }
            else{
                if (bestPenalty){
                    cout << "      (penalty=" << bestPenalty << ')';
                }
            }
            cout << endl;
            if (lik+pen>bestLikelihood+bestPenalty){
                chosenIndice = i;
                bestLikelihood = lik;
                bestPenalty = pen;
                bestNode = (OptimizerNode*)allowedBranches[i];
                getBranchVector( bestBranches );
                if (optimizeModel){
                    pmodel->getAllParameters( bestModelParameters );
                }
                if (optimizeModelPenalty){
                    pmodel->getAllPenaltyParameters( bestModelPenalty );
                }
            }

            //roll back for the next try
            removal( node, newParent, (OptimizerNode*)(allowedBranches[i]) );
            setBranchVector( oldBranches );
            if (optimizeModel){
                pmodel->setAllParameters( oldModelParameters );
            }
            if (optimizeModelPenalty){
                pmodel->setAllPenaltyParameters ( oldModelPenalty );
            }
            if (optimizeModel || optimizeModelPenalty){
                pmodel->validChange();
            }
        }
        assert(bestNode);
        insertion(node, newParent, bestNode);
        setBranchVector( bestBranches );
        insertionValidate( node, newParent, bestNode );
        if (optimizeModel){
            pmodel->setAllParameters( bestModelParameters );
        }
        if (optimizeModelPenalty){
            pmodel->setAllPenaltyParameters( bestModelPenalty );
        }
        if (optimizeModel || optimizeModelPenalty){
            pmodel->validChange();
        }
    }
    //if only one possibility
    else{
        insertion(node, newParent, (OptimizerNode*)allowedBranches[0]);
        cout << toString(true) << endl;
        Optimise::optimiseQuasiNewton( *this, *pmodel, optimizeModel, empiricalFreqs, optimizeModelPenalty ? 1e-9 : 1e-8, false, false );
        insertionValidate( node, newParent, (OptimizerNode*)(allowedBranches[0]) );
        double lik = loglikelihood();
        double pen = pmodel->getLnPenalty();
        if(optimizeModelPenalty){
            cout << "optimized value=" << lik + pen << "     (";
        }
        cout << "maxLikelihood=" << lik;
        if(optimizeModelPenalty){
            cout << ", maxPenalty=" << pen << ')';
        }
        else{
            if (pen){
                cout << "      (penalty=" << pen << ')';
            }
        }
        cout << endl;
    }
    ++index;
    return true;
}



void SearchTreeStepwiseAdd::printResults( ostream& outputStream ){
    assert(index == additionOrder.size());
    outputStream << "Model found:" << endl;
    pmodel->printParameters( outputStream );
    
    outputStream << "Tree found:" << endl;
    outputStream << toString() << ';' << endl;

    outputStream.setf(ios::fixed);
    outputStream << setprecision(4);
    double lik = loglikelihood();
    double pen = 0.0;
    if (optimizeModelPenalty){
        pen = pmodel->getLnPenalty();
    }
    if(optimizeModelPenalty){
        outputStream << "optimized value=" << lik + pen << "     (";
    }
    outputStream << "maxLikelihood=" << lik;
    if(optimizeModelPenalty){
        outputStream << ", maxPenalty=" << pen << ')';
    }
    else{
        if (pen){
            outputStream << "      (penalty=" << pen << ')';
        }
    }
    outputStream << endl;
}

void SearchTreeStepwiseAdd::printTree( ostream& outputStream ){
    // Save a tree to a file on its own
    outputStream << toString(false) << ';' << endl;
}
