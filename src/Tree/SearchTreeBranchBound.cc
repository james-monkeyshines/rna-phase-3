#include "Tree/SearchTreeBranchBound.h"

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



//Probably one of the hardest (and worse) algorithm in phase
//because of the clusters and because the root is dynamically
//changed during the process
//to check it uncomment following flag, and modify returnValue


#ifdef BRANCHBOUNDDEBUG1
ofstream testVisited;
stack< string > treeSave;

double returnValue(InferenceTree* tree){
    cout << "WARNING! DEBUG MODE" << endl;
    switch( tree->getNumberBranches() ){
        case 9:
            //heuristic result
            if (tree->toString(true) == "(90CF11697,93TH057,(U455,(Q23,(BRU,NDK))))") return -7000.0;
            //same heuristic tree but in a different form since not rerooted yet
            if (tree->toString(true) == "(U455,(90CF11697,93TH057),(Q23,(BRU,NDK)))") return -7000.0;
            //only (90CF11697,NDK,(U455,(BRU,Q23))) and (90CF11697,U455,(Q23,(BRU,NDK)))    will pass the -7000 bound -> 14 trees
            //the 6 suboptimal trees near the heuristic tree
            if (tree->toString(true) == "(90CF11697,U455,(Q23,(NDK,(93TH057,BRU))))") return -7300.0;
            if (tree->toString(true) == "(90CF11697,(93TH057,U455),(Q23,(BRU,NDK)))") return -7300.0;
            if (tree->toString(true) == "(90CF11697,U455,(93TH057,(Q23,(BRU,NDK))))") return -7300.0;
            if (tree->toString(true) == "(90CF11697,U455,((93TH057,Q23),(BRU,NDK)))") return -7600.0;
            if (tree->toString(true) == "(90CF11697,U455,(Q23,(93TH057,(BRU,NDK))))") return -7600.0;
            if (tree->toString(true) == "(90CF11697,U455,(Q23,(BRU,(93TH057,NDK))))") return -7600.0;
            //the optimal tree (not rerooted yet)
            if (tree->toString(true) == "(NDK,(90CF11697,93TH057),(U455,(BRU,Q23)))") return -6000.0;
            //the 6 suboptimal trees near the optimal tree
            if (tree->toString(true) == "(90CF11697,(93TH057,NDK),(U455,(BRU,Q23)))") return -6500.0;
            if (tree->toString(true) == "(90CF11697,NDK,(93TH057,(U455,(BRU,Q23))))") return -6500.0;
            if (tree->toString(true) == "(90CF11697,NDK,((93TH057,U455),(BRU,Q23)))") return -6500.0;
            if (tree->toString(true) == "(90CF11697,NDK,(U455,(93TH057,(BRU,Q23))))") return -7500.0;
            if (tree->toString(true) == "(90CF11697,NDK,(U455,(Q23,(93TH057,BRU))))") return -7500.0;
            if (tree->toString(true) == "(90CF11697,NDK,(U455,(BRU,(93TH057,Q23))))") return -7500.0;

            if (tree->toString(true) == "(U455,(BRU,Q23),(NDK,(90CF11697,93TH057)))") return -6000.0;
            if (tree->toString(true) == "((90CF11697,NDK),(93TH057,U455),(BRU,Q23))") return -6500.0;
            if (tree->toString(true) == "(U455,(BRU,Q23),(93TH057,(90CF11697,NDK)))") return -6500.0;
            if (tree->toString(true) == "(U455,(BRU,Q23),(90CF11697,(93TH057,NDK)))") return -6500.0;
            if (tree->toString(true) == "(U455,(90CF11697,NDK),(93TH057,(BRU,Q23)))") return -7500.0;
            if (tree->toString(true) == "(U455,(90CF11697,NDK),(BRU,(93TH057,Q23)))") return -7500.0;
            if (tree->toString(true) == "(U455,(90CF11697,NDK),(Q23,(93TH057,BRU)))") return -7500.0;

            cerr << "TREE NOT RECOGNIZED DURING DEBUG" << endl;
            return -8000.0;
        break;
        case 7: //only (90CF11697,NDK,(Q23,U455)) and (90CF11697,U455,(NDK,Q23)) -> 10 tree
            //from optimal root
            if (tree->toString(true) == "(90CF11697,NDK,(U455,(BRU,Q23)))") return -5500.0; //opt
            if (tree->toString(true) == "(NDK,(90CF11697,BRU),(Q23,U455))") return -8000.0;
            if (tree->toString(true) == "(90CF11697,(BRU,NDK),(Q23,U455))") return -8000.0;
            if (tree->toString(true) == "(90CF11697,NDK,(Q23,(BRU,U455)))") return -8000.0;
            if (tree->toString(true) == "(90CF11697,NDK,(BRU,(Q23,U455)))") return -8000.0;
            //from heuristic root
            if (tree->toString(true) == "(90CF11697,U455,(Q23,(BRU,NDK)))") return -4500.0; //heur
            if (tree->toString(true) == "(90CF11697,(BRU,U455),(NDK,Q23))") return -8000.0;
            if (tree->toString(true) == "(90CF11697,U455,(BRU,(NDK,Q23)))") return -8000.0;
            if (tree->toString(true) == "(U455,(90CF11697,BRU),(NDK,Q23))") return -8000.0;
            if (tree->toString(true) == "(90CF11697,U455,(NDK,(BRU,Q23)))") return -8000.0;

            if (tree->toString(true) == "(U455,(90CF11697,NDK),(BRU,Q23))") return -5500.0;
            if (tree->toString(true) == "(Q23,(90CF11697,NDK),(BRU,U455))") return -8000.0;
            if (tree->toString(true) == "(Q23,U455,(BRU,(90CF11697,NDK)))") return -8000.0;
            if (tree->toString(true) == "(Q23,U455,(90CF11697,(BRU,NDK)))") return -8000.0;
            if (tree->toString(true) == "(Q23,U455,(NDK,(90CF11697,BRU)))") return -8000.0;

         //   if (tree->toString(true) == "(Q23,U455,(BRU,(90CF11697,NDK)))") return -8000.0;



            cerr << "TREE NOT RECOGNIZED DURING DEBUG" << endl;
            return -8000.0;
        break;
        case 5:
            //optimal root
            if (tree->toString(true) == "(90CF11697,NDK,(Q23,U455))") return -4000.0;
            //heuristic root
            if (tree->toString(true) == "(90CF11697,U455,(NDK,Q23))") return -4500.0;
            //over_bound
            if (tree->toString(true) == "(NDK,U455,(90CF11697,Q23))") return -8000.0;
            cerr << "TREE NOT RECOGNIZED DURING DEBUG" << endl;
            return -8000.0;
        break;
        case 3: //U455 NDK and 90CF11697
            return -3000.0;
        break;
        default: assert(0);
    }
}
#endif


using namespace std;

//register the tree to the tree factory with the name Branch-and-bound
SearchTreeBranchBound SearchTreeBranchBound::prototype("Branch-and-bound");


SearchTreeBranchBound::SearchTreeBranchBound( const string & registrationName ) :
SearchTreeInsertionBasic(registrationName){
}

SearchTreeBranchBound::SearchTreeBranchBound(ParametersSet& parameters) :
InferenceTree(),OptimizerTree(),SearchTreeInsertionBasic(parameters){
    //instantiate internally a Heuristic stepwise addition ML tree in order
    //to find an initial good tree whose likelihood can be used as an acceptable
    //bound
    Singleton < Factory<SearchTree> > & treeFactory = Singleton < Factory<SearchTree> >::instance();
    //transmit parameters (ie, clusters at the moment) to the heuristic tree
    heuristictree = treeFactory.create( "Stepwise addition", parameters );
}

SearchTreeBranchBound::~SearchTreeBranchBound(){
}

bool SearchTreeBranchBound::initialisation( SequenceTable * ptable, Model * pmodel ){
    SearchTreeInsertionBasic::initialisation( ptable, pmodel );

    //find the good upper bound first
    cout << "Use the heuristic stepwise addition to find a good bound" << endl;
    heuristictree->optimizeFlag( optimizeModel );
    if (heuristictree->initialisation( ptable, pmodel )){
        cout << heuristictree->toString(false) << endl;
        cout << "--------------------------------------------------------------------------------" << endl;
    }
#ifndef BRANCHBOUNDDEBUG2
    while(heuristictree->processNext()){
        cout << "-------------------------" << endl;
        cout << heuristictree->toString(false) << endl;
        cout << "--------------------------------------------------------------------------------" << endl;
    }
    bestTree = heuristictree->toString();
#else
    bestTree = "Not visited in debug mode";
#endif

    if ( optimizeModel ){
        pmodel->getAllParameters( bestModelParameters );
    }
    if ( pmodel->getNumberPenaltyParameters() ){
        pmodel->getAllPenaltyParameters( bestModelPenalty );
    }
    bestLnPenalty = pmodel->getLnPenalty();
    bestLnLik = heuristictree->loglikelihood();
    upperBound = bestLnLik + bestLnPenalty;
#ifdef BRANCHBOUNDDEBUG2
    cout << "Hijacking the upper-bound for debugging" << endl;
    upperBound = -7000.0;
#endif
#ifdef BRANCHBOUNDDEBUG1
    //save all the tree visited
    testVisited.open("visited.trees");
#endif
    cout << "Done... initial bound=" << upperBound;
    if (bestLnPenalty){
        cout << "   (likelihood=" << bestLnLik << ", penalty=" << bestLnPenalty << ')';
    }
    cout << endl;
    cout << "*******************************************************" << endl;

    //random order for the addition of species to the tree
    additionOrder.resize(ptable->getNumberSpecies());
    for ( unsigned int i = 0 ; i < additionOrder.size() ; ++i ) {
        additionOrder[i] = i ;
    }
#ifndef BRANCHBOUNDDEBUG2
    additionOrder = statlib::random_permutation( additionOrder );
#else
    additionOrder[0]=3;
    additionOrder[1]=4;
    additionOrder[2]=0;
    additionOrder[3]=2;
    additionOrder[4]=5;
    additionOrder[5]=1;
#endif

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

    //don't display this initial state
    return false;
}



OptimizerTree* SearchTreeBranchBound::clone(ParametersSet& parameters) const{
    return new SearchTreeBranchBound(parameters);
}

bool SearchTreeBranchBound::backTrack(){

    //we need to save the old insertion point in case the backtrack is a "serious" one (ie, not just trying the next branch)
    OptimizerNode* oldInsert = *(branchIter.top());
    //remove the last added node
    removal( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                           oldInsert );
    ++(branchIter.top());

#ifdef BRANCHBOUNDDEBUG1
    assert(branches.size());
    setBranchVector( branches.top() );
    assert(treeSave.size());
    if(toString()!=treeSave.top()){
        cerr << "ACTUAL TREE=" << toString() << endl;
        cerr << "EXPECTED=" << treeSave.top() << endl;
        assert(0);
    }
#endif

    //and continue to remove until we reach a subtree that leads to a new unexplored portion of the topology space
    while(branchIter.top()==branchesStack.top().end()){
        //remove node, reassign clusters
        branchIter.pop();
        branchesStack.pop();
        branches.pop();
        if (optimizeModel){
            modelParameters.pop();
            if (optimizeModelPenalty){
                modelPenalty.pop();
            }
        }
        //finished?
        if (branchIter.empty()){
            return false;
        }
        //the node is removed definitely
        removalValidate( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                                *(branchIter.top()) );
        removal( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                                *(branchIter.top()) );
#ifdef BRANCHBOUNDDEBUG1
        treeSave.pop();
        assert(branches.size());
        setBranchVector( branches.top() );
        assert(treeSave.size());
        if(toString()!=treeSave.top()){
            cerr << "ACTUAL TREE=" << toString() << endl;
            cerr << "EXPECTED=" << treeSave.top() << endl;
            assert(0);
        }
#endif
        ++(branchIter.top());
    }
    //try next, restore the state as it was
    if (optimizeModel){
        pmodel->setAllParameters( modelParameters.top() );
    }
    if (optimizeModelPenalty){
        pmodel->setAllPenaltyParameters( modelPenalty.top() );
    }
    if ( (optimizeModel) || (optimizeModelPenalty) ){
        pmodel->validChange();
    }
    setBranchVector( branches.top() );
    //insert the next node (id=insertionPoints.size()) to the next possible branch
    insertion( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                        *(branchIter.top()) );
    return true;
}

bool SearchTreeBranchBound::buildRec(){

    //optimise the tree as it is
    cout << "optimizing: " << toString(true) << ';' << endl;
#ifndef BRANCHBOUNDDEBUG2
    Optimise::optimiseQuasiNewton( *this, *pmodel, optimizeModel, empiricalFreqs, optimizeModelPenalty ? 1e-9 : 1e-8, false, false );
#else
    vector<double> bv;
    getBranchVector(bv);
    for (unsigned int i = 0; i<bv.size(); ++i){
        bv[i]=(double)rand()/(double)RAND_MAX;
    }
    setBranchVector(bv);
#endif
    cout << toString() << ';' << endl;

#ifdef BRANCHBOUNDDEBUG1
    testVisited << toString(true) << endl;
    double lik = returnValue(this);
#else
    double lik = loglikelihood();
#endif
    double pen = pmodel->getLnPenalty();
    if(optimizeModelPenalty){
        cout << "optimized value=" << lik + pen << "     (";
    }
    cout << "maxLikelihood = " << lik;
    if(optimizeModelPenalty){
        cout << ", maxPenalty=" << pen << ')';
    }
    else{
        //even if we do not optimize penalty parameters, there can be a penalty to output
        if (pen){
            cout << "      (penalty=" << pen << ')';
        }
    }
    cout << endl;
    if (lik+pen<upperBound){
        //upperBound crossed, time for a backTrack
        if (branchesStack.size()+3 == ptable->getNumberSpecies()){
            cout << "tree refused..." << endl;
        }
        else{
            cout << "bound reached, back-tracking..." << endl;
        }
        cout << "--------" << endl;
        bool cont = backTrack();
        //if this is not finish
        if (cont){
            return buildRec();
        }
        else{
            return false;
        }
    }
    else{
        //candidate tree found ?
        if (branchesStack.size()+3 == ptable->getNumberSpecies()){
            cout << "better tree found" << endl;
            upperBound=lik+pen;
            cout << "better tree found, new bound=" << upperBound << endl;
            bestLnLik = lik;
            bestLnPenalty = pen;
            if (optimizeModel){
                pmodel->getAllParameters(  bestModelParameters );
            }
            if (pmodel->getNumberPenaltyParameters()){
                pmodel->getAllPenaltyParameters( bestModelPenalty );
            }
            //validate the insertion of the last node
            insertionValidate( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                               *(branchIter.top()) );
            bestTree = toString();
            //return this partial result (the search will restart by a call to processNext)
            return true;
        }
        //usual case, a partial tree was found and we did not cross the upperBound yet
        else{
            cout << "continue search..." << endl;
            cout << "--------" << endl;
            //validate the insertion of the previous node (if it was not the beginning, ie the third node)
            if (branchesStack.size()){
                insertionValidate( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                                    *(branchIter.top()) );
            }
           //push the current parameter in the stacks (useful in case of backtrack)
            branches.push( vector<double>() );
            getBranchVector( branches.top() );
            if (optimizeModel){
                modelParameters.push( vector<double>() );
                pmodel->getAllParameters( modelParameters.top() );
            }
            if (pmodel->getNumberPenaltyParameters()){
                modelPenalty.push( vector<double>() );
                pmodel->getAllPenaltyParameters( modelPenalty.top() );
            }
#ifdef BRANCHBOUNDDEBUG1
            cout << "PUSH" << endl;
            treeSave.push(toString());
#endif

            //and now let's continue
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
            insertion( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                                *(branchIter.top()) );
            return buildRec();
        }
    }
}

bool SearchTreeBranchBound::processNext(){
    //the algorithm:
    //The tree is constructed by successive addition of species and backtracking
    //if the likelihood of a subtopology is found to be lower than the stored
    //bound then we backtrack and we avoid the useless optimisation of many
    //sub-optimal topologies (because the likelihood MUST decrease when a species is
    //added to an existing tree so it is not possible to find something better)
    //if at some point all leaves have been added and the likelihood is found
    //to be better than the stored bound, then the bound is updated and
    //processNext returns true
    //processNext returns false once he finished visiting the topology space
    //(minus the tree avoided with the help of the upperBound)

/*
    //A QUICK TEST

    vector<double> bi;
    bi.resize(3);

    bi[0]=3;
    bi[1]=5;
    bi[2]=7;
    setBranchVector(bi);

    stack<string> treeF;

    cout << toString() << endl;
    vector< InferenceNode* > allowedBranches1;
    retrievePossibleBranches( allowedBranches1, additionOrder[3] );
    treeF.push( toString() );
    for (unsigned int i = 0; i < allowedBranches1.size(); ++i ){
        insertion( insertNodes[0].second, insertNodes[0].first, (OptimizerNode*)allowedBranches1[i] );
        cout << "BEF1=" << toString() << endl;
        insertionValidate( insertNodes[0].second, insertNodes[0].first, (OptimizerNode*)allowedBranches1[i] );
        cout << "AFT1=" << toString() << endl;
        treeF.push( toString() );
        vector< InferenceNode* > allowedBranches2;
        retrievePossibleBranches( allowedBranches2, additionOrder[4] );
        for (unsigned int j = 0; j < allowedBranches2.size(); ++j ){
            insertion( insertNodes[1].second, insertNodes[1].first, (OptimizerNode*)allowedBranches2[j] );
            cout << "BEF2=" << toString() << endl;
            insertionValidate( insertNodes[1].second, insertNodes[1].first, (OptimizerNode*)allowedBranches2[j] );
            cout << "AFT2=" << toString() << endl;
            treeF.push( toString() );
            vector< InferenceNode* > allowedBranches3;
            retrievePossibleBranches( allowedBranches3, additionOrder[5] );
            for (unsigned int k = 0; k < allowedBranches3.size(); ++k ){
                insertion( insertNodes[2].second, insertNodes[2].first, (OptimizerNode*)allowedBranches3[k] );
                cout << "FINAL=" << toString() << endl;
                removal( insertNodes[2].second, insertNodes[2].first, (OptimizerNode*)allowedBranches3[k] );
                assert( treeF.top() == toString() );
            }
            treeF.pop();
            cout << "RETOURAFT2=" << toString() << endl;
            removalValidate( insertNodes[1].second, insertNodes[1].first, (OptimizerNode*)allowedBranches2[j] );
            cout << "RETOURBEF2=" << toString() << endl;
            removal( insertNodes[1].second, insertNodes[1].first, (OptimizerNode*)allowedBranches2[j] );
            assert( treeF.top() == toString() );
        }
        treeF.pop();
        cout << "RETOURAFT1=" << toString() << endl;
        removalValidate( insertNodes[0].second, insertNodes[0].first, (OptimizerNode*)allowedBranches1[i] );
        cout << "RETOURBEF1=" << toString() << endl;
        removal( insertNodes[0].second, insertNodes[0].first, (OptimizerNode*)allowedBranches1[i] );
        assert( treeF.top() == toString() );
    }
    treeF.pop();
    exit(EXIT_FAILURE);
*/

    bool contin = true;

    //each time processNext is called a backtrack is needed.
    //(each time but the first one)
    if ( !branchIter.empty() ){
        //the last accepted leaf in the last accepted tree was validated (ie, the tree was reconstructed)...unvalidate
        removalValidate( insertNodes[branchesStack.size()-1].second, insertNodes[branchesStack.size()-1].first,
                               *(branchIter.top()) );
        //backTrack destroys the current tree (it should be a complete one) and
        //return to the next state we have to try.
        contin = backTrack();
    }
    //buildRec continues the building procedure at the point designed by the current stacks if it is not finish
    if (contin==true){
        contin=buildRec();
    }
    return contin;
}

void SearchTreeBranchBound::printResults( ostream& outputStream ){
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
        //we did not optimize penalty parameters, nevertheless there is can be a best penalty to output
        if (bestLnPenalty){
            outputStream << "      (penalty=" << bestLnPenalty << ')';
        }
    }
    outputStream << endl;
}

void SearchTreeBranchBound::printTree( ostream& outputStream ){
    // Save a tree to a file on its own
    outputStream << bestTree << ';' << endl;
}
