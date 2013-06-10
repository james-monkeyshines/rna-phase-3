#include "Tree/OptimizerTreeBasic.h"

#include "PatternDesign/Singleton.h"
#include "PatternDesign/Factory.h"

#include "Util/FileParser.h"
#include "Util/ParametersSet.h"

#include "Tree/TreeMap.h"
#include "Tree/OptimizerNode.h"
#include "Models/Model.h"

OptimizerTreeBasic::OptimizerTreeBasic():UnrootedTree(){
}

OptimizerTreeBasic::OptimizerTreeBasic(const string & registrationName) : UnrootedTree(){
    Singleton < Factory<OptimizerTree> > & treeFactory = Singleton < Factory<OptimizerTree> >::instance();
    treeFactory.subscribe( this, registrationName );
}

OptimizerTreeBasic::OptimizerTreeBasic(ParametersSet& parameters):UnrootedTree(parameters){
}

void OptimizerTreeBasic::constructFromString( string stringTree ){
    assert(ptable);
    // Create a depth map of the tree
    TreeMap map( stringTree );
    // Create the nodes
    root = new OptimizerNode( map, this );
    createIndex();
    if (clustersTree){
        assignClusters();
    }
    if( unroot() ){
        cerr << "WARNING: your tree has been unrooted" << endl;
    }
    
    //checkBinary();
    createOutgroup();
    //in case the given tree is not fully resolved...
    lookupResize();
}

OptimizerTreeBasic::~OptimizerTreeBasic(){
}



void OptimizerTreeBasic::diffloglikelihood(vector < double > & gradientVector) {
    assert ( pmodel && ptable );

    gradientVector.resize( getNumberBranches() );
    for ( unsigned int i = 0; i < getNumberBranches(); ++i ){
        gradientVector[i] = 0.0;
    }

    //BEWARE, the order of processing is important
    //compute the children likelihood before
    for ( vector< BasicNode*>::reverse_iterator iter = nodeRefVector.rbegin();
          iter != nodeRefVector.rend(); ++iter ){
        if (!(*iter)->isLeaf()){
            ((InferenceNode*)(*iter))->likelihood();
        }
    }
    //fill belief and prior for the nodes in the tree (recursive call)
    ((OptimizerNode*)root)->computeBeliefPrior();

    // calculate the gradients
    unsigned int numberCategories = ptable->getNumberCategories();
    vector < BasicNode * >::const_iterator iter = nodeRefVector.begin();
    ++iter;
    unsigned int i = 0;
    while ( iter != nodeRefVector.end() ){
        for (unsigned int cat = 0; cat < numberCategories; ++cat){
            int numberRatesCategories = pmodel->getNumberRatesCategories(cat);
            int numberStates = pmodel->getNumberStates(cat);
            unsigned int sequencesLength = ptable->getSequencesLength(cat);
            unsigned int invariantsLength = ptable->getInvariantsLength(cat);
            unsigned int totalLength = sequencesLength + invariantsLength;
            // Initialise lookup tables (probability and diffProb)
            for (int rate = 0; rate < numberRatesCategories; ++rate){
                for ( int initialState = 0; initialState < numberStates;
                      ++initialState ){
                    for ( int finalState = 0; finalState < numberStates;
                          ++finalState ){
                        lookup[0][cat]( rate, initialState, finalState ) =
                            ((OptimizerNode*)(*iter))->getModel()->probability( initialState, finalState,
                                           (*iter)->getParentDistance(), rate, cat );
                        lookup[1][cat]( rate, initialState, finalState ) =
                            ((OptimizerNode*)(*iter))->getModel()->diffProbability( initialState, finalState,
                                            (*iter)->getParentDistance(), rate, cat );
                    }
                }
            } //end initialisation two lookup tables

            for ( unsigned int site = 0; site < totalLength; ++site){
                double siteLikelihood = 0.0;
                double siteGradient = 0.0;
                for ( int rate = 0; rate < numberRatesCategories; ++rate ){
                    double sumdiff = 0.0;
                    double suml = 0.0;
                    if ( !(*iter)->isLeaf() ){
                        for ( int initialState = 0; initialState < numberStates;
                          ++initialState ){
                            for ( int finalState = 0; finalState < numberStates;
                                  ++finalState ){
                                double priorXpartialLikelihood =
                                     ((OptimizerNode*)(*iter))->prior[cat](rate, initialState, site) *
                                     (*((OptimizerNode*)(*iter))->partialLikelihood[cat])( site, finalState, rate );
                                sumdiff += lookup[1][cat](rate, initialState, finalState) * priorXpartialLikelihood;
                                suml += lookup[0][cat](rate, initialState, finalState) * priorXpartialLikelihood;
                            } // end final state
                        } // end initial state
                    }
                    //if the child is a leaf use the unique partial likelihood
                    else{
                        for ( int initialState = 0; initialState < numberStates;
                          ++initialState ){
                            for ( int finalState = 0; finalState < numberStates;
                                  ++finalState ){
                                double priorXpartialLikelihood =
                                     ((OptimizerNode*)(*iter))->prior[cat](rate, initialState, site) *
                                     (*((OptimizerNode*)(*iter))->partialLikelihood[cat])( site, finalState, 0 );
                                sumdiff += lookup[1][cat](rate, initialState, finalState) * priorXpartialLikelihood;
                                suml += lookup[0][cat](rate, initialState, finalState) * priorXpartialLikelihood;
                            } // end final state
                        } // end initial state
                    }
                    siteGradient += pmodel->getRateCategoryProbability(rate, cat) * sumdiff;
                    siteLikelihood += pmodel->getRateCategoryProbability(rate, cat) * suml;
#ifdef DEBUG1
                    assert( !isnan( siteGradient ) );
                    assert( !isinf( siteGradient ) );
                    assert( !isnan( siteLikelihood ) );
                    assert( !isinf( siteLikelihood ) );
#endif
                } // end rate
                if (site<sequencesLength){
                    gradientVector[i] += ( 1.0 / siteLikelihood * siteGradient );
                }
                else{
                    gradientVector[i] += (1.0 / siteLikelihood * siteGradient) *
                        ptable->getInvariantBases(cat)[site-sequencesLength].second;
                }
            } // site
        } // end symbol category
        ++iter;
        ++i;
    } //end branch


#ifdef DEBUG1
    for (unsigned int i = 0; i < gradientVector.size(); ++i){
        assert( !isnan( gradientVector[i] ) );
        assert( !isinf( gradientVector[i] ) );
    }
#endif
}

void OptimizerTreeBasic::loadDataAndModel( SequenceTable * psequenceTable,
Model * pmodel ) {
    //complete the initialisation for OptimizerNode
    UnrootedTree::loadDataAndModel( psequenceTable, pmodel );
    for ( vector< BasicNode*>::reverse_iterator iter = nodeRefVector.rbegin();
          iter != nodeRefVector.rend(); ++iter ){
        if ( ((*iter)==root) || (*iter)->getNumberChildren() ){
            ((OptimizerNode*)(*iter))->initialisationOpti( false );
        }
        else{
            ((OptimizerNode*)(*iter))->initialisationOpti( true );
        }
    }
}
