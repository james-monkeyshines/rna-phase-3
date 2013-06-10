#include "Tree/OptimizerNode.h"

#include "Tree/OptimizerTree.h"
#include "Tree/TreeMap.h"
#include "Models/Model.h"

#include <assert.h>
#include <math.h>

OptimizerNode::OptimizerNode( OptimizerTree* ptree ):InferenceNode(ptree){
    if ( ptree->getTable() ){
        initialisationOpti( false );
    }
}

OptimizerNode::OptimizerNode( const string& label, OptimizerTree* ptree ):
        InferenceNode(label,ptree){
    if ( ptree->getTable() ){
        initialisationOpti( true );
    }
}

OptimizerNode::OptimizerNode( const TreeMap& treeMap, OptimizerTree* ptree ):
InferenceNode(){
    //if this node is an internal node
    if ( treeMap.getNumberChildren() ){
        cnumber = -1;
        this->ptree = ptree;
        //create the children
        for ( unsigned int i = 0; i < treeMap.getNumberChildren() ; ++i ){
            addChild( new OptimizerNode(treeMap.getChildMap(i), ptree),
                      treeMap.getDistance(i) );
        }
        updateNodesCount();
        //launch the initialisation if the table/model was loaded
        if ( ptree->getTable() ){
            assert( ptree->getModel() );
            initialisation( false );
            initialisationOpti( false );
        }
    }
    //if this node is a leaf
    else{
        //the basic node constructor used was for an internal node
        //correct this
        nbInternalNodes = 0;
        nbLeaves = 1;
        cnumber = -1;
        this->label = treeMap.getLabel();
        this->ptree = ptree;
        //launch the initialisation if the table/model was loaded
        if ( ptree->getTable() ){
            assert( ptree->getModel() );
            initialisation( true );
            initialisationOpti( true );
        }
    }
}

OptimizerNode::~OptimizerNode(){
}


void OptimizerNode::initialisationOpti( bool leafNode ) {
    unsigned int numberCategories = ptree->getTable()->getNumberCategories();
    if(!leafNode){
        belief.resize( numberCategories );
    }
    prior.resize( numberCategories );
    for (unsigned int cat = 0; cat < numberCategories; ++cat){
        unsigned int totalLength = ptree->getTable()->getSequencesLength(cat) +
                                ptree->getTable()->getInvariantsLength(cat);
        int numberStates = pnodeModel->getNumberStates(cat);
        int numberRatesCategories = pnodeModel->getNumberRatesCategories(cat);
        if(!leafNode){
            belief[cat].resize( numberRatesCategories,
                                numberStates, totalLength );
        }
        prior[cat].resize( numberRatesCategories,
                            numberStates, totalLength );
    }
}


void OptimizerNode::computeBeliefPrior(){
    unsigned int numberCategories = ptree->getTable()->getNumberCategories();

    SequenceTable* ptable = ptree->getTable();

    //if the node is NOT the root fill lookup table
    if(getParent()){
        for( unsigned int cat = 0; cat < numberCategories; ++cat ){
            int numberStates = pnodeModel->getNumberStates(cat);
            int numberRatesCategories = pnodeModel->getNumberRatesCategories(cat);

            // Initialize lookup table
            for ( int initialState = 0; initialState < numberStates; ++initialState ){
                for ( int finalState = 0; finalState < numberStates; ++finalState ){
                    for ( int rate = 0; rate < numberRatesCategories; ++rate ) {
                        ptree->lookup[0][cat](rate, initialState, finalState ) =
                            pnodeModel->probability( initialState, finalState,
                                getParentDistance(), rate, cat );
                    }
                }
            }
        }
    }


    // Calculate the belief if internal node and propagate the call down the
    // tree. During this recursive call parents initialise their belief array
    // so that children can use it and children initialise the prior array
    // so that parent can use it when the call returns.
    if(!isLeaf()){
        for( unsigned int cat = 0; cat < numberCategories; ++cat ){
            lookup[0] = &(ptree->lookup[0][cat]);
            unsigned int sequencesLength = ptable->getSequencesLength(cat);
            unsigned int invariantsLength = ptable->getInvariantsLength(cat);
            unsigned int totalLength = sequencesLength + invariantsLength;
            int numberStates = pnodeModel->getNumberStates(cat);
            int numberRatesCategories = pnodeModel->getNumberRatesCategories(cat);

            //if the node is the root
            if(!getParent()){
                for ( int rate = 0; rate < numberRatesCategories; ++rate ) {
                    for ( unsigned int site = 0; site < totalLength; ++site ) {
                        for ( int state = 0; state < numberStates; ++state ) {
                           belief[cat]( rate, state, site ) =
                                ptree->getModel()->getFrequency( state, rate, cat) *
                                (*partialLikelihood[cat])(site, state, rate );
#ifdef DEBUG1
                            assert( !isnan(
                              belief[cat]( rate, state, site ) ) ) ;
                            assert( !isinf(
                              belief[cat]( rate, state, site ) ) ) ;
#endif
                        }
                    }
                } //end rate
            } // end root
            else{
                //initialise temp
                double temp[numberStates];

                for ( int rate = 0; rate < numberRatesCategories; ++rate ) {
                    for ( unsigned int site = 0; site < totalLength; ++site ){
                        for ( int initialState = 0; initialState < numberStates;
                              ++initialState ){
                            double contrib = 0.0;
                            for ( int finalState = 0; finalState < numberStates;
                                  ++finalState ){
                                contrib +=
                                    (*lookup[0])(rate, initialState, finalState ) *
                                    (*partialLikelihood[cat])( site, finalState, rate);
                            }
                            //use belief of the parent
                            if (((OptimizerNode*)getParent())->belief[cat](rate, initialState, site)==0.0){
                                temp[initialState] = 0.0;
                            }
                            else{
                                temp[initialState] =
                                ((OptimizerNode*)getParent())->belief[cat](rate,initialState, site)/
                                                   contrib;
                            }
#ifdef DEBUG1
                            assert( !isnan( temp[initialState] ) ) ;
                            assert( !isinf( temp[initialState] ) ) ;
#endif
                        } //end initialState

                        for (int finalState=0 ; finalState<numberStates; ++finalState){
                            double contrib = 0.0;
                            for (int initialState=0 ; initialState<numberStates;
                                 ++initialState){
                                contrib += (*lookup[0])(rate, initialState, finalState ) *
                                           temp[initialState] ;
                            }
                            //fill the belief array of this node
                            belief[cat](rate, finalState, site) = contrib *
                               (*partialLikelihood[cat])(site, finalState, rate);
#ifdef DEBUG1
                            assert(!isnan(belief[cat](rate, finalState, site)));
                            assert(!isinf(belief[cat](rate, finalState, site)));
#endif
                        } // finalState
                    } //sites
                } //rates
            } //belief internal node
        } //for each category


        //propagate computeBeliefPrior to the children
        for( list<BasicNode*>::iterator childrenIter = getChildrenList().begin();
                     childrenIter != getChildrenList().end(); ++childrenIter ){
            ((OptimizerNode*)(*childrenIter))->computeBeliefPrior();
        }
    } // if !isLeaf()

    //time to compute the priors now (not for the root)
    if (getParent()){

        for( unsigned int cat = 0; cat < numberCategories; ++cat ){
            int numberStates = pnodeModel->getNumberStates(cat);
            int numberRatesCategories = pnodeModel->getNumberRatesCategories(cat);
            unsigned int sequencesLength = ptable->getSequencesLength(cat);
            unsigned int invariantsLength = ptable->getInvariantsLength(cat);
            unsigned int totalLength = sequencesLength + invariantsLength;

            // WE HAVE TO RECOMPUTE THE LOOKUP table because it was erased
            // during the recursive call
            for ( int initialState = 0; initialState < numberStates; ++initialState ){
                for ( int finalState = 0; finalState < numberStates; ++finalState ){
                    for ( int rate = 0; rate < numberRatesCategories; ++rate ) {
                        ptree->lookup[0][cat](rate, initialState, finalState ) =
                            pnodeModel->probability( initialState, finalState,
                            getParentDistance(), rate, cat );
                    }
                }
            }

            lookup[0] = &(ptree->lookup[0][cat]);

            for ( int rate = 0; rate < numberRatesCategories; ++rate ){
                // if a children is a leaf, rates are not in the partial
                // Likelihood array and endRate = 0. Otherwise dim
                // partialLikelihood[site] = (numberStates, numberRates)
                // and "endRate = startRate"
                int endRate = isLeaf() ? 0 : rate;
                for ( unsigned int site = 0; site < totalLength; ++site ){
                    for ( int initialState = 0; initialState < numberStates;
                          ++initialState ) {
                        double contrib = 0.0;
                        for ( int finalState = 0; finalState < numberStates;
                              ++finalState ){
                            contrib += (*lookup[0])(rate, initialState, finalState ) *
                                   (*partialLikelihood[cat])(site, finalState, endRate);
                         } //finalState
                         if( ((OptimizerNode*)getParent())->belief[cat](rate, initialState, site) == 0.0 ){
                             prior[cat](rate, initialState, site) = 0.0;
                         }
                         else{
                             prior[cat](rate, initialState, site) =
                                ((OptimizerNode*)getParent())->belief[cat](rate, initialState, site)/contrib;
                             assert(!isnan(prior[cat](rate, initialState, site)));
                             assert(!isinf(prior[cat](rate, initialState, site)));
                         }
                    } //initial state
                } //site
            } //rate
        } // priors for each symbol category
    } //end if !root
}
