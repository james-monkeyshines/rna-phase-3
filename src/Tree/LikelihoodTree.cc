#include "Tree/LikelihoodTree.h"

#include "Models/Model.h"

#include <assert.h>
#include <math.h>

#include <iostream>

using namespace std;

LikelihoodTree::LikelihoodTree( string stringTree ): InferenceTree(stringTree){
    for ( vector< BasicNode*>::const_iterator iter = nodeRefVector.begin();
          iter != nodeRefVector.end(); ++iter ){
        if ( (*iter!=root) && ((*iter)->getParentDistance() < 0.0) ){
            cerr << "Error: forgotten or negative branch lengths during the initialization of the tree" << endl;
            exit(EXIT_FAILURE);
        }
    }
}

void LikelihoodTree::getAncestralStatesRec( InferenceNode* node, InferenceNode* oldFather, list< pair<InferenceNode*, vector< array2D<double> > > > & ancestralStates ){

    //node is supposed to be a temporary "root" when the algorithm is called but root will
    //not change accordingly and stay to the same point
    assert(node->getParent()==NULL);
    //recompute partialLikelihood at this node, only the root will escape this since there is no point to invalidate it (yet)
    if (node!=root){
        node->partialLikelihood=node->partialLikelihoodWork;
        //increase lookup table by 1 since there is a new son
        node->lookup.push_back(NULL);
        node->likelihood();
    }
    
    ancestralStates.push_back( pair<InferenceNode*, vector< array2D<double> > > (node,vector< array2D<double> >()) );
    vector< array2D<double> > & ancestralSeqProbs = ancestralStates.back().second;
    ancestralSeqProbs.resize( ptable->getNumberCategories() );

    for( unsigned int cat = 0; cat < ptable->getNumberCategories(); ++cat ){
        unsigned int numberRateCategories = pmodel->getNumberRatesCategories(cat);
        double rateProb[numberRateCategories];
        for ( unsigned int rateCat = 0; rateCat < numberRateCategories; ++rateCat ) {
            rateProb[rateCat] = pmodel->getRateCategoryProbability( rateCat, cat );
        }
        unsigned int numberStates = pmodel->getNumberStates(cat);
        double frequency[numberStates][numberRateCategories];
        for ( unsigned int state = 0; state < numberStates; ++state ) {
            for ( unsigned int rateCat = 0; rateCat < numberRateCategories; ++rateCat ) {
                frequency[state] [rateCat] = pmodel->getFrequency( state,
                                                         rateCat, cat );
            }
        }
        unsigned int sequenceLength = ptable->getSequencesLength(cat);
        unsigned int invariantLength = ptable->getInvariantsLength(cat);
        ancestralSeqProbs[cat].resize( sequenceLength + invariantLength, numberStates );
        
        for ( unsigned int site = 0; site < sequenceLength; ++site ){
            double tot = 0.0;
            for ( unsigned int nucleotide = 0; nucleotide < numberStates; ++nucleotide ) {
                ancestralSeqProbs[cat](site,nucleotide) = 0.0;
                for ( unsigned int rate = 0; rate < numberRateCategories; ++rate ) {
                    ancestralSeqProbs[cat](site,nucleotide) += rateProb[rate] *
                        (*(node->partialLikelihood[cat]))(site, nucleotide, rate) *
                         frequency[nucleotide][rate];
                }
                tot += ancestralSeqProbs[cat](site,nucleotide);
            } // end rate
            for ( unsigned int nucleotide = 0; nucleotide < numberStates; ++nucleotide ){
                ancestralSeqProbs[cat](site,nucleotide) /= tot;
            }
        }//end site

        //invariant
        const vector< pair< string, unsigned int > > &inv = ptable->getInvariantBases(cat);
        for ( unsigned int invId = 0; invId < inv.size(); ++invId ) {
            double tot = 0.0;
            for ( unsigned int nucleotide = 0; nucleotide < numberStates;
                  ++nucleotide ) {
                ancestralSeqProbs[cat](sequenceLength+invId,nucleotide) = 0.0;
                for ( unsigned int rate = 0; rate < numberRateCategories; ++rate ) {
                    ancestralSeqProbs[cat](sequenceLength+invId,nucleotide) +=
                            rateProb[rate] *
                        (*(node->partialLikelihood[cat]))(sequenceLength+invId, nucleotide, rate) *
                         frequency[nucleotide] [rate];
                }
                tot += ancestralSeqProbs[cat](sequenceLength+invId,nucleotide);
            }
            for ( unsigned int nucleotide = 0; nucleotide < numberStates; ++nucleotide ) {
                ancestralSeqProbs[cat](sequenceLength+invId,nucleotide) /= tot;
            }
        } // end invId
    } // end for all categories
    
    //now the tricky stuff, transmit call back down to other internal nodes
    //a nice use of partialLikelihoodWork and partialLikelihoodSave to do
    //the marginal reconstruction quickly in only 2 passes over the tree
    //using what is already implemented
    
    //if not already done invalidate yourself (save will keep the good computation)
    //only the root should invalid itself at that point
    if (node==root){
        assert(node->partialLikelihood==node->partialLikelihoodSave);
        node->partialLikelihood = node->partialLikelihoodWork;
    }
    else{
        assert(node->partialLikelihood==node->partialLikelihoodWork);
    }

    list<BasicNode*> listChildren = node->getChildrenList();
    for( list<BasicNode*>::iterator iter = listChildren.begin();
         iter != listChildren.end(); ++iter ){
        //stop the recursive call if a leaf is reached
        //and do not send back the call to the previous father which is now a child
        if ( (!(*iter)->isLeaf()) && (*iter!=oldFather) ){
            double dist = (*iter)->getParentDistance();
            node->removeChild((*iter));
            //recompute a partial likelihood for node, it's old father (if any) is supposed
            //to have done that for itself before the recursive call
            node->likelihood();
            (*iter)->addChild( node, dist );
            getAncestralStatesRec( (InferenceNode*)(*iter), node, ancestralStates );
            (*iter)->removeChild( node );
            node->addChild( *iter, dist );
        }
    }
    //recover the old partial computations, so that the father can recompute a proper likelihood when it
    //will put himself under its other children
    node->partialLikelihood = node->partialLikelihoodSave;
}


void LikelihoodTree::loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel ){
    InferenceTree::loadDataAndModel( psequenceTable, pmodel );
    //in case the given tree is not fully resolved...
    lookupResize();
}

void LikelihoodTree::getAncestralStates( list< pair<InferenceNode*, vector< array2D<double> > > > & ancestralStates ){
    assert(ptable);
    //We can compute directly the ancestral sequence at the root and we have to "send back
    //the partialLikelihood array back to the tips" (so that each node becomes root at its turn)

    //save the current computations in partialLikelihood
    //(for all nodes, partialLikelihoodWork will switch to the other workspace
    //and partialLikelihood will be equal to partialLikelihoodSave
    saveNodes(-1);
    assert(((InferenceNode*)root)->partialLikelihood==((InferenceNode*)root)->partialLikelihoodSave);
    
    //becareful with the size of lookup...
    lookup.resize(lookup.size()+1);
    lookup.back() = new array3D<double> [ptable->getNumberCategories()];
    for( unsigned int cat = 0; cat < ptable->getNumberCategories(); ++cat ){
        lookup.back()[cat].resize(
                pmodel->getNumberRatesCategories(cat),
                pmodel->getNumberStates(cat),
                pmodel->getNumberStates(cat) );
    }
    getAncestralStatesRec( (InferenceNode*)root, NULL, ancestralStates );

}


void LikelihoodTree::getRateCategory( vector< array2D<double> >& rateCat ){
    rateCat.resize( ptable->getNumberCategories() );

    for( unsigned int cat = 0; cat < ptable->getNumberCategories(); ++cat ){
        int sequenceLength = ptable->getSequencesLength(cat);
        int numberRateCategories = pmodel->getNumberRatesCategories(cat);
        double rateProb[numberRateCategories];
        for ( int rate = 0; rate < numberRateCategories; ++rate ) {
            rateProb[rate] =
            pmodel->getRateCategoryProbability( rate, cat );
        }
        int numberStates = pmodel->getNumberStates(cat);
        double frequency[numberStates][numberRateCategories];
        for ( int state = 0; state < numberStates; ++state ) {
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                frequency[state] [rate] = pmodel->getFrequency( state,
                                                         rate, cat );
            }
        }

        //invariant
        const vector< pair< string, unsigned int > > &inv =
            ptable->getInvariantBases(cat);

        rateCat[cat].resize(sequenceLength + inv.size(), numberRateCategories);

        for ( int site = 0; site < sequenceLength; ++site ) {
            double tot = 0.0;
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                double tempn = 0.0;
                for ( int nucleotide = 0; nucleotide < numberStates;
                      ++nucleotide ) {
                    tempn += frequency[nucleotide][rate]*
                        (*(((InferenceNode*)root)->
                           partialLikelihood[cat]))(site, nucleotide, rate);
                }
                rateCat[cat](site, rate) = tempn * rateProb[rate];
                tot += rateCat[cat](site, rate);
            } // end rate
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                rateCat[cat](site, rate) /= tot;
            }
        }//end site


        for ( unsigned int invId = 0; invId < inv.size(); ++invId ) {
            double tot = 0.0;
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                double tempn = 0.0;
                for ( int nucleotide = 0; nucleotide < numberStates;
                      ++nucleotide ) {
                    tempn += frequency[nucleotide] [rate]  *
                        (*(((InferenceNode*)root)->
                           partialLikelihood[cat]))(sequenceLength+invId, nucleotide, rate);
                }
                rateCat[cat](sequenceLength+invId, rate) = tempn * rateProb[rate];
                tot += rateCat[cat](sequenceLength+invId, rate);
            }
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                rateCat[cat](sequenceLength+invId, rate) /= tot;
            }
        } // end invId
    } // end for all categories
}

void LikelihoodTree::getSiteSpecificLikelihood( vector< vector<double> >& siteSpecificLik ){
    siteSpecificLik.resize( ptable->getNumberCategories() );
    for( unsigned int cat = 0; cat < ptable->getNumberCategories(); ++cat ){
        siteSpecificLik[cat].clear();
        int sequenceLength = ptable->getSequencesLength(cat);
        int numberRateCategories = pmodel->getNumberRatesCategories(cat);
        double rateProb[numberRateCategories];
        for ( int rateCat = 0; rateCat < numberRateCategories; ++rateCat ) {
            rateProb[rateCat] =
            pmodel->getRateCategoryProbability( rateCat, cat );
        }
        int numberStates = pmodel->getNumberStates(cat);
        double frequency[numberStates][numberRateCategories];
        for ( int state = 0; state < numberStates; ++state ) {
            for ( int rateCat = 0; rateCat < numberRateCategories; ++rateCat ) {
                frequency[state] [rateCat] = pmodel->getFrequency( state,
                                                         rateCat, cat );
            }
        }

        for ( int site = 0; site < sequenceLength; ++site ) {
            double tempr = 0.0;
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                double tempn = 0.0;
                for ( int nucleotide = 0; nucleotide < numberStates;
                      ++nucleotide ) {
                    tempn += frequency[nucleotide][rate]*
                        (*(((InferenceNode*)root)->
                           partialLikelihood[cat]))(site, nucleotide, rate);
                }
                tempr += ( rateProb[rate] * tempn );
            } // end rate
            siteSpecificLik[cat].push_back( log( tempr ) );
        }//end site

        //invariant
        const vector< pair< string, unsigned int > > &inv =
            ptable->getInvariantBases(cat);
        for ( unsigned int invId = 0; invId < inv.size(); ++invId ) {
            double tempr = 0.0;
            for ( int rate = 0; rate < numberRateCategories; ++rate ) {
                double tempn = 0.0;
                for ( int nucleotide = 0; nucleotide < numberStates;
                      ++nucleotide ) {
                    tempn += frequency[nucleotide] [rate]  *
                    (*(((InferenceNode*)root)->
                           partialLikelihood[cat]))
                                  (sequenceLength+invId, nucleotide, rate);
                }
                tempr += ( rateProb[rate] * tempn );
            }
            for (unsigned int nb=0; nb<inv[invId].second; ++nb){
                siteSpecificLik[cat].push_back( log( tempr ) );
            }
        } // end invId
    }
}
