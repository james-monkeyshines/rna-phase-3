#ifndef LIKELIHOODTREE_H
#define LIKELIHOODTREE_H

#include "Tree/InferenceTree.h"

#include <vector>

#include "Util/array2D.h"

using namespace std;

/** ****************************************************************************
 * LikelihoodTree
 * A direct descendant of InferenceTree to implement methods related
 * to rate category estimatation and ancestral sequence reconstruction
 * (the aim is just to remove these methods from the InferenceTree class
 * and keep the methods related to post-analysis together)
 * LikelihoodTree might not be a good name but this tree is used by the main
 * program likelihood so...
 **************************************************************************** */
 
class LikelihoodTree : public InferenceTree {
public:
    /** ************************************************************************
     * LikelihoodTree
     * @input       Tree in Newick format (with branch lengths a priori)
     * @semantics   Constructor from a string
     ************************************************************************ */
    LikelihoodTree( string stringTree );

    /** ************************************************************************
     * getAncestralStates
     * @return          vector< Node*, ancestralFreq[modelId](site, probState) > 
     * @preconditions   likelihood was called before to initialise the root
     *                  node partial likelihood array. BEWARE!
     * @semantics       Marginal reconstruction. return for each internal node, 
     *                  and each site the values P(state|X,M) (site in the 
     *                  special sequence provided by sequenceTable)
     ************************************************************************ */
    void getAncestralStates( list< pair<InferenceNode*, vector< array2D<double> > > > & ancestralStates );
    void getAncestralStatesRec( InferenceNode* node, InferenceNode* oldFather, list< pair<InferenceNode*, vector< array2D<double> > > > & ancestralStates );

    /** ************************************************************************
     * getRateCategory
     * @return          rateCat[modelId](site, probRate)
     * @preconditions   likelihood was called before to initialise the root
     *                  node partial likelihood array. BEWARE!
     * @semantics       return for each site the values P(rate|X,M)
     ************************************************************************ */
    void getRateCategory( vector< array2D<double> >& rateCat );

    /** ************************************************************************
     * getSiteSpecificLikelihood
     * @return          siteSpecificLik[modelId][site]
     * @preconditions   likelihood was called before to initialise the root
     *                  node partial likelihood array. BEWARE!
     * @semantics       return for each site the values P(X|M)
     ************************************************************************ */
    void getSiteSpecificLikelihood( vector< vector<double> >& siteSpecificLik );
 
 
    /** ************************************************************************
     * loadDataAndModel
     * @input           the sequences and the model to use with this tree
     * @preconditions   no data or model already loaded
     * @sematics        redefined only to change the size of lookup
     *                  once the model is loaded
     ************************************************************************ */
    virtual void loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel );
    
};
    
#endif //LIKELIHOODTREE_H
