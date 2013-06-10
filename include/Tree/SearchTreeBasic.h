#ifndef SEARCHTREEBASIC_H
#define SEARCHTREEBASIC_H

#include "Tree/SearchTree.h"
#include "Tree/OptimizerTreeBasic.h"
#include "Tree/ClustersTree.h"

#include <map>
#include <deque>

class OptimizerNode;

/** ****************************************************************************
 * Basic methods to perform a ML search.
 * This class extends OptimizerTreeBasic. Its children define methods to
 * perform the search into the tree space but mlphase has to rely on
 * the class optimise to optimise the likelihood for each visited tree.
 * (since the optimisation depends on the substitution model too)
 **************************************************************************** */

/*because of its inheritance, we must explicitely define the method
 *from OptimizerTree in this class */
class SearchTreeBasic : virtual public InferenceTree, virtual public OptimizerTree, public SearchTree, public OptimizerTreeBasic{
protected:
    /** ************************************************************************
     * SearchTreeBasic
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< MLSearchTree >
     ************************************************************************ */
    SearchTreeBasic( const string & registrationName );

    /** ************************************************************************
     * SearchTreeBasic
     * @semantics  constructor of an empty MLSearchTree according to the
     *             parameters, this constructor should be called by all
     *             descendant
     ************************************************************************ */
    SearchTreeBasic( ParametersSet& treeParameters );


    bool optimizeModel;
    bool optimizeModelPenalty;
    bool empiricalFreqs;


public:
    /** ************************************************************************
     * ~SearchTreeBasic
     * @semantics  destructor of a ML (search) inference tree
     ************************************************************************ */
    virtual ~SearchTreeBasic();

    /** ************************************************************************
     * initialisation
     * @input          the sequences and the model to use with this tree
     * @preconditions  no data or model already loaded
     * @semantics      call OptimizerTreeBasic::loadDataAndModel
     * @return         true if there is an initial tree
     ************************************************************************ */
    virtual bool initialisation( SequenceTable * ptable, Model * pmodel );

    /** ************************************************************************
     * diffloglikelihood
     * @semantics  compute the gradient of the likelihood function
     *             at the actual point
     ************************************************************************ */
    virtual void diffloglikelihood(vector < double > & gradientVector);

    /** ************************************************************************
     * optimizeFlag
     * @semantics  load the model and
     *             set the model optimization flag
     ************************************************************************ */
    inline virtual void optimizeFlag( bool flag ){
        optimizeModel = flag;
    }

    /** ************************************************************************
     * empiricalFreqsFlag
     * @semantics  load the model and
     *             set the empirical frequencies flag
     ************************************************************************ */
    inline virtual void empiricalFreqsFlag( bool flag ){
        empiricalFreqs = flag;
    }

    /** ************************************************************************
     * processNext
     * @semantics  go to the next step of the search
     * @return     bool false if the search is over.
     ************************************************************************ */
    virtual bool processNext() = 0;

    /** ************************************************************************
     * printResults
     * @semantics  print the results
     ************************************************************************ */
    virtual void printResults( ostream& outputStream ) = 0;

};

#endif //SEARCHTREEBASIC_H
