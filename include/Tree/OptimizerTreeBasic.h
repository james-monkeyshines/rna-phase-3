#ifndef OPTIMIZERTREEBASIC_H
#define OPTIMIZERTREEBASIC_H

#include "Tree/UnrootedTree.h"
#include "Tree/OptimizerTree.h"
#include "PatternDesign/Singleton.h"
#include "Util/randombox.h"

/** ****************************************************************************
 * Basic methods to optimise a tree, actually this class is almost the tree used
 * in optimise. It contains method to initialise the tree and a method to compute
 * the derivative of the loglikelihood wrt branch lengths. Optimise uses
 * standard methods from the tree class to optimise the tree otherwise.
 **************************************************************************** */

class OptimizerTreeBasic : virtual public InferenceTree, virtual public OptimizerTree, public UnrootedTree{

protected:
    /** ************************************************************************
     * OptimizerTreeBasic
     * @semantics  empty constructor to bypass the registration to
     *             Factory< OptimizerTree >
     ************************************************************************ */
    OptimizerTreeBasic();

    /** ************************************************************************
     * OptimizerTreeBasic
     * @semantics  primitive for the constructor of the prototypes used
     *             in ML inference, fully functionnal descendants call
     *             this constructor with their registration name to the
     *             unique Factory< OptimizerTree >
     ************************************************************************ */
    OptimizerTreeBasic( const string & registrationName );

    /** ************************************************************************
     * OptimizerTreeBasic
     * @semantics  constructor of an empty MLTreeBasic according to the
     *             parameters, this constructor should be call by all
     *             descendant
     ************************************************************************ */
    OptimizerTreeBasic( ParametersSet& treeParameters );

    
public:


    /** ************************************************************************
     * ~OptimizerTreeBasic
     * @semantics  destructor of a ML inference tree
     ************************************************************************ */
    ~OptimizerTreeBasic();

    /** ************************************************************************
     * constructFromString
     * @input      the string representation of the tree (Newick)
     * @semantics  override the default behaviour defined in InferenceTree
     *             and unrooted tree in order to construct a tree with MLNodes
     ************************************************************************ */
     virtual void constructFromString( string stringTree );

    /** ************************************************************************
     * diffloglikelihood
     * @semantics  compute the gradient of the likelihood function
     *             at the actual point
     ************************************************************************ */
    virtual void diffloglikelihood(vector < double > & gradientVector);

    /** ************************************************************************
     * loadDataAndModel
     * @input          the sequences and the model to use with this tree
     * @preconditions  no data or model already loaded
     ************************************************************************ */
    virtual void loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel );

};

#endif //OPTIMIZERTREEBASIC_H




