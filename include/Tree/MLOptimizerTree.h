#ifndef MLOPTIMIZERTREE_H
#define MLOPTIMIZERTREE_H

#include "Tree/InferenceTree.h"
#include "Tree/OptimizerTreeBasic.h"
#include "Tree/UnrootedTree.h"

#include "Models/PerturbatorHelper.h"

/** ****************************************************************************
 * The tree used in optimise.
 **************************************************************************** */
class MLOptimizerTree : virtual public InferenceTree, public OptimizerTreeBasic{
protected:
    /** ************************************************************************
     * MLOptimizerTree
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< MLTree >
     ************************************************************************ */
    MLOptimizerTree( const string & registrationName );

public:
    /** ************************************************************************
     * MLOptimizerTree
     * @semantics  constructor of the class
     ************************************************************************ */
    MLOptimizerTree(ParametersSet& parameters);

    /** ************************************************************************
     * ~MLOptimizerTree
     * @semantics  virtual destructor of the class
     ************************************************************************ */
    virtual ~MLOptimizerTree();

    /** ************************************************************************
     * clone
     * @semantics  the clone method for the factory<MLTree>
     ************************************************************************ */
    virtual OptimizerTree* clone(ParametersSet& parameters) const;

    /** ************************************************************************
     * initialisation
     * @semantics  initialise the tree (an underlying call to loadDataAndModel
     *             is performed)
     * @return     true if there is an initial tree
     ************************************************************************ */
    virtual bool initialisation( SequenceTable * ptable, Model * pmodel );

private:
    static MLOptimizerTree prototype;

    string optimizedTree;

};

#endif //MLOPTIMISERTREE_H
