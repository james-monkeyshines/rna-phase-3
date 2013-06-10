#ifndef OPTIMIZERTREE_H
#define OPTIMIZERTREE_H

#include "Tree/InferenceTree.h"


/** ****************************************************************************
 * the interface used by the program optimise
 **************************************************************************** */
class ParametersSet;

class OptimizerTree : virtual public InferenceTree{
public:

    /** ************************************************************************
     * ~OptimizerTree
     * @semantics  virtual destructor needed for inheritance
     ************************************************************************ */
    virtual ~OptimizerTree(){};

    /** ***********************************************************************
     * clone
     * @semantics  function used by the Factory<OptimizerTree> to clone the prototype
     *             and create a real model.
    ************************************************************************ */
    virtual OptimizerTree * clone( ParametersSet & parameters ) const = 0;

    /** ************************************************************************
     * initialisation
     * @semantics  initialise the tree (an underlying call to loadDataAndModel
     *             is performed)
     * @return     true if there is an initial tree
     ************************************************************************ */
    virtual bool initialisation( SequenceTable * ptable, Model * pmodel ) = 0;

    /** ************************************************************************
     * diffloglikelihood
     * @semantics  compute the gradient of the likelihood function
     *             at the actual point
     ************************************************************************ */
    virtual void diffloglikelihood(vector < double > & gradientVector) = 0;

};

#endif //OPTIMIZERTREE_H
