#ifndef OPTIMISE_H
#define OPTIMISE_H

#include <math.h>
#include <vector>

class OptimizerTree;
class Model;


using namespace std;


class Optimise{

public:
    /** ************************************************************************
     * optimiseQuasiNewton
     * @input      tree, the phylogenetic tree
     * @input      model, the model
     * @input      modelOpt, whether the model parameters should be initialised as well
     * @input      empFreqs, whether empirical frequencies should be used
     * @input      tol, tolerance used as a step size
     * @input      display, yes to print likelihood,... at each iter
     * @input      displayPoint, yes to print a '.' on the screen at each iteration
     * @return     a vector with the optimized values for branches and
     *             optionnaly the model
     * @semantics  Optimise the branches (and optionally the model parameters)
     *             of a phylogenetic tree.
     ************************************************************************ */
    static vector<double> optimiseQuasiNewton( OptimizerTree & tree, Model & model,
    bool modelOpt, bool empFreqs, double tolerance, bool display, bool displayPoint = true );

  
private:

    /** ************************************************************************
     * diffLoglikelihoodSqrt
     * @input          tree, the phylogenetic tree
     * @input          branchesLengthSqrt, the square root of the current
     *                 branch lengths
     * @inout          a gradient vector to fill
     * @return         the partial derivatives (for each sqrt(branch length))
     *                 at the current point in the likelihood curve
     * @semantics      use the difflikelihood function of the OptimizerTree to compute
     *                 the gradient for each branch, branchLengthSqrt is used
     *                 to perform the conversion of the gradient
     *                 dLikelihood/dsqrt(branch) =
     *                 2 sqrt(branch) * dLikelihood/dbranch
     * @preconditions  branchesLengthSqrt is equal to
     *                 sqrt(tree->getBranchVector())
     * @preconditions  branchGradientVector.size == branchesLengthSqrt.size
     ************************************************************************ */
     static void diffLoglikelihoodSqrt( OptimizerTree & tree ,
                                  const vector<double> & branchesLengthSqrt,
                                  vector<double> & branchGradientVector);

    /** ************************************************************************
     * modelPartialGradients
     * @input          tree, the phylogenetic tree
     * @input          model, the substitution model
     * @input          modelParameters, the parameters of this substitution
     *                 model. The function will modify it internally but it
     *                 will be restored before leaving
     * @input          optParameters, the indices of the modelParameters
     *                 vector that should be optimised - eg if empirical
     *                 frequencies are used, then those parameters are not
     *                 optimised.
     * @input          step, the small step used to compute the gradient
     *                 for each parameter
     * @inout          a gradient vector to fill
     * @semantics      Compute numerically the partial derivatives of the
     *                 likelihood for the model parameters. The parameters of
     *                 the current model are modified during execution but
     *                 actual parameters are restored on exit.
     * @preconditions  modelParameters == model.getAllParameters()
     *                 the function ask for this redondant vector as input
     *                 to avoid a memory allocation
     * @preconditions  gradientVector.size == modelParameters.size
     ************************************************************************ */
    static void modelPartialGradients( OptimizerTree & tree, Model & model,
                 vector<double> & modelParameters,
                 vector<unsigned int> &optParameters, double step,
                 vector<double> &gradientVector );
                 
    /** ************************************************************************
     * modelPenaltyPartialGradients
     * @input          model, the substitution model
     * @input          modelPenalty, the parameters of the penalized approach
     * @input          step, the small step used to compute the gradient
     *                 for each parameter
     * @inout          a gradient vector to fill
     * @semantics      Compute numerically the partial derivatives of the
     *                 likelihood for the free penalty parameters. The parameters
     *                 the current model are modified during execution but
     *                 actual parameters are restored on exit.
     * @preconditions  modelPenalty == model.getAllPenaltyParameters()
     *                 the function ask for this redondant vector as input
     *                 to avoid a memory allocation
     * @preconditions  gradientVector.size == modelPenalty.size
     ************************************************************************ */
    static void modelPenaltyPartialGradients( Model & model,
        vector < double > & modelPenalty, double step,
        vector<double> &gradientVector );
    
    
    /** ************************************************************************
     * logOptValue
     * @input      tree, the phylogenetic tree
     * @input      sqrtInitialBranches, square root of the branches of tree
     * @input      newBranches, junk vector
     * @input      initialModelPenalty, free parameters of the penalized
     *             approach
     * @input      newModelPenalty, junk vector
     * @input      the direction we are aiming too
     * @input      magnitude, the factor used to modify the starting point
     *             towards the given direction
     * @return     the loglikelihood of the tree with the new branches length
     *             defined by [ old branches ] + magnitude*direction
     * @semantics  BEWARE! the tree is modified during the call
     *             and TREE IS NOT RESTORED (idem for the penalization parameters)
     *             This function is called during the lineSearch, it returns
     *             the value we are trying to optimise
     *             (ie logLikelihood+logPenalty) at the point defined by:
     *             current point + magnitude * direction_vector
     * @preconditions  sqrtInitialBranches is equal to
     *                 sqrt(tree->getBranchVector())
     ************************************************************************ */
    static double logOptValue( OptimizerTree & tree, Model & model,
       const vector<double>& sqrtInitialBranches, vector<double>& newBranches,
       const vector<double>& initialModelPenalty, vector<double>& newModelPenalty,  
       vector<double>& direction, double magnitude );
    
    /** ************************************************************************
     * logOptValue
     * @input      tree, the phylogenetic tree
     * @input      model, the model of evolution used
     * @input      sqrtInitialBranches, square root of the branches of tree
     * @input      newBranches, junk vector
     * @input      initialModelPenalty, free parameters of the penalized
     *             approach
     * @input      newModelPenalty, junk vector
     * @input      initialParameters, parameters of model
     * @input      newParameters, junk vector
     * @input      the direction we are aiming too
     * @input      magnitude, the factor used to modify the starting point
     *             towards the given direction
     * @return     the loglikelihood of the tree with the new branches length
     *             and new model parameters defined by
     *             [ old branches, old parameters ] + magnitude*direction
     * @semantics  BEWARE! TREE and MODEL are modified during the call
     *             and they are NOT RESTORED (idem for the penalization parmaeters)
     * @preconditions  sqrtInitialBranches is equal to
     *                 sqrt(tree->getBranchVector())
     * @preconditions  initialParameters == model.getAllParameters()
     *                 the function ask for this redondant vector as input
     *                 to avoid a memory allocation
     ************************************************************************ */
    static double logOptValue( OptimizerTree & tree, Model & model,
       const vector<double>& sqrtInitialBranches, vector<double>& newBranches,
       const vector<double>& initialModelPenalty, vector<double>& newModelPenalty,  
       const vector<double>& initialParameters, vector<double>& newParameters,
       vector<double>& direction, double magnitude );

    /** ************************************************************************
     * lineSearch
     * @input          tree, the phylogenetic tree
     * @input          model, the model of evolution used
     * @input          branchVector, the current branch lengths in tree
     * @input          sqrtBranchVector, the square root of the current
     *                 branch lengths in tree
     * @input          modelVector, the current parameters in the model
     *                 (size can be 0 if modelOpt = FALSE)
     * @input          direction of the line search
     * @input          the opposite of the gradient vector
     *                 (used to compute the slope)
     * @input          modelOpt, true if the model parameters are in the search
     *                 space
     * @semantics      the approximate line search algorithm to be used with the
     *                 quasi-newton optimizations
     * @preconditions  branchVector == tree.getBranchVector() and
     *                 sqrtBranchVector == sqrt(branchVector)
     * @preconditions  if modelOpt, modelVector ==  model.getAllParameters()
     * @preconditions  direction.size() == branchVector.size() +
     *                   ( modelOpt ? modelVector.size() : 0 )
     ************************************************************************ */
    static double lineSearch( OptimizerTree & tree, Model & model,
        const vector<double> & branchVector,
        const vector<double> & sqrtBranchVector,
        const vector < double > & modelPenaltyVector,
        const vector<double> & modelVector,
        vector<double> & direction, vector<double> & oppGradient, bool modelOpt);

};

#endif //OPTIMISE_H




