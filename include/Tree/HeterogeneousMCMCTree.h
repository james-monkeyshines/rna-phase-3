#ifndef HETEROGENEOUSMCMCTREE_H
#define HETEROGENEOUSMCMCTREE_H

#include "Tree/InferenceTree.h"
#include "Tree/MCMCTreeBasic.h"
#include "Tree/RootedTree.h"

#include "Models/Heterogeneous.h"
#include "Models/PerturbatorHelper.h"


#include <vector>
#include <map>

class HeterogeneousMCMCTree : virtual public InferenceTree, public MCMCTreeBasic, public RootedTree {
private:
    /** ************************************************************************
     * HeterogeneousMCMCTree
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< MCMCTree >
     ************************************************************************ */
    HeterogeneousMCMCTree( const string & registrationName );

public:
    /** ************************************************************************
     * HeterogeneousMCMCTree
     * @semantics  constructor of the class, u
     ************************************************************************ */
    HeterogeneousMCMCTree(ParametersSet& parameters);

    /** ************************************************************************
     * ~HeterogeneousMCMCTree
     * @semantics  virtual destructor of the class
     ************************************************************************ */
    virtual ~HeterogeneousMCMCTree();

    /** ************************************************************************
     * clone
     * @semantics  the clone method for the factory<MCMCTree>
     ************************************************************************ */
    virtual MCMCTree* clone(ParametersSet& parameters) const;

    /** ************************************************************************
     * globalTopologyChange
     * @return      ln(Hasting ratio) of the perturbation
     * @semantics   perturb the topology of the tree
     ************************************************************************ */
    virtual double globalTopologyChange();

    /** ************************************************************************
     * changeBranchLength
     * @return      ln(Hasting ratio) of the perturbation (ie, 0.0)
     * @semantics   perturb the length of one edge
     ************************************************************************ */
    virtual double changeBranchLength();

     /** ************************************************************************
     * validateTopologyChange
     * @input       true if the validation is accepted, false to rollback
     * @return      true
     * @semantics   validate the change in the tree topology
     ************************************************************************ */
    virtual bool validateTopologyChange( bool validation );

    /** ************************************************************************
     * validateBranchLength
     * @input       true if the validation is accepted, false to rollback
     * @return      true
     * @semantics   validate the new length for the edge
     ************************************************************************ */
    virtual bool validateBranchLength( bool validation );

    /** ************************************************************************
     * initialiseMCMC
     * @input       parameters to initialise the perturbator
     * @semantics   initialise the perturbator with the fields in parameters
     ************************************************************************ */
    virtual void initialiseMCMC(ParametersSet& parameters);

    /** ************************************************************************
     * initSampling
     * @input       a set of parameters, gives the name for the files,
     *              see MCMCTreeBasic::initSampling( ParametersSet& ).
     *              contains the base name used for the file that stores
     *              the model id of each branch
     * @input       overwrite, rewrite or append (in case of restore) ?
     * @semantics   initialise what is necessary to do the sample
     ************************************************************************ */
    virtual void initSampling(ParametersSet& parameters, bool overwrite);

    /** ************************************************************************
     * sample
     * @semantics   store the current state on files
     ************************************************************************ */
    virtual void sample();

    /** ************************************************************************
     * printPerturbationParameters
     * @input       the output stream to output the parameters
     * @semantics   after a MCMC run, output some stats about the changes
     *              during the sampling period
     ************************************************************************ */
    virtual void printPerturbationParameters(ostream& outputStream);

    /** ************************************************************************
     * getAllPerturbationParameters
     * @semantics      get all the parameters used in the perturbator
     *                 to save its state
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    virtual void getAllPerturbationParameters( vector<double>& params ) const;

    /** ************************************************************************
     * setAllPerturbationParameters
     * @semantics      restore the state of the perturbator
     ************************************************************************ */
    virtual void setAllPerturbationParameters( const vector<double>& params);

    /** ************************************************************************
     * getNumberPerturbationParameters
     * @semantics      return the number of parameters used in the perturbator
     * @preconditions  the perturbator should have been initialised
     ************************************************************************ */
    virtual unsigned int getNumberPerturbationParameters() const;

    /** ************************************************************************
     * getAllPriorParameters
     * @semantics      get all the parameters used in the priors
     *                 to save its state
     * @preconditions  the prior should have been initialised
     ************************************************************************ */
    virtual void getAllPriorParameters( vector<double>& params ) const;

    /** ************************************************************************
     * setAllPriorParameters
     * @semantics      restore the state of the prior
     ************************************************************************ */
    virtual void setAllPriorParameters( const vector<double>& params);

    /** ************************************************************************
     * getNumberPriorParameters
     * @semantics      return the number of free parameters used in the prior
     * @preconditions  the prior should have been initialised
     ************************************************************************ */
    virtual unsigned int getNumberPriorParameters() const;

    /** ************************************************************************
     * getLnPrior
     * @input      none
     * @semantics  this function returns ln(prior) according to the user prior
     *             specification (for MCMC runs)
     ************************************************************************ */
    virtual double getLnPrior() const;

    /** ************************************************************************
     * chooseNNISwapNode
     * @input       a node, not the leaf, not the root
     * @semantics   choose a child of this node and choose a node adjacent to
     *              its parent in order to swap them on a future NNI.
     *              the two selected nodes are pushed into pertSave
     ************************************************************************ */
    void chooseNNISwapNode( InferenceNode* nodeChild );

    /** ************************************************************************
     * NNI
     * @input       a branch ID adjacent with two INTERNAL nodes
     * @semantics   NNI perturbation, two nodes linked to this internal edge
     *              are swapped
     ************************************************************************ */
    void NNI( InferenceNode* node1, InferenceNode* node2,
                InferenceNode* swapNode1, InferenceNode* swapNode2 );

    /** ************************************************************************
     * insertBranch
     * @input       branchMovedId, the ID of the branch to move
     * @input       branchInsertId, the ID of the branch where the branch moved
     *              is to be inserted
     * @input       frac, a double between 0.0 and 1.0 to specify the insertion
     *              point in branchInsertId
     * @return      length(branchInsertId)/length(length(2 branches linked to
     *              movedId before the move), this is the hasting ratio
     * @semantics   SPR perturbation, detach
     ************************************************************************ */
    double insertBranch( InferenceNode* insert1, InferenceNode* insert2,
                         InferenceNode* moved1, InferenceNode* moved2,
                         double frac );


    /** ************************************************************************
     * stopBurn
     * @semantics  end of the burnin
     ************************************************************************ */
    virtual void stopBurn();

    /** ************************************************************************
     * loadDataAndModel
     * @input          the sequences and the model to use with this tree
     * @precondition   no data or model already loaded
     * @semantics      this method is redefined (from InferenceTree).
     *                 The aim is to check the model is an heterogeneous one
     *                 and to initiate the extra cooperation between tree
     *                 and model in the heterogeneous framework.
     ************************************************************************ */
    virtual void loadDataAndModel( SequenceTable * psequenceTable, Model * pmodel );

     /** ************************************************************************
     * invalidateNodes
     * @input       cat, the category (cf. InferenceTree::invalidateNodes)
     * @input       node, a node information because sometimes not all the
     *              nodes are to be invalidated when the model is perturbed
     * @semantics   switch to the working space for the category cat and for
     *              the nodes designed by node. Currently, node is the value
     *              of the heterogeneous model modified
     ************************************************************************ */
    //virtual void invalidateNodesAux( int cat, int node );

    /** ***********************************************************************
     * update
     * @input     message from the substitution model
     ************************************************************************ */
     virtual void update( UpdateMessage* message );


     /** ************************************************************************
     * retrieveNodes
     * @input       cat, the category (cf. InferenceTree::retrieveNodes)
     * @input       node, a node information because sometimes not all the
     *              nodes are to be invalidated when the model is perturbed
     * @semantics   come back to the last saved computation for the the
     *              category cat for the nodes designed by
     *              node. Currently, node is the value of the heterogeneous
     *              model modified
     ************************************************************************ */
    virtual void retrieveNodesAux();

    /** ************************************************************************
     * saveNodes
     * @input       cat, the category (cf. InferenceTree::saveNodes)
     * @input       node, a node information because sometimes not all the
     *              nodes are to be invalidated when the model is perturbed
     * @semantics   save the computation for the category cat for the nodes
     *              designed by node. Currently, node is the value of the
     *              heterogeneous model modified
     ************************************************************************ */
    virtual void saveNodesAux();


    /** ************************************************************************
     * getAllParameters
     * @semantics    return all the parameters of the tree (in order to
     *               save/restore the state)
     *               at that stage of the hierarchy it means basically the
     *               branch lengths
     ************************************************************************ */
    virtual void getAllParameters( vector<double>& params ) const;

    /** ************************************************************************
     * setAllParameters
     * @input        the parameters to restore the state
     * @semantics    restore all the parameters of the tree with
     *               the given parameters
     ************************************************************************ */
    virtual void setAllParameters( const vector<double>& params);

    /** ************************************************************************
     * getNumberTreeParameters
     * @input      none
     * @semantics  return the number of free parameters in the tree (branch
     *             lengths, model for each branch, ...)
     ************************************************************************ */
     virtual unsigned int getNumberTreeParameters() const;

    /** ************************************************************************
     * initialiseNodeModel
     * @semantics  once the heterogeneous model has been loaded and the tree
     *             constructed, assign a random model to each node
     ************************************************************************ */
    void initialiseNodeModel();

    /** ************************************************************************
     * constructRandomly
     * @input      the maximum length for the branches during the creation
     * @semantics  create a random initial tree for the species in table
     *             redefined to initialise the model for each node
     ************************************************************************ */
    virtual void constructRandomly( double maxBranchLength );

     /** ************************************************************************
     * constructFromString
     * @input          the string representation of the tree
     * @semantics      initialise the tree with a user-specified string
     *                 redefined to initialise the model for each node
     ************************************************************************ */
    virtual void constructFromString( string stringTree );




    enum{
        SPR_PROP,
        NNI_PROP,
        SWAP_PROP,
        INVALID_GLOBAL
    };
    unsigned int lastGlobalChange;
    unsigned int numberSPRPerturbation;
    unsigned int acceptedSPRPerturbation;
    unsigned int numberNNIPerturbation;
    unsigned int acceptedNNIPerturbation;

    enum{
        LENGTH,
        LOCAL_NNI,
        INVALID_LOCAL
    };
    unsigned int lastLocalChange;
    int modifiedBranchId;
    unsigned int numberLocalNNI;
    unsigned int acceptedLocalNNI;

    //output file stream to sample states
    ofstream branchModelsFile;

    //some parameters stored to retrieve the old state if the new state is
    //not accepted during a MCMC cycle
    vector< InferenceNode* > pertSave;
    double oldFrac;

private:
    static HeterogeneousMCMCTree prototype;
    map< Model*, int> modelMap;
    vector< pair<InferenceNode*, Model*> > modelSave;

    PerturbatorBase* swapPerturbator;
    PerturbatorHelper* branchPerturbator;

    Heterogeneous * pheterogeneous;
    unsigned int numberModels;
};

#endif //HETEROGENEOUSMCMCTREE_H
