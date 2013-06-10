#ifndef MOLCLKMCMCTREE_H
#define MOLCLKMCMCTREE_H

#include "Tree/InferenceTree.h"
#include "Tree/MCMCTreeBasic.h"
#include "Tree/MolClkTree.h"

#include "Models/PerturbatorHelper.h"

class PerturbatorParameter;

class MolClkMCMCTree : virtual public InferenceTree, public MCMCTreeBasic, public MolClkTree {
protected:
    /** ************************************************************************
     * MolClkMCMCTree
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< MCMCTree >
     ************************************************************************ */
    MolClkMCMCTree( const string & registrationName );

public:
    /** ************************************************************************
     * MolClkMCMCTree
     * @input      parametersSet, parameters transmitted to MolClkTree
     *             and MCMCTreeBasic
     * @semantics  constructor of the class
     ************************************************************************ */
    MolClkMCMCTree(ParametersSet& parameters);

    /** ************************************************************************
     * ~MolClkMCMCTree
     * @semantics  virtual destructor of the class
     ************************************************************************ */
    virtual ~MolClkMCMCTree();

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
     * NNI
     * @input      two descendant of the same node a grandchild of the closer
     *             child and the farer child
     * @semantics  NNI perturbation, the farerChild is swapped with the
     *             child of closerChild, distance to the root are preserved
     ************************************************************************ */
    void NNI( InferenceNode* closerChildChild, InferenceNode* farerChild );

    /** ************************************************************************
     * NNI
     * @input      slideChild, the child of the branch that will move up or down
     * @input      destChild, the child of the destination branch
     * @semantics  second NNI perturbation designed for topology changes,
     *             when clusters or tip dates are defined in the tree
     *             destChild is either (slideChild->parent->parent)
     *             the child of the branch just above the slided one
     *             or one of the child of the two branches below
     *             (slideChild->parent->second_child->child1_or_child2)
     ************************************************************************ */
    void slideNNI( InferenceNode* slideChild, InferenceNode* destChild, double dist );


    /** ************************************************************************
     * sweeping
     * @input      the node to "sweep"
     * @semantics  the subtree defined by the node is detached from the tree
     *             and is reattached between newBrother and its parent (the
     *             parent of node is inserted at distance dist from the parent).
     *             Branch lengths are conserved during the move to ensure
     *             minimal computation, consequently the distance from
     *             node to its parent vary depending on newBrother and dist
     ************************************************************************ */
    void sweeping( InferenceNode* node, InferenceNode* newBrother, double dist );

     /** ************************************************************************
     * constructFromString
     * @input          the string representation of the tree
     * @semantics      initialise the tree with a user-specified string
     * @preconditions  stringTree is a rooted strictly binary tree
     ************************************************************************ */
    virtual void constructFromString( string stringTree );

    
    /** ************************************************************************
     * stopBurn
     * @semantics  end of the burnin
     ************************************************************************ */
    virtual void stopBurn();

    enum{
        SWEEP_PROP,
        NNI_PROP,
        NNI_SLIDE,
        SWAP_PROP, //dirty but useful for the descendant HeterogeneousMolClk
        INVALID_GLOBAL
    };
    unsigned int lastGlobalChange;
    unsigned int numberSweepingPerturbation;
    unsigned int acceptedSweepingPerturbation;
    unsigned int numberNNIPerturbation;
    unsigned int acceptedNNIPerturbation;

    enum{
        LENGTH,
        NNI_LOCAL,
        TREE_HEIGHT,
        INVALID_LOCAL
    };
    unsigned int lastLocalChange;
    unsigned int numberLocalNNI;
    unsigned int acceptedLocalNNI;

    PerturbatorHelper* treeHeightPerturbator;
    PerturbatorHelper* branchLengthPerturbator;

    //some parameters stored to retrieve the old state if the new state is
    //not accepted during a MCMC cycle
    vector< InferenceNode* > pertSave;
    double oldDist;
    InferenceNode* internalNode;

private:
    /** ************************************************************************
     * getCloseFarChildren
     * @semantics  return the two children of internalNode
     ************************************************************************ */
    void getCloseFarChildren( InferenceNode* internalNode,
                              InferenceNode*& closerChild,
                              InferenceNode*& farerChild );



    static MolClkMCMCTree prototype;

};

#endif
