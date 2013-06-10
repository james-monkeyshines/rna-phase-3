#ifndef UNROOTEDMCMCTREE_H
#define UNROOTEDMCMCTREE_H

#include "Tree/InferenceTree.h"
#include "Tree/MCMCTreeBasic.h"
#include "Tree/UnrootedTree.h"

#include "Models/PerturbatorHelper.h"

class UnrootedMCMCTree : virtual public InferenceTree, public MCMCTreeBasic, public UnrootedTree {
private:
    /** ************************************************************************
     * UnrootedMCMCTree
     * @semantics  constructor of the class for the prototype, called once,
     *             with the name used to register to the unique (Singleton)
     *             Factory< MCMCTree >
     ************************************************************************ */
    UnrootedMCMCTree( const string & registrationName );
    
public:
    /** ************************************************************************
     * UnrootedMCMCTree
     * @semantics  constructor of the class, u
     ************************************************************************ */
    UnrootedMCMCTree(ParametersSet& parameters);
                
    /** ************************************************************************
     * ~UnrootedMCMCTree
     * @semantics  virtual destructor of the class
     ************************************************************************ */
    virtual ~UnrootedMCMCTree();
    
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
                            InferenceNode* swapNode1, InferenceNode* swapNode2 );;
            
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
     * constructFromString
     * @input          the string representation of the tree
     * @semantics      initialise the tree with a user-specified string
     * @preconditions  one leaf label in stringTree is the outgroup
     * @preconditions  stringTree is a strictly binary tree (rooted or not)
     ************************************************************************ */
    virtual void constructFromString( string stringTree );
            
    /** ************************************************************************
     * stopBurn
     * @semantics  end of the burnin
     ************************************************************************ */
    virtual void stopBurn();    
    
        
    enum{
        SPR_PROP,
        NNI_PROP,
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
    double maxBranchLength;
    unsigned int numberLocalNNI;
    unsigned int acceptedLocalNNI;
    
    //some parameters stored to retrieve the old state if the new state is
    //not accepted during a MCMC cycle
    vector< InferenceNode* > pertSave;
    double oldFrac;
    
private:
    static UnrootedMCMCTree prototype;

    PerturbatorHelper* branchPerturbator;            
};

#endif //UNROOTEDMCMCTREE_H
