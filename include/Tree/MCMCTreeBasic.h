#ifndef MCMCTREEBASIC_H
#define MCMCTREEBASIC_H


#include "Tree/MCMCTree.h"
#include "Tree/InferenceTree.h"
#include "PatternDesign/Singleton.h"
#include "Util/randombox.h"


class SequenceTable;
class Perturbator;

class MCMCTreeBasic : virtual public InferenceTree, public MCMCTree {
protected:
    /** ************************************************************************
     * MCMCTreeBasic
     * @semantics  primitive for the constructor of the prototypes used
     *             in MCMC inference, fully functionnal descendants call
     *             this constructor with their registration name to the
     *             unique Factory< MCMCTree >
     ************************************************************************ */
    MCMCTreeBasic( const string & registrationName );

    /** ************************************************************************
     * MCMCTreeBasic
     * @semantics  constructor of an empty MCMCTreeBasic according to the
     *             parameters, this constructor should be call by all
     *             descendant
     ************************************************************************ */
    MCMCTreeBasic( ParametersSet& treeParameters );

public:

   /** ************************************************************************
     * ~MCMCTree
     * @semantics  destructor of an unrooted inference tree
     ************************************************************************ */
    ~MCMCTreeBasic();



    /** ************************************************************************
     * initialiseMCMC
     * @semantics   initialise the MCMC system (priors and proposals)
     ************************************************************************ */
    virtual void initialiseMCMC( ParametersSet& parameters );

    /** ************************************************************************
     * perturb
     * @return      ln(Hasting ratio) of the perturbation
     * @semantics   perturb the tree
     ************************************************************************ */
    virtual double perturb();

    /** ************************************************************************
     * initSampling
     * @input       a set of parameters, gives the name for the files
     *              that store branch lengths and topologies and specify
     *              the character to put at the end of each tree
     *              (';' : PHYLIP, '' : BAMBE)
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
     * @semantics   output the parameters of the pertubation
     ************************************************************************ */
    virtual void printPerturbationParameters(ostream& outputStream);


    /** ************************************************************************
     * globalTopologyChange
     * @return      ln(Hasting ratio) of the perturbation
     * @semantics   perturb the topology of the tree
     ************************************************************************ */
    virtual double globalTopologyChange() = 0;

    /** ************************************************************************
     * changeBranchLength
     * @return      ln(Hasting ratio) of the perturbation
     * @semantics   perturb the length of one edge
     ************************************************************************ */
    virtual double changeBranchLength() = 0;

     /** ************************************************************************
     * validateTopologyChange
     * @input       true if the validation is accepted, false to rollback
     * @return      true
     * @semantics   validate the change in the tree topology
     ************************************************************************ */
    virtual bool validateTopologyChange( bool validation ) = 0;

    /** ************************************************************************
     * validateBranchLength
     * @input       true if the validation is accepted, false to rollback
     * @return      true
     * @semantics   validate the new length for the edge
     ************************************************************************ */
    virtual bool validateBranchLength( bool validation ) = 0;

    /** ************************************************************************
     * validatePerturbation
     * @input      boolean, yes if the perturbation is accepted
     * @semantics  return the tree to its initial state if the perturbation
     *             was not accepted
     ************************************************************************ */
    virtual bool validatePerturbation( bool validation );


    /** ************************************************************************
     * stopBurn
     * @semantics  end of the burnin
     ************************************************************************ */
    virtual void stopBurn();


    /** ************************************************************************
     * update
     * @input       an update message from the observed model
     ************************************************************************* */
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

protected:
    Singleton<randombox> & randBox;

    //InferenceTree copy;

    double treePriority;
    double branchPriority;
    double hyperPriority;

    //output file stream to sample states
    ofstream samplesFile;
    ofstream branchLengthsFile;
    ofstream hyperParametersFile;
    string terminator;

    /*necessary for eventual hyperpriors */
    Perturbator* hyperPerturbator;

private:
    enum{
        TOPOLOGY,
        LOCAL,
        HYPER_PRIORS,
        INVALID
    };
    unsigned int lastPerturbation;

};

#endif //MCMCTREEBASIC_H
