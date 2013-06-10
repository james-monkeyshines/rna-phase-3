#ifndef MCMCPHASE_H
#define MCMCPHASE_H

#include "Phase.h"

class Model;
class MCMCTree;
class SequenceTable;



class MCMCphase : public Phase {

public:

    /** ***********************************************************************
     * runMCMC
     * @input      startingTree, initialTree
     * @input      startingModel, initialModel
     * @input      perturbationTree, perturbationModel, priority for the tree
     *             and the model
     * @input      burnin and sampling cycles, nuber of cycles in each period
     * @input      sampling period, number of cycles between two samples
     *             likelihood is sampled since the beginning of the run,
     *             sampling of other values starts with the sampling period
     * @input      terminator, a character to close the string representation
     *             of the tree in the best tree file (';' or '')
     * @semantics  launch the MCMC chain
     ************************************************************************ */
    void runMCMC( MCMCTree * startingTree, Model * startingModel,
                  unsigned int perturbationTree, unsigned int perturbationModel,
                  int burninCycles, int samplingCycles, int samplingPeriod,
                  int currentIteration, string terminator );

    /** ***********************************************************************
     * randomTree
     * @input      table, the sequence table used to provide the labels
     * @input      maxBranchLengths, the biggest distance allowed between 2
     *             nodes
     * @input      outgroup, the outgroup used to standardize the tree before
     *             the method returns
     * @semantics  create a random tree, the internal branch are more likely
     *             to be shorter than the terminal ones
     ************************************************************************ */
    InferenceTree* randomTree( SequenceTable * table, double maxBranchLengths,
                               const string& outgroup );

    /** ***********************************************************************
     * run
     * @input          argc, number of arguments in argv[]
     * @input          argv, an array of string arguments
     * @preconditions  argc is equal to 1 and argv[0] is the location
     *                 of a 'control-file' for a MCMC inference
     ************************************************************************ */
    int run( int argc, char * argv[] );

   /* some ofstream for different outputs */
    ofstream plotfile;
    ofstream output;
    ofstream bestTreeFile;
    ofstream bestModelParametersFile;
    ofstream mpFile;
    ofstream priorFile;
    ofstream rescue;
    ofstream rescueBest;

   /* some information about the best state visited by the chain during
    * the run */
    string bestStringTree;
    vector< double > bestModelParams;
    double bestLikelihood;

   /* time structure to record the begining and the end of the simulation */
    time_t start;
    time_t end;

    //store information about the acceptance rates
    unsigned int acceptedModelPerturbation;
    unsigned int numberModelPerturbation;
    unsigned int acceptedTreePerturbation;
    unsigned int numberTreePerturbation;
};
#endif //MCMCPHASE_H




