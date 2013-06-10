#ifndef SIMULATIONTREE_H
#define SIMULATIONTREE_H

#include "Tree/Tree.h"

#include "Models/Model.h"

#include "Util/randombox.h"
#include "PatternDesign/Singleton.h"

class SimulationNode;

#include <set>
using namespace std;

class SimulationTree : public Tree {

private:
    SimulationTree();

protected:
    void constructFromString( string stree );

    /* ************************************************************************
     * generate
     * @include     a set of parameters
     * @semantics   construct a random topology according to the rule
     *              specified in parameters
     ************************************************************************ */
    void generate( ParametersSet& parameters );

    /* ************************************************************************
     * generateBetaSplitTopology
     * @include     the beta-splitting parameters
     * @include     numberLeaves, the number of leaves in the tree
     * @include     numberLeaves2, if specified create a bifurcating topology
     *              with numberLeaves species on one side and numberLeaves2
     *              on the other side.
     * @semantics   construct a random topology with the given number of
     *              leaves.
     *              beta = -2->comb
     *              beta = 0 -> pure birth/Yule/random-joining/coalescent
     *              beta = -3/2 -> uniform
     *              beta = -1 -> Aldous special case, arguably better fit
     *              beta = inf -> symetric tries
     ************************************************************************ */
    void generateBetaSplitTopology( double beta, unsigned int numberLeaves,
                                    unsigned int numberLeaves2 = 0 );
    void recursiveConstruct( unsigned int cladeSize,
                             SimulationNode* parent, double beta );

    
    /* ************************************************************************
     * assignDates
     * @include     the node from where to start
     * @include     an ordered set of speciation date, the last (biggest)
     *              value being the age of the given node
     * @semantics   dates of speciation are transformed into branch lengths
     *              leaves are assigned the time t=0
     ************************************************************************ */
    void assignDates( BasicNode* node, const set<double> & speciationDates );
  
    /* ************************************************************************
     * createTopology
     * @include     the desired TopoID
     * @include     the set of parameters to get the parameter necessary
     *              for the execution of topoId
     * @include     numberSpecies/numberOutgroupSpecies, total number
     *              of species and eventually the number of species
     *              in the outgroup
     * @semantics   call generateBetaTopology with the appropriate parameters
     *              to create the tree with a uniform, Yule or beta
     *              distributed topology
     ************************************************************************ */
    void createTopology( unsigned int topoId, ParametersSet& parameters,
                         unsigned int numberSpecies,
                         unsigned int numberOutgroupSpecies=0 );
    
    /* ************************************************************************
     * createBranchLengths
     * @include     the desired BranchId
     * @include     the set of parameters to get the parameter necessary
     *              for the execution of branchId
     * @include     numberSpecies/numberOutgroupSpecies, total number
     *              of species and eventually the number of species
     *              in the outgroup
     * @semantics   initialise the branches according to the specified
     *              distribution (uniform, exponential) or from the selected
     *              process.
     *              Transform time into evolutionary distance if necessary
     *              and scale the tree
     ************************************************************************ */
    void createBranchLengths( unsigned int branchId, ParametersSet& parameters,
                              unsigned int numberSpecies,
                              unsigned int numberOutgroupSpecies=0 );
    
public:
    /* ************************************************************************
     * SimulationTree
     * @include     a set of parameters
     * @semantics   the standard constructor with a set of parameters
     ************************************************************************ */
    SimulationTree( ParametersSet& parameters );

    /* ************************************************************************
     * loadModel
     * @input         model, an initialized model
     * @semanctics    there is no data anymore in a simulation tree since we
     *                have to generate then. However there is still a model.
     ************************************************************************ */
    void loadModel( Model * pmodel  );
    
    /* ************************************************************************
     * SimulationTree
     * @include     the number of symbols for each class of the model
     * @semantics   generate the sequences according to the loaded model
     ************************************************************************ */
    void simulate( const vector<unsigned int> length );
    
    /* ************************************************************************
     * printSequenceTable
     * @input           stream to output the sequences
     * @semantics       output the sequences
     * @preconditions   sequences have been gnerated by a call to simulate()
     ************************************************************************ */
    void printSequenceTable( ostream & out );

private:
    vector< vector<unsigned int> > categories;
    
    enum TopoID{
        TOPO_YULE,
        TOPO_BETA,
        TOPO_UNI,
        NB_TOPO_MODEL
    };
    
    enum BranchID{
        BRANCH_PURE_BIRTH,
        BRANCH_BIRTH_DEATH,
        BRANCH_EXP,
        BRANCH_UNI,
        NB_BRANCH_MODEL
    };
        
    BasicNode* outgroup;
    
public:
    Model * pmodel;

};

#endif //SIMULATIONTREE_H




