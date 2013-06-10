#ifndef HETEROGENEOUSCONSENSUSTREE_H
#define HETEROGENEOUSCONSENSUSTREE_H


class Cluster;
class ClustersSet;
class SequenceTable;
class Heterogeneous;

#include "Tree/ConsensusTree.h"

class HeterogeneousConsensusTree : public ConsensusTree{

public:
    /** ************************************************************************
     * ConsensusTree
     * @input      a registration name to save prototype in the ConsensusTree
     *             factory
     * @semantics  save prototypes of this class and of its descendant in the
     *             factory
     ************************************************************************ */
    HeterogeneousConsensusTree( const string & registrationName );

    /** ************************************************************************
     * HeterogeneousConsensusTree
     * @input      parameters, a ParametersSet containing information about how
     *             the sampling was done
     * @semantics  initialise the consensus tree with the same parameters used
     *             to initialise the sampling
     ************************************************************************ */
    HeterogeneousConsensusTree( ParametersSet& parameters );

    /** ************************************************************************
     * ~HeterogeneousConsensusTree
     * @semantics  virtual destructor
     ************************************************************************ */
    virtual ~HeterogeneousConsensusTree();
        
    /** ************************************************************************
     * clone
     * @input      a ParametersSet used to create a new ConsensusTree
     ************************************************************************ */
    virtual ConsensusTree* clone( ParametersSet& parameters ) const;

     /** ************************************************************************
     * process
     * @return     the number of processed tree
     * @semantics  launch the consensus processus
     ************************************************************************ */
    virtual unsigned int process();
    
    
    
protected:    
    /** ************************************************************************
     * prepareWriting
     * @semantics  once the trees are ready save them to a file and prepare
     *             the consensus file (called by writeConsensusFile()
     ************************************************************************ */
    virtual void prepareWriting();

     /** ************************************************************************
     * readBranchModels
     * @return     the next vector of models in the .bm  file
     * @semantics  parse the next line of branch models in branchModelsFile
     ************************************************************************ */
    void readBranchModels( istream& inputFile,
                             vector<unsigned int> & branchModels,
                             unsigned int expectedNumber );
                             
     /** ************************************************************************
     * readModelParameters
     * @return     the next vector of model parameters in the .mp file, ie
     *             the average substitution rate and the parameters for each
     *             model of the heterogeneous bundle.
     * @semantics  parse the next line of model parameters in modelParametersFile
     ************************************************************************ */     
    void readModelParameters( istream& inputFile,
                              vector<double> & averageSubstitutionRate,
                              vector< vector<double> > & modelParameters );
    
    /** ************************************************************************
     * produceTree
     * @input      a list of clades to build the tree, key = number species
     * @return     build the tree according to the included clades
     * @semantics  called by produce consensus once clades have been selected
     ************************************************************************ */
    virtual void produceTree( multimap< unsigned int, Cluster* > & treeClades );
    
    
    /** ************************************************************************
     * createCluster
     * @input      a node of a sampled tree
     * @return     create an empty cluster for that node
     * @semantics  create the type of cluster specific for the consensus tree
     *             and fill its parameters part.
     ************************************************************************ */
    virtual Cluster* createCluster(BasicNode* node);
            
    
     /** ************************************************************************
     * processModelClades
     * @input      treeClades.begin(), model parameters, branchModels.begin()
     *             sampledTree
     * @return     modified clusters
     * @semantics  add the model parameters for each cluster (add_on to
     *             processClades)
     *             NB: recursive call
     ************************************************************************ */
    void processModelClades( vector< Cluster* >::iterator& cl,
                       const vector < double > & averageSubstitutionRate,
                       const vector < vector<double> > & modelParameters,
                       vector< unsigned int >::iterator& bm, Tree& sampledTree,
                       BasicNode* startPoint = NULL );
    
        
     /** ************************************************************************
     * processModelClades
     * @input      a cluster and the id of a parameter 
     * @return     the mean of the parameter i-1 for the given cluster
     * @semantics  similar to ConsensusTree::getLength, used to return a value
     *             for a parameter (use i==0 for average substitution rate)
     ************************************************************************ */
    double getValue( Cluster* child, unsigned int i );
    
    
    //output file stream where states have been sampled (model for each branch)
    ifstream branchModelsFile;

    //output file stream where substitution model parameters have been
    //sampled (there are more than one model and average substitution rates
    //might differ)
    ifstream modelParametersFile;
        
    Heterogeneous* model;
    
    vector< Tree > parametersTree;
    
    static HeterogeneousConsensusTree prototype;
};

#endif //HETEROGENEOUSCONSENSUSTREE_H
