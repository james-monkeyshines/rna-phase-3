#ifndef CONSENSUSTREE_H
#define CONSENSUSTREE_H

#include "Util/ParametersSet.h"

#include "Tree/Tree.h"

#include <vector>
#include <map>
#include <fstream>
using namespace std;

class Cluster;
class ClustersSet;
class BasicNode;
class SequenceTable;

#define BLFILE_LINELENGTH 20000
#define MPFILE_LINELENGTH 20000

class ConsensusTree : public Tree{

public:
    /** ************************************************************************
     * ConsensusTree
     * @input      a registration name to save prototype in the ConsensusTree
     *             factory
     * @semantics  save prototypes of this class and of its descendant in the
     *             factory
     ************************************************************************ */
    ConsensusTree( const string & registrationName );

    /** ************************************************************************
     * ConsensusTree
     * @input      parameters, a ParametersSet containing information about how
     *             the sampling was done
     * @semantics  initialise the consensus tree with the same parameters used
     *             to initialise the sampling
     ************************************************************************ */
    ConsensusTree( ParametersSet& parameters );

    /** ************************************************************************
     * ~ConsensusTree
     * @semantics  virtual destructor
     ************************************************************************ */
    virtual ~ConsensusTree();

    /** ************************************************************************
     * clone
     * @input      a ParametersSet used to create a new ConsensusTree
     ************************************************************************ */
    virtual ConsensusTree* clone( ParametersSet& parameters ) const;

     /** ************************************************************************
     * process
     * @return     the number of processed tree
     * @semantics  create the sets of cluster found at least one in all the
     *             samples
     ************************************************************************ */
    virtual unsigned int process();

    /** ************************************************************************
     * produceConsensus
     * @return     the consensus tree (topology + branch lengths) + files
     * @semantics  perform the consensus with all the clades sampled
     ************************************************************************ */
    virtual void produceConsensus();

    /** ************************************************************************
     * writeConsensusFile
     * @semantics  output the results in the consensus file (trees are in
     *             another file)
     ************************************************************************ */
    virtual void writeConsensusFile();
    
protected: 
        
    /** ************************************************************************
     * prepareWriting
     * @semantics  once the trees are ready save them to a file and prepare
     *             the consensus file (called by writeConsensusFile()
     ************************************************************************ */
    virtual void prepareWriting();
    
    /** ************************************************************************
     * produceTree
     * @input      a list of clades to build the tree, key = number species
     * @return     build the tree according to the included clades
     * @semantics  called by produce consensus once clades have been selected
     ************************************************************************ */
    virtual void produceTree( multimap< unsigned int, Cluster* > & treeClades );


    /** ************************************************************************
     * readTree
     * @return     the next string tree in samples file
     * @semantics  read the next tree, EXIT if something goes wrong
     ************************************************************************ */
    void readTree( istream& inputFile, string & stringTree );

    /** ************************************************************************
     * readBranchLengths
     * @return     the next vector of lengths in the .bl  file
     * @semantics  parse the next line of branch lengths in branchLengthsFile
     ************************************************************************ */
    void readBranchLengths( istream& inputFile,
                             vector<double> & branchLengths,
                             unsigned int expectedNumber );
    
    /** ************************************************************************
     * processClades
     * @input      the sampled tree topology
     * @return     a vector of cluster (specific to the consensus tree)
     * @semantics  retrieve the clades in one sampled tree
     *             NB: recursive call
     ************************************************************************ */
    virtual Cluster* processClades( const Tree& sampledTree,
                                vector< Cluster* > & treeClades,
                                BasicNode* startPoint = NULL );

    
    
    /** ************************************************************************
     * createCluster
     * @input      a node of a sampled tree
     * @return     create an empty cluster for that node
     * @semantics  create the type of cluster specific for the consensus tree
     *             and fill its parameters part.
     ************************************************************************ */
    virtual Cluster* createCluster(BasicNode* node);
        
    
    /** ************************************************************************
     * getDist
     * @input      two clusters
     * @return     the distance to put between father and child
     * @semantics  each kind of cluster contains information to set the length
     *             in the final consensus tree. getDist is virtual too allow
     *             modification by descendant class
     ************************************************************************ */
    virtual double getDist( Cluster* father, Cluster* child );
    
    
    //output file stream where states have been sampled
    ifstream samplesFile;
    ifstream branchLengthsFile;

    ofstream consensusFile;
    ofstream consensusTreeFile;
    
    ParametersSet consensusFileParameters;
    
    //the number of samples according to the control file
    unsigned int expectedNumberSamples;
    unsigned int foundNumberSamples;
    //a pointer to the sequence table (to match number<->species' names)
    SequenceTable* seq;
    
    ClustersSet* clades;
    
    Tree supportTree;
    
    static ConsensusTree prototype;
};;

#endif //CONSENSUSTREE_H
