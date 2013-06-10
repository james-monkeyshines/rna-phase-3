#include "Tree/MolClkConsensusTree.h"

#include "Sequence/SequenceTable.h"

#include "Tree/Cluster.h"
#include "Tree/LengthCluster.h"

//Convention used in consensus : the name must be the name of the 
//corresponding MCMC tree with consensus instead of MCMC
MolClkConsensusTree MolClkConsensusTree::prototype("Rooted consensus tree with molecular clock");

MolClkConsensusTree::MolClkConsensusTree( const string & registrationName ):
ConsensusTree(registrationName){};

MolClkConsensusTree::MolClkConsensusTree( ParametersSet& parameters ):
ConsensusTree(parameters){};


ConsensusTree* MolClkConsensusTree::clone( ParametersSet& parameters ) const{
    return new MolClkConsensusTree(parameters);
}
        
        
Cluster* MolClkConsensusTree::createCluster(BasicNode* node){
    double distanceLeaf = 0.0;
    while(!node->isLeaf()){
        node = node->getChild(0);
        distanceLeaf += node->getParentDistance();
    }
    return new LengthCluster( seq->getNumberSpecies(),
                              distanceLeaf );
}

double MolClkConsensusTree::getDist( Cluster* father, Cluster* child ){
    return ( ((LengthCluster*)father)->getLength()/
                  (double)((LengthCluster*)father)->getNumber()  -
             ((LengthCluster*)child)->getLength()/
                  (double)((LengthCluster*)child)->getNumber() );
}
