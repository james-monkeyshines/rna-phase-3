#include "Tree/HeterogeneousMolClkConsensusTree.h"

#include "Sequence/SequenceTable.h"

#include "Tree/Cluster.h"
#include "Tree/LengthCluster.h"
#include "Tree/HeterogeneousCluster.h"

//Convention used in consensus : the name must be the name of the 
//corresponding MCMC tree with consensus instead of MCMC
HeterogeneousMolClkConsensusTree HeterogeneousMolClkConsensusTree::prototype
("Heterogeneous consensus tree with molecular clock");

HeterogeneousMolClkConsensusTree::HeterogeneousMolClkConsensusTree( const string & registrationName ):
HeterogeneousConsensusTree(registrationName){};

HeterogeneousMolClkConsensusTree::HeterogeneousMolClkConsensusTree( ParametersSet& parameters ):
HeterogeneousConsensusTree(parameters){};


ConsensusTree* HeterogeneousMolClkConsensusTree::clone( ParametersSet& parameters ) const{
    return new HeterogeneousMolClkConsensusTree(parameters);
}

Cluster* HeterogeneousMolClkConsensusTree::createCluster(BasicNode* node){
    double distanceLeaf = 0.0;
    while(!node->isLeaf()){
        node = node->getChild(0);
        distanceLeaf += node->getParentDistance();
    }
    return new HeterogeneousCluster( seq->getNumberSpecies(),
                              distanceLeaf );
}

double HeterogeneousMolClkConsensusTree::getDist( Cluster* father, Cluster* child ){
    return ( ((LengthCluster*)father)->getLength()/
                  (double)((LengthCluster*)father)->getNumber()  -
             ((LengthCluster*)child)->getLength()/
                  (double)((LengthCluster*)child)->getNumber() );
}
