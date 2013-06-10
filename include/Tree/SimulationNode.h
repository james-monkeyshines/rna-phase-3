#ifndef SIMULATIONNODE_H
#define SIMULATIONNODE_H

#include "Tree/BasicNode.h"

class SimulationTree;

#include <assert.h>

#include <vector>

using namespace std;

class SimulationNode : public BasicNode{
public:
    /** ***********************************************************************
     * SimulationNode
     * constructors
     *********************************************************************** */
    SimulationNode( SimulationTree* ptree );
    SimulationNode( const string& label, SimulationTree* ptree );
    SimulationNode( const TreeMap& map, SimulationTree* ptree );
    
    /** ***********************************************************************
     * simulate
     * @input        a vector of the mixture category (usually the gamma
     *               category) assumed for each site (one vector per model)
     * @input        the state generated at the parent (one vector of state
     *               per model)
     * @semantics    generate a sequence for this node and transmit
     *               the recursive call to children
     *********************************************************************** */
    void simulate( const vector< vector<unsigned int> > & categories,
                   const vector< vector<unsigned int> > & parent );
                   
    /** ***********************************************************************
     * getSequence
     * @input          modelId, the desired model
     * @return         the generated sequences for the model at this node
     * @preconditions  simulate was successfully executed before
     *********************************************************************** */
    inline const vector<unsigned int> & getSequence(unsigned int modelId) const{
        assert(modelId<sequence.size());
        return sequence[modelId];
    }
    
protected:
    SimulationTree* ptree;
    vector< vector<unsigned int> > sequence;
};

#endif  //SIMULATIONNODE_H
