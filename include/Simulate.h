#ifndef SIMULATE_H
#define SIMULATE_H

#include "Phase.h"

class SimulationTree;

class Simulate: public Phase {
public:
    int run( int argc, char * argv[] );
};

#endif //SIMULATE_H



