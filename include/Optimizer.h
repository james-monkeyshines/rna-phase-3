#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include "Phase.h"

class Optimizer : public Phase{
public:
    int run( int argc, char * * argv );        
};

#endif //OPTIMIZER_H
