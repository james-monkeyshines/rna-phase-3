#ifndef PHASE_H
#define PHASE_H

#include <math.h>
#include <stdio.h>
#include <iomanip>
#include <time.h>
// #include <sunmath.h>
// #include <sunperf.h>

// Molecular sequence table
#include "Sequence/SequenceTable.h"

// substitution model prototypes
#include "Models/Model.h"

// Phylogenetic Tree
#include "Tree/InferenceTree.h"

// Statistics library
#include "Util/statlib.h"



class Phase {
public:

    Phase();

    int run( int argc, char * argv[] );
};


#endif //PHASE_H




