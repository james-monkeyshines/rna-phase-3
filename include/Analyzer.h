#ifndef ANALYZER_H
#define ANALYZER_H

#include "Phase.h"


class Analyzer : public Phase {

public:


    /** ***********************************************************************
     * run
     * @input          argc, number of arguments in argv[]
     * @input          argv, an array of string arguments
     * @return         the return of this program, 0 if no error
     * @preconditions  argc is equal to 1 and argv[0] is the location
     *                 of a 'control-file' for a MCMC inference
     ************************************************************************ */
    int run( int argc, char * argv[] );

};
#endif //ANALYZER_H




