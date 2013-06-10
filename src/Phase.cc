#include "Phase.h"
#include <locale.h>

Phase::Phase() {
}

int Phase::run( int, char *[] ) {
    //to avoid some conflict with international systems
    //(when using static build and other softs)
    //output (and input) files in phase will respect a common standard
    setlocale(LC_ALL, "POSIX");
    return EXIT_SUCCESS;
}
