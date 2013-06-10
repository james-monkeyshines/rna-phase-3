#ifndef CONFIGFIX
#define CONFIGFIX

//Mac OS X users
#if defined(__APPLE__)
#include <cmath>
    extern "C" int isnan (double);
    extern "C" int isinf (double);
#endif

//Sun
#if defined (__sun)
#include <ieeefp.h>
#include <math.h>
    inline int isinf(double x) { return !finite(x) && x==x; }
#endif

#endif
