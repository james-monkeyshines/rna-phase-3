#include "Util/randombox.h"

#include <assert.h>
#include <math.h>
#include <cstdlib>

#define repeat for(;;)
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.0e-14
#define RNMX (1.0-EPS)


//to be removed
#include <iostream>

// Constructor (Initializes random number generator with seed ==1)
randombox::randombox() {
    int j;
    long int k;

    idnum = 1;
    // Load shuffle table after 8 warmups
    for (j = NTAB + 7; j >= 0; j--) {
        k = (idnum / IQ);
        idnum = IA * (idnum - k * IQ) - IR * k;
        if (idnum <= 0) idnum += IM;
        if (j < NTAB) iv[j] = idnum;
    }
    iy = iv[0];

    calls = 0;
    iset = 0;
    aa = aaa = 0.0;
    qg[0] = 0.6931471805599453;
    qg[1] = 0.9333736875190459;
    qg[2] = 0.9888777961838675;
    qg[3] = 0.9984959252914960;
    qg[4] = 0.9998292811061389;
    qg[5] = 0.9999833164100727;
    qg[6] = 0.9999985691438767;
    qg[7] = 0.9999998906925558;
    qg[8] = 0.9999999924734159;
    qg[9] = 0.9999999995283275;
    qg[10] = 0.9999999999728814;
    qg[11] = 0.9999999999985598;
    qg[12] = 0.9999999999999289;
    qg[13] = 0.9999999999999968;
    qg[14] = 0.9999999999999999;
    qg[15] = 1.0000000000000000;

    a1 = 0.3333333;
    a2 = -0.250003;
    a3 = 0.2000062;
    a4 = -0.1662921;
    a5 = 0.1423657;
    a6 = -0.1367177;
    a7 = 0.1233795;
    e1 = 1.0;
    e2 = 0.4999897;
    e3 = 0.166829;
    e4 = 0.0407753;
    e5 = 0.010293;
    q1 = 0.04166669;
    q2 = 0.02083148;
    q3 = 0.00801191;
    q4 = 0.00144121;
    q5 = -7.388e-5;
    q6 = 2.4511e-4;
    q7 = 2.424e-4;
    sqrt32 = 5.656854;
}



// Initializes random number generator with seed
void randombox::setSeed(long int seed) {
    int j;
    long int k;

    idnum = abs(seed);
    if (idnum == 0) idnum = 1;
    // Load shuffle table after 8 warmups
    for (j = NTAB + 7; j >= 0; j--) {
        k = (idnum / IQ);
        idnum = IA * (idnum - k * IQ) - IR * k;
        if (idnum <= 0) idnum += IM;
        if (j < NTAB) iv[j] = idnum;
    }
    iy = iv[0];
}


long int randombox::number_of_calls() {
    return calls;
}

// Return uniformly distributed random deviate
double randombox::ran() {
    int j;
    long int k;
    double temp;
    k = idnum / IQ;
    idnum = IA * (idnum - k * IQ) - IR * k;
    if (idnum <= 0) idnum += IM;
    j = iy / NDIV;
    iy = iv[j];
    iv[j] = idnum;

    calls++;
    temp = AM * iy;
    if (temp > RNMX) return RNMX;
    return temp;
}

// Return normally distributed random deviate
double randombox::gasdev() {
    double fac, rsq, v1, v2;
    ++calls;
    if (iset == 0) {
        do {
            v1 = (2.0 * ran()) - 1.0;
            v2 = (2.0 * ran()) - 1.0;
            rsq = (v1 * v1) + (v2 * v2);
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return (v2 * fac);
    }
    else {
        iset = 0;
        return gset;
    }
}


// return exponentialy distributed random deviate
double randombox::sexp() {
    double a, u, ustar, umin;
    int i;

    a = 0.0;
    u = ran();
    for ( ; ; ) {
        u = u + u;
        if (u > 1.0) break;
        a = a + qg[0];
    }
    u = u - 1.0;

    if (u <= qg[0]) return a + u;

    i = 0;
    ustar = ran();
    umin = ustar;
    do {
        ustar = ran();
        if (ustar < umin) umin = ustar;
        i = i + 1;
    }
    while (u > qg[i]);
    return a + umin * qg[0];
}


double randombox::sgamma(double a, double scale) {
    double ret_val;

    if (a < 1.0) {
    /* alternate method for parameters a below 1 */

    /* 0.36787944117144232159 = exp(-1) */

        aa = 0.0;
        b = 1.0 + 0.36787944117144232159 * a;
        repeat {
            p = b * ran();
            if (p >= 1.0) {
                ret_val = -log((b - p) / a);
                if (sexp() >= (1.0 - a) * log(ret_val)) break;
            } else {
                ret_val = exp(log(p) / a);
                if (sexp() >= ret_val) break;
            }
        }
        return (scale * ret_val);
    }

  /* Step 1: Recalculations of s2, s, d if a has changed */

    if (a != aa) {
        aa = a;
        s2 = a - 0.5;
        s = sqrt(s2);
        d = sqrt32 - s * 12.0;
    }

  /* Step 2: t = standard normal deviate, */

  /* x = (s,1/2)-normal deviate. */

  /* immediate acceptance (i) */

    t = gasdev();
    x = s + 0.5 * t;
    ret_val = x * x;
    if (t >= 0.0) return (scale * ret_val);

  /* Step 3: u = 0,1 - uniform sample. squeeze acceptance (s) */

    u = ran();
    if (d * u <= t * t * t) {
        return (scale * ret_val);
    }

  /* Step 4: recalculations of q0, b, si, c if necessary */

    if (a != aaa) {
        aaa = a;
        r = 1.0 / a;
        q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r + q2) * r + q1) * r;

    /* Approximation depending on size of parameter a */

    /* The constants in the expressions for b, si and */

    /* c were established by numerical experiments */

        if (a <= 3.686) {
            b = 0.463 + s + 0.178 * s2;
            si = 1.235;
            c = 0.195 / s - 0.079 + 0.16 * s;
        } else if (a <= 13.022) {
            b = 1.654 + 0.0076 * s2;
            si = 1.68 / s + 0.275;
            c = 0.062 / s + 0.024;
        } else {
            b = 1.77;
            si = 0.75;
            c = 0.1515 / s;
        }
    }

  /* Step 5: no quotient test if x not positive */

    if (x > 0.0) {
    /* Step 6: calculation of v and quotient q */

        v = t / (s + s);
        if (fabs(v) <= 0.25) q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
        else q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);


    /* Step 7: quotient acceptance (q) */

        if (log(1.0 - u) <= q) return (scale * ret_val);
    }

  /* Step 8: e = standard exponential deviate */

  /* u= 0,1 -uniform deviate */

  /* t=(b,si)-double exponential (laplace) sample */

    repeat {
        e = sexp();
        u = ran();
        u = u + u - 1.0;
        if (u < 0.0) t = b - si * e;
        else t = b + si * e;

    /* Step  9:  rejection if t < tau(1) = -0.71874483771719 */

        if (t >= -0.71874483771719) {
      /* Step 10:  calculation of v and quotient q */

            v = t / (s + s);
            if (fabs(v) <= 0.25)
                q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
            else q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);

      /* Step 11:  hat acceptance (h) */

      /* (if q not positive go to step 8) */

            if (q > 0.0) {
                if (q <= 0.5) w = ((((e5 * q + e4) * q + e3) * q + e2) * q + e1) * q;
                else w = exp(q) - 1.0;

	/* if t is rejected */

	/* sample again at step 8 */

                if (c * fabs(u) <= w * exp(e - 0.5 * t * t)) break;
            }
        }
    }
    x = s + 0.5 * t;
    return scale * x * x;
}


// Sample from a beta distribution
double randombox::sbeta( double a, double b ) {
    double x1 = sgamma( a, 1.0);
    double x2 = sgamma( b, 1.0);
    return x1 / (x1+x2);
}


// Samples from an n-dimensional Dirichlet distribution
void randombox::sdirichlet(const vector < double > & parameters, vector<double>& dir) {
    double gammas[parameters.size()];
    dir.resize(parameters.size());

    double sum = 0.0;
    for (unsigned int i = 0; i < parameters.size(); ++i) {
        gammas[i] = sgamma(parameters[i], 1.0);
        sum += gammas[i];
    }
    for (unsigned int i = 0; i < parameters.size(); ++i) {
        dir[i] = gammas[i] / sum;
    }
}

void randombox::getState( vector<double>& params ) const{
    params.resize(37);
    for (unsigned int i = 0; i < 32; ++i){
        params[i] = (double)iv[i];
    }
    params[32] = (double)idnum;
    params[33] = (double)iy;
    params[34] = (double)iset;
    params[35] = gset;
    params[36] = (double)calls;
}

void randombox::setState( const vector<double>& params ){
    assert (params.size() == 37);
    for (unsigned int i = 0; i < 32; ++i){
        iv[i] = (long int)params[i];
    }
    idnum = (long int)params[32];
    iy = (long int)params[33];
    iset = (long int)params[34];
    gset = params[35];
    calls = (long int)params[36];
}
