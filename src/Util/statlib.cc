/**
 * Some algorithms included here are NOT covered by PHASE license
 * They are copyrighted by the Royal Statistical Society,
 * permission for distribution is granted provided that no fee is charged
 */
 
#include "Util/statlib.h"

#include <math.h>


using namespace std;

#include <iostream>

#include "PatternDesign/Singleton.h"
#include "Util/randombox.h"

#ifndef _REENTRANT
double lgamma_r(double a, int * sgamma) {
    * sgamma = 1;
    return lgamma(a);;
}

#endif

extern "C" {
    double dstev_( char * jobz, int* n, double* d,
                    double* e, double* z,int* ldz, double* work, int* info);
}




// returns the log dirichlet probability density given a set of n values and
// n dirichlet parameters
double statlib::ldirichlet(const vector < double > & values, const vector < double > & parameter) {
    double sum1, sum2, sum3;
    int sgamma;

    sum1 = sum2 = sum3 = 0.0;
    for (unsigned int i = 0; i < values.size(); i++) {
        sum1 += parameter[i];
        sum2 += lgamma_r(parameter[i], & sgamma);
        sum3 += ((parameter[i] - 1.0) * log(values[i]));
    }

    return (sum3 - sum2 + lgamma_r(sum1, & sgamma));
};


vector < unsigned int > statlib::random_permutation(const vector < unsigned int > & thearray) {

    Singleton<randombox>& box = Singleton<randombox>::instance();

    vector <unsigned int> ans = thearray;

    for (vector<unsigned int>::iterator i = ans.begin(); i != ans.end(); ++i) {
        vector<unsigned int>::iterator select = ans.begin() + (unsigned int)(box.ran() * (double)thearray.size());
        swap( *i, *select );
    }

    return ans;
};

void statlib::gamma_percentile_means(vector < double > & ans, double alpha, double beta, int nsections) {
    double sector[nsections - 1];
    double lower, upper, size = 1.0 / (double)nsections;
    int fault;

    ans.resize(nsections);
    //fill the sector vector with the cutting point.
    for (int i = 0; i < nsections - 1; i++) {
        sector[i] = ppchi2((i + 1.0) * size, 2.0 * alpha, fault);
        sector[i] /= (2.0 * beta);
    }

    //use the sector vector to retrieve the means of each category
    upper = gammad(sector[0] * beta, alpha + 1.0, fault);
    ans[0] = ((alpha / beta) * upper * nsections);

    for (int i = 1; i < nsections - 1; i++) {
        fault = -1;
        lower = upper;
        upper = gammad(sector[i] * beta, alpha + 1.0, fault);
        ans[i] = ((alpha / beta) * (upper - lower) * (double)nsections);
    }

    ans[nsections - 1] = ((alpha / beta) * (1.0 - upper) * (double)nsections);
}

// Finds the value which corresponds to the p*100th percentile of a chi-squared distribution
// with v degrees of freedom
double statlib::ppchi2(double p, double v, int & errorcode) {
//Algorithm AS 91   Appl. Statist. (1975) Vol.24, P.35
    int sgamma;
    double ppnd(double, int &); // Calculates the cutting points of the standard normal distribution
    // double gammad(double , double , int &) ; // Calculates the incomplete gamma ratio given x and p

    double g = lgamma_r(v / 2.0, & sgamma), pp;
    int maxit = 20, i, if1;

    double aa = 0.6931471806, e = 0.5e-6;
    double zero = 0.0, half = 0.5, one = 1.0, two = 2.0, three = 3.0, six = 6.0;
    double pmin = 0.000002, pmax = 0.999998;

    double c[38] = { 0.01, 0.222222 , 0.32 , 0.4 , 1.24, 2.2,
                   4.67, 6.66, 6.73, 13.32, 60.0, 70.0,
                   84.0, 105.0, 120.0, 127.0, 140.0, 175.0,
                  210.0, 252.0, 264.0, 294.0, 346.0, 420.0,
                  462.0, 606.0, 672.0, 707.0, 735.0, 889.0,
                  932.0, 966.0, 1141.0, 1182.0, 1278.0, 1740.0 , 2520.0, 5040.0};
    double s[6];
    double a, b, c1, ch, p1, p2, q, t, x, xx;

    g = sgamma * g;

    // Test arguments and initialize
    pp = -one;
    if ((p < pmin) || (p > pmax)) return 1;
    if (v <= zero) return 2;

    xx = half * v;
    c1 = xx - one;

    // Starting approximation for small chi-squared
    if (!(v >= (-c[4] * log(p)))) {
        ch = pow((p * xx * exp(g + xx * aa)), (one / xx));
        if (ch < e) { // Goto 6
            pp = ch;
            errorcode = 0;
            return pp;
        }
        // Goto 4
    }
    else { // Goto 1

        if (!(v > c[2])) {
            ch = c[3];
            a = log(one - p);
            do { // 2
                q = ch;
                p1 = one + ch * (c[6] + ch);
                p2 = ch * (c[8] + ch * (c[7] + ch));
                t = -half + (c[6] + two * ch) / p1 - (c[8] + ch * (c[9] + three * ch)) / p2;
                ch = ch - (one - exp(a + g + half * ch + c1 * aa) * p2 / p1) / t;
            } while (fabs(q / ch - one) > c[0]);
        }
        else { // 3
            x = statlib::ppnd(p, if1);

            // starting approximation using Wilson and Hilferty estimate
            p1 = c[1] / v;
            ch = v * pow((x * sqrt(p1) + one - p1), 3.0);

            // starting approximation for p tending to 1

            if (ch > c[5] * v + six) ch = -two * (log(one - p) - c1 * log(half * ch) + g);

        }
    }

    // 4
    //  Call to algorithm AS 239 and calculation of seven term Taylor series
    i = 0;
    do {
        q = ch;
        p1 = half * ch;
        p2 = p - gammad(p1, xx, if1);
        if (if1 != 0) {
            cerr << "Gamma Error = " << if1 << endl;
            errorcode = 3;
            return 0.0;
        }
        t = p2 * exp(xx * aa + g + p1 - c1 * log(ch));
        b = t / ch;
        a = half * t - b * c1;
        s[0] = (c[18] + a * (c[16] + a * (c[13] + a * (c[12] + a * (c[11] + c[10] * a))))) / c[23];
        s[1] = (c[23] + a * (c[28] + a * (c[31] + a * (c[32] + c[34] * a)))) / c[36];
        s[2] = (c[18] + a * (c[24] + a * (c[27] + c[30] * a))) / c[36];
        s[3] = (c[19] + a * (c[26] + c[33] * a) + c1 * (c[21] + a * (c[29] + c[35] * a))) / c[37];
        s[4] = (c[12] + c[20] * a + c1 * (c[17] + c[25] * a)) / c[36];
        s[5] = (c[14] + c1 * (c[22] + c[15] * c1)) / c[37];
        ch = ch + t * (one + half * t * s[0] - b * c1 * (s[0] - b * (s[1] - b * (s[2] - b * (s[3] - b * (s[4] - b * s[5]))))));

        if (fabs(q / ch - one) > e) {
            pp = ch;
            errorcode = 0;
            return (pp);
        }
        i++;
    } while (i < maxit);

    errorcode = 4;
    pp = ch;
    return (pp);
}


// Finds the lower tail end which corresponds to the p*100th percentile of the standard normal
// distribution
double statlib::ppnd(double p, int & errorcode) {
//ALGORITHM AS 111, APPL.STATIST., VOL.26, 118-121, 1977.
    double split = 0.42, cp;
    double q, r;
    double a[4] = { 2.50662823884 , -18.61500062529 , 41.39119773534 , -25.44106049637};
    double b[4] = {-8.47351093090 ,  23.08336743743 ,-21.06224101826 , 3.13082909833};
    double c[4] = {-2.78718931138 , -2.29796479134  , 4.85014127135  , 2.32121276858};
    double d[2] = {3.54388924762  , 1.63706781897};

    q = p - 0.5;
    errorcode = 0;

    if (!(fabs(q) > split)) {
        r = q * q;
        cp = q * (((a[3] * r + a[2]) * r + a[1]) * r + a[0]) / ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0);
        return cp;
    }
    r = p;
    if (q > 0.0) r = 1.0 - p;

    if (r <= 0.0) {
        cp = 0.0;
        errorcode = 1;
        return cp;
    }

    r = sqrt(-log(r));
    cp = (((c[3] * r + c[2]) * r + c[1]) * r + c[0]) / ((d[1] * r + d[0]) * r + 1.0);
    if (q < 0.0) cp = -cp;
    return (cp);


}


// calculates the incomplete gamma ratio given x and p
double statlib::gammad(double x, double p, int & ifault) {
//algorithm AS239  APPL. STATIST. (1988) VOL. 37, NO. 3
    double alnorm(double, bool);
    int sgamma;
    double result;
    double pn[6];
    double tol = 1.0e-14, oflo = 1.0e37;
    double xbig = 1.0e8, arg, c, rn, a, b, one = 1.0, zero = 0.0, an, two = 2.0;
    double elimit = -88.0, plimit = 1000.0, three = 3.0, nine = 9.0;

    result = 0.0;

    if (p <= 0.0 || x < 0.0) {
        ifault = 1;
        return result;
    }

    ifault = 0;
    if (x == 0.0) return result;

    // Use a normal appromimation if p > plimit
    if (p > plimit) {
        pn[0] = three * sqrt(p) * (pow((x / p), (one / three)) + one / (nine * p) - one);
        result = statlib::alnorm(pn[0], false);
        return result;
    }

    // If x is extremely large compared to p then return 1
    if (x > xbig) {
        result = one;
        return (result);
    }

    // Pearson's series expansion
    if ((x <= one) || (x < p)) {
        arg = p * log(x) - x - lgamma_r(p + one, & sgamma);
        c = one;
        result = one;
        a = p;
        do { // 40
            a += one;
            c = c * x / a;
            result = result + c;
        } while (c > tol);
        arg = arg + log(result);
        result = zero;

        if (arg >= elimit) {
            // elmit exceeded
            result = exp(arg);
        };
    }
    else {

        // Use continued fraction
        arg = p * log(x) - x - lgamma_r(p, & sgamma);
        a = one - p;
        b = a + x + one;
        c = zero;
        pn[0] = one;
        pn[1] = x;
        pn[2] = x + one;
        pn[3] = x * b;
        result = pn[2] / pn[3];
        do { // 60
            //  cout << "looping 60" << endl ;
            a += one;
            b += two;
            c += one;
            an = a * c;
            pn[4] = b * pn[2] - an * pn[0];
            pn[5] = b * pn[3] - an * pn[1];
            if (fabs(pn[5]) > zero) {
                rn = pn[4] / pn[5];
                if ( (fabs(result - rn) < tol) &&
                     (fabs(result - rn) < (tol * rn) ) ) {
                    // 80
                    arg += log(result);
                    result = 1.0;

                    if (arg >= elimit) result = one - exp(arg);
                    return (result);

                }
                result = rn;
            }
            pn[0] = pn[2];
            pn[1] = pn[3];
            pn[2] = pn[4];
            pn[3] = pn[5];
            if (fabs(pn[4]) >= oflo) { // Rescale
                pn[0] = pn[0] / oflo;
                pn[1] = pn[1] / oflo;
                pn[2] = pn[2] / oflo;
                pn[3] = pn[3] / oflo;
            }
        } while (true); // End do
    }


    return (result);

}


double statlib::alnorm(double x, bool upper) {
//algorithm AS66  Applied Statistics (1973) vol22 no.3
    double alnorm1, zero = 0.0, one = 1.0, half = 0.5;
    double con = 1.28, ltone = 7.0, utzero = 18.66, z, y;
    double p = 0.398942280444, q = 0.39990348504, r = 0.398942280385;
    double a[3] = {5.75885480458  , 2.62433121679 , 5.92885724438};
    double b[2] = {-29.8213557807 ,48.6959930692};
    double c[6] = {-3.8052e-8 , 3.98064794e-4 , -0.151679116635 , 4.8385912808 , 0.742380924027 , 3.99019417011};
    double d[5] = {1.00000615302 , 1.98615381364 , 5.29330324926 , -15.1508972451 , 30.789933034};
    bool up = upper;

    z = x;

    if (z < zero) {
        up = !up;
        z = -z;
    }
    // 10
    if (!((z <= ltone) || up && z <= utzero)) {

        alnorm1 = zero;
        if (!up) alnorm1 = one - alnorm1;
        return (alnorm1);
    }

    // 20
    y = half * z * z;
    if (z > con) {
        // 30
        alnorm1 = r * exp(-y) / (z + c[0] + d[0] /
           (z + c[1] + d[1] / (z + c[2] + d[2] / (z + c[3] + d[3] / (z + c[4] + d[4] / (z + c[5]))))));
        // 40
        if (!up) alnorm1 = one - alnorm1;
        return (alnorm1);
    }
    alnorm1 = half - z * (p - q * y / (y + a[0] + b[0] / (y + a[1] + b[1] / (y + a[2]))));

    // 40
    if (!up) alnorm1 = one - alnorm1;
    return (alnorm1);

}

void statlib::gammaLaguerre(vector < double > & rates, vector < double > & weigths, double alpha,
        double beta, int nsections){

    rates.resize( nsections );
    weigths.resize( nsections );

    //coef of the Laguerre polynomial
    double a[nsections];
    double b[nsections];

    //w(x) = (beta*x)^alpha * exp (-beta*x)
    //apha parameter for the Laguerre polynome and Yang's gamma distribution
    //do not match perfectly     (alpha_Laguerre = alpha - 1)

    //recursive coef for the laguerre polynomial
    for (int i = 0; i< nsections; ++i){
        a[i] = alpha + (double)( 2 * i );
    }
    for (int i = 0; i< nsections-1; ++i){
        b[i] = -sqrt((i+1)*(i+alpha));
    }

//    cc = exp(lgamma_r( alpha + 1.0, &mom )) * product;
    //find the eigenvalues and first components
    //of the eigen vector of the symetric tridiagonal matrix
    //on input rates are the diagonal element and b the subdiagonal elements
    // weights is [1 0 0 0 0 0 0 ...]_n
    // on output, rates are the eigen values in ascending order
    //            weights are the first components of the orthonormal eigen
    //            vectors

    char V = 'V';
    int info;
    double z[nsections][nsections];
    double workspace[2*nsections-2];
    dstev_( &V, &nsections, a, b, (double*)z, &nsections, workspace, &info);

//    double tot1=0.0;
//    double tot2=0.0;
    for ( int i = 0; i < nsections; ++i){
        rates[i] = a[i] / beta;
        weigths[i] = z[i][0] *  z[i][0];
//        tot1 += rates[i] * weigths[i];
//        tot2 += (rates[i]-alpha/beta)  * (rates[i]-alpha/beta) * weigths[i];
    }
//    cerr << "average rate calc " << tot1 << " vs " << alpha/beta  << endl;
//    cerr << "variance calc " << tot2 << " vs " << alpha/(beta*beta)  << endl;
}
