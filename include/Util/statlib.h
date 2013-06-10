/**
 * Some algorithms included here are NOT covered by PHASE license
 * They are copyrighted by the Royal Statistical Society,
 * permission for distribution is granted provided that no fee is charged
 */

#ifndef STATLIB_H
#define STATLIB_H

#include <vector>


#ifndef M_1_SQRT_2
#define M_1_SQRT_2      0.707106781186547524400844362105
#endif
/*
#ifndef M_PI
#define M_PI            3.141592653589793238462643383279
#endif
*/
#ifndef M_PI_half
#define M_PI_half       1.570796326794896619231321691640
#endif

#ifndef M_SQRT_PI
#define M_SQRT_PI       1.772453850905516027298167483341
#endif

#ifndef M_1_SQRT_2PI
#define M_1_SQRT_2PI    0.398942280401432677939946059934
#endif


#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI   0.918938533204672741780329736406
#endif

using namespace std;

class statlib{
public:
    static double ldirichlet(const vector <double> &, const vector <double> &);
    static vector<unsigned int> random_permutation(const vector <unsigned int> &);
    static void gamma_percentile_means(vector <double> &,
            double alpha , double beta , int nsections);
    static double ppchi2(double p , double v , int &errorcode);
    static double gammad(double x , double p , int &ifault) ;
    static double ppnd(double p, int & errorcode);
    static double alnorm(double x, bool upper);
    static void gammaLaguerre(vector < double > & rates, vector < double > & weigths, double alpha,
        double beta, int nsections);
private:
    static double rootLaguerre( double approxRoot, unsigned int order,
                              double alpha, double& dValueLx,
                              double& valueLy,
                              const double* recB, const double* recC );
};

#endif //STATLIB_H
