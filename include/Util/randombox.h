#ifndef RANDOMBOX_H
#define RANDOMBOX_H

#include <math.h>
#include <vector>

#ifndef M_1_SQRT_2
#define M_1_SQRT_2      0.707106781186547524400844362105
#endif

#ifndef M_PI
#define M_PI            3.141592653589793238462643383279
#endif

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

class randombox {
protected:
    // Constructors protected to be sure only Singleton instance can be made
    randombox();
    randombox(const randombox&){};

public:
    void setSeed( long int seed ) ;

    // Methods (ge)
    /* ************************************************************************
     * ran/gasdev/sexp/sgamma/sbeta/sdirichlet
     * @semantics   a set of methods to generate samples from standard
     *              probability dnsity function
     *********************************************************************** */
    double ran() ; // returns uniform random deviate r such that 0.0 < r < 1.0
    double gasdev() ; // returns gaussian random deviate with mean 0 and variance 1
    double sexp(); // returns exponential dist. random deviate
    double sgamma( double , double ) ; // returns gamma dist. random deviate
    double sbeta( double , double ) ; // returns beta dist. random deviate
    void sdirichlet( const vector < double > &, vector < double >& dir ) ;
    
    long int number_of_calls() ;
    
    void getState( vector<double>& params ) const;
    void setState( const vector<double>& params );
    
    inline unsigned int getNumberParams() const{
        return 37;
    }
    
    // Methods query
protected:
    long int calls ;
    long int idnum ;
    long int iy ;
    long int iv[32] ;

    // Needed for sampling from the normal distribution
    int iset ;
    double gset ;

    // Needed for sampling from the gamma distribution
    double aa;
    double aaa;
    double b;
    double c;
    double d;
    double e;
    double p;
    double q;
    double r;
    double s;
    double t;
    double u;
    double v;
    double w;
    double x;
    double q0;
    double s2;
    double si;

  // Constructs
  double a1 , a2 , a3 , a4 , a5 , a6 , a7 ;
  double e1 , e2 , e3 , e4 , e5;
  double q1 , q2 , q3 , q4 , q5 , q6 , q7;
  double sqrt32 ;
 
  double qg[16] ;
};

#endif //RANDOMBOX_H
