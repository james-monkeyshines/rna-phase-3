#include "Util/matrixmath.h"

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <float.h>
#include <cstdio>

#include <algorithm>
#include <numeric>

#include <functional>

using namespace std;


extern "C" {
  void MAIN__(){};
}

// Lapack svd routine
 extern "C" {
    int dgesvd_(char * , char * , int * , int *,  double * , int * , double * ,
                double * , int * , double * , int * , double * , int * , int *) ;
  }

// Lapack eigensystem routine
 extern "C" {
  int dgeev_(char *jobvl, char *jobvr, int *n, double *a, int *lda, double *wr, double *wi,
              double *vl, int *ldvl, double *vr, int *ldvr, double *work, int *lwork, int *info);
  }

// Lapack LU inversion
 extern "C" {
  int dgetri_( int *n, double *a, int* lda, int* ipiv, double* work, int* lwork, int *info);
  }
 extern "C" {
  int dgetrf_(int *m, int *n, double *a, int* lda, int* ipiv, int *info);
  }

 extern "C" {
  int ilaenv_(int *ispec, char *name, char* opts, int* n1, int* n2, int *n3, int *n4);
  }

array2D<double> transpose(array2D<double> &matrix) {
  int nrows = matrix.numberRows() , ncol = matrix.numberColumns() ;
  int i , j ;
  array2D<double>ans(ncol , nrows) ;

  for (i=0 ; i<nrows ; i++)
    for (j=0 ; j<ncol ; j++)
      ans(j,i) = matrix(i,j) ;

  return ans ;
}

// Matrix * scalar
array2D<double> operator* (const array2D<double> & matrix , double num) {
  int rows1 = matrix.numberRows() , col1 = matrix.numberColumns() ;
  int i , j;

  array2D<double>ans(rows1 , col1) ;

  for (i=0 ; i<rows1 ; i++)
      for (j=0 ; j<col1 ; j++)
          ans(i,j) = matrix(i,j) * num ;

  return ans ;
}

// Matrix * scalar
array2D<double> operator* (double num , const array2D<double> & matrix) {
  int rows1 = matrix.numberRows() , col1 = matrix.numberColumns() ;
  int i , j;

  array2D<double>ans(rows1 , col1) ;

  for (i=0 ; i<rows1 ; i++)
    for (j=0 ; j<col1 ; j++)
       ans(i,j) = matrix(i,j) * num ;

  return ans ;
}

// Matrix * matrix
array2D<double> operator* (const array2D<double> & matrix1 , const array2D<double> & matrix2) {
  unsigned int rows1 = matrix1.numberRows() , col1 = matrix1.numberColumns() ;
  unsigned int col2 = matrix2.numberColumns() ;
  unsigned int i , j , k ;
  double sum ;
  array2D<double> ans(rows1 , col2) ;

#ifdef DEBUG1
    assert(col1==matrix2.numberRows());
#endif

  for (i=0 ; i<rows1 ; i++)
    for (k=0 ; k<col2 ; k++) {
      sum = 0 ;
      for (j=0 ; j<col1 ; j++)
	sum += matrix1(i,j) * matrix2(j , k) ;
      ans(i,k) = sum ;
    }

  return ans ;
}

array2D<double> outer_product(vector <double> &a, vector <double> &b) {

    array2D<double> ans(a.size() , b.size()) ;

    for (unsigned int i=0 ; i < a.size() ; ++i){
        for (unsigned int j=0 ; j < b.size() ; ++j){
            ans(i,j) = a[i] * b[j];
        }
    }
    return (ans) ;
}



// Finds the inverse of a square matrix i.e N*N using svd (Singular value decomposition
array2D<double> inverse(array2D<double> &matrix) {
    // SVD routine from LAPACK

    int i , j , rows = matrix.numberRows() , cols = matrix.numberColumns() ;
    array2D<double> A = transpose(matrix);
    array2D<double> U(rows , cols) ;
    array2D<double> VT(rows, cols) ;
    array2D<double> S(rows , cols) ;

    // A =  matrix ;
    for (i=0 ; i<rows ; i++){
        for (j=0 ; j<cols ; j++) {
            U(i,j) = VT(i,j) = S(i,j) = 0.0 ;
        }
    }

    int Lwork = rows * 5 ;  // Workspace size for the svd routine
    int info ;
    double * Work = new double[Lwork] ;  // Workspace for the svd routine
    double * diag = new double[rows] ;
    dgesvd_("A" , "A" , &rows , &rows , &(A(0,0)) , &rows , diag , &(U(0,0)) , &rows , &(VT(0,0)) , &rows , Work , &Lwork , &info) ;
    assert(!info);


    // U has effectively already been transposed since lapack treats matrices as column major
    // ditto with V**T
    for (i=0 ; i<rows ; i++){
        for (j=0 ; j<rows ; j++){
            S(i,j) = 1/diag[i] * U(i,j) ;
        }
    }
    delete[] diag ;
    delete[] Work ;
    return (VT * S) ;
}

//inverse a SYMMETRIC definite positive matrix using LU decomposition, return its det
double inverseSPD(array2D<double> &matrix, int order){
    int info;
    double logdet = 0.0;
    int ipiv[order];

    dgetrf_( &order, &order, &(matrix(0,0)), &order, ipiv, &info);
    if(info){
        assert(info>=1);
        fprintf(stderr,"WARNING, badly conditionned matrix in the gaussian process (U(i,i)=%f)\n",matrix(info-1,info-1));
        matrix(info-1,info-1) = 0.0;
        logdet = -DBL_MAX;
    }

    int count = 1;
    //compute the determinant ( sign[ipiv] * diagprod(matrix) )
    for (unsigned int i = 0; i< (unsigned int)order; ++i){
        if (ipiv[i]!=(int)i+1){
            count=-count;
        }
        if (matrix(i,i) < 0){
            count = -count;
        }
        logdet+=log(fabs(matrix(i,i)));
    }
    if(count!=1){
        fprintf(stderr,"WARNING, badly conditionned matrix in the gaussian process (computation of det(C))\n");
        logdet = -DBL_MAX;
    }
    int one = 1;
    int minusone = -1;
    char opts[2];
    opts[0] = ' ';
    opts[1] = 0;
    char dgetriName[7] = "DGETRI";
    int optSize = ilaenv_( &one, dgetriName, opts, &order, &minusone, &minusone, &minusone );
    optSize*=order;
    double workspace[optSize];
    dgetri_( &order, &(matrix(0,0)), &order, ipiv, workspace, &optSize, &info);
    assert(!info);
    return logdet;
}


//inverse a SYMMETRIC definite positive matrix using SVD, return its det
double inverseSPD2(array2D<double> &matrix, int order){

    vector<double> junk(order);
    array2D<double> junkmatrix(order,order);
    vector<double> eigenValues(order);

    // Find the eigenValues and eigenVectors of the rate matrix
    RightEigenSystem( matrix, junkmatrix, eigenValues, junk );
    array2D<double> ret = inverse(matrix);
    matrix = ret;
    std::transform(eigenValues.begin(), eigenValues.end(), eigenValues.begin(), (double (*) (double))fabs);
    std::transform(eigenValues.begin(), eigenValues.end(), eigenValues.begin(), (double (*) (double))log);
    return std::accumulate<vector<double>::iterator, double>(eigenValues.begin(), eigenValues.end(), 0);
}





// Returns true if the matrix given is a unit square matrix
bool unit(array2D<double> &A) {
  int size = A.numberRows() , i , j ;

  for (i=0 ; i<size ; i++) {
    if ((A(i,i) != 1.0))
      return false ;
    for (j=i+1 ; j<size ; j++) {
      if (A(i,j) != 0.0)
	return false;
    }
  }

  return true ;
}

// Determines the "left" eigen decomposition of a general matrix
int LeftEigenSystem(array2D<double> &A , array2D<double>  &Evectors ,  vector <double> &Reigenvalues,
                    vector <double> &Ieigenvalues ) {
    int N=A.numberRows();
    char rvecs = 'N' ;  // Don't get the right eigen vectors
    char lvecs = 'V' ;  // Do calculate the left eigen vectors
    int LWork=24*N;
    if (N>400)
       LWork=4*N;
    double Work[LWork];
    int rinfo , LDVL=1;

    array2D<double> A_t = transpose(A);

    dgeev_(&lvecs, &rvecs, &N, &(A_t(0,0)), &N, &(Reigenvalues[0]) , &(Ieigenvalues[0]), &(Evectors(0,0)), &N,
	   NULL, &LDVL, Work, &LWork, &rinfo);

    // left eigen vectors are stored as the ROWS of the eigen matrix
    return(rinfo) ;  // 0 indicates success
}


// Determines the "right" eigen decomposition of a general matrix
int RightEigenSystem(array2D<double> &A , array2D<double> &Evectors ,  vector <double> &Reigenvalues,
                    vector <double> &Ieigenvalues ) {
    int N=A.numberRows();
    char rvecs = 'V' ;  // Don't get the right eigen vectors
    char lvecs = 'N' ;  // Do calculate the left eigen vectors
    int LWork=12*N;
    if (N>400)
       LWork=4*N;
    double Work[LWork];
    int rinfo , LDVL=1;

    array2D<double> A_t = transpose(A) ;
    dgeev_(&lvecs, &rvecs, &N, &(A_t(0,0)) , &N, &Reigenvalues[0], &Ieigenvalues[0], NULL, &LDVL,
      &(Evectors(0,0)), &N, Work, &LWork, &rinfo);

    // right eigen vectors are stored as the COLUMNS of the eigen matrix
    Evectors = transpose(Evectors) ;
    return(rinfo) ;  // 0 indicates success
}
