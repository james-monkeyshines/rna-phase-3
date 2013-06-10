#ifndef MATRIXMATH_H
#define MATRIXMATH_H

#include <vector>

#include "array2D.h"

//using namespace TNT ;    // Template Numerical Toolkit

using namespace std;

// Finds the outer product of two vectors (arrays)
array2D<double> outer_product(vector <double> &a , vector <double> &b);

// Matrix operators
array2D<double> operator* (const array2D<double> & matrix1 , const array2D<double> & matrix2);
array2D<double> operator* (const array2D<double> & matrix , double num);
array2D<double> operator* (double num , const array2D<double> & matrix );
array2D<double> transpose(array2D<double> &mat);


array2D<double> inverse(array2D<double> &matrix);
double inverseSPD(array2D<double> &matrix, int order);
double inverseSPD2(array2D<double> &matrix, int order);
bool unit(array2D<double> & matrix) ;    

// Determines the "left" eigen decomposition of a general matrix 
int LeftEigenSystem(array2D<double> &A , array2D<double> &Evectors ,  vector <double> &Reigenvalues, vector <double> &Ieigenvalues );
// Finds the "right" eigen decomposition of a general matrix
int RightEigenSystem(array2D<double> &A , array2D<double> &Evectors ,  vector <double> &Reigenvalues, vector <double> &Ieigenvalues );


#endif //MATRIXMATH_H
