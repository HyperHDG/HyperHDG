/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "LapackWrapper.hxx"
#include "HyAssert.hxx"



extern "C" {
  void daxpy_(int* n,double* alpha,double* dx,int* incx,double* dy,int* incy);
  double dnrm2_(int* n,double* x, int* incx);

  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetrs_(char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
  void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
}


void lapack_solve(int system_size, double *mat_a, double *rhs_b)
{
  int one = 1, info = -1;
  int ipiv[system_size];
  dgesv_(&system_size, &one, mat_a, &system_size, ipiv, rhs_b, &system_size, &info);
  hy_assert( info == 0 ,
             "LAPACK's solve failed and the solution of the local problem might be inaccurate." );
}
