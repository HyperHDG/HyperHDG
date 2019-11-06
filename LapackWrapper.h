/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 * 
 * Providing lapack_solve(..) function to solve mat_a^T x = rhs_b, where x will be written to rhs_b and
 * info provides information whether everything worked out (info == 0).
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef LAPACKWRAPPER_H
#define LAPACKWRAPPER_H

extern "C" {
  void daxpy_(int* n,double* alpha,double* dx,int* incx,double* dy,int* incy);
  double dnrm2_(int* n,double* x, int* incx);

  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetrs_(char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
  void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

  void lapack_solve(int system_size, double *mat_a, double *rhs_b, int *info)
  {
    int one = 1;
    int ipiv[system_size];
    dgesv_(&system_size, &one, mat_a, &system_size, ipiv, rhs_b, &system_size, info);
  }
}

#endif
