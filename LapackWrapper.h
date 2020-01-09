/*!*************************************************************************************************
 * @file    LapackWrapper.h
 * @brief   This file provides the function lapack_solve.
 *
 * This is a wrapper file to provide LAPACK based functions that have the ability to solve (dense)
 * local systems of linear equations in an efficient way. The functions @c daxpy_ and @c dnrm2_ are
 * needed to provide the three functions that solve local systems of equations. From these three
 * functions, one is chosen to be used in the remainder of the code (i.e., @c dgesv_ cf. LAPACK
 * manual for further details).
 *
 * @authors   Guido Kanschat, University of Heidelberg, 2019.
 * @authors   Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************/

#ifndef LAPACKWRAPPER_H
#define LAPACKWRAPPER_H

extern "C" {
  void daxpy_(int* n,double* alpha,double* dx,int* incx,double* dy,int* incy);
  double dnrm2_(int* n,double* x, int* incx);

  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  void dgetrs_(char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
  void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);

  /*!***********************************************************************************************
   * @brief   Solve local system of equations.
   *
   * Solve linear (dense) system of equations @f$Ax=b@f$, where @$A@$ is an @f$n \times n@f$ square
   * matrix, which enters as a pointer to @c double, @f$n@f$ is provided via @c system_size, and the
   * @c double pointer @c rhs_b is both, input (i.e., @f$b@f$) and output (i.e., @f$x@f$) of the
   * function.
   *
   * Independent of @c const expressions of functions using @c lapack_solve one should not use
   * @c mat_a after calling this function. The input that has been within @c rhs_b will have been
   * replaced by the solution of the system of equations if @c info comprises 0. Otherwise, the
   * local solve is likely to have failed and @c rhs_b contains no valuable information.
   *
   * @param  system_size  Size of the system of equations.
   * @param  mat_a        Pointer to the matrix describing the linear system of equations.
   * @param  rhs_b        Pointer to the right-hand side of the system.
   * @retval rhs_b        If @c info is 0, then this is a pointer to the solution of the system of
   *                      equations.
   * @retval info         Describes whether the algorithm has been able to compute the solution.
   *                      0 = success.
   ************************************************************************************************/
  void lapack_solve(int system_size, double *mat_a, double *rhs_b, int *info)
  {
    int one = 1;
    int ipiv[system_size];
    dgesv_(&system_size, &one, mat_a, &system_size, ipiv, rhs_b, &system_size, info);
  }
}

#endif
