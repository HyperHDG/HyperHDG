/*!*************************************************************************************************
 * \file    LapackWrapper.hxx
 * \brief   This file provides the function lapack_solve.
 *
 * This is a wrapper file to provide LAPACK based functions that have the ability to solve (dense)
 * local systems of linear equations in an efficient way. The functions \c daxpy_ and \c dnrm2_ are
 * needed to provide the three functions that solve local systems of equations. From these three
 * functions, one is chosen to be used in the remainder of the code (i.e., \c dgesv_ cf. LAPACK
 * manual for further details).
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019.
 * \authors   Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************/

#ifndef LAPACKWRAPPER_HXX
#define LAPACKWRAPPER_HXX

#include <array>

extern "C"
{
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, University of Heidelberg, 2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2020.
   ************************************************************************************************/
  void daxpy_(int* n,double* alpha,double* dx,int* incx,double* dy,int* incy);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, University of Heidelberg, 2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2020.
   ************************************************************************************************/
  double dnrm2_(int* n,double* x, int* incx);

  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, University of Heidelberg, 2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2020.
   ************************************************************************************************/
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, University of Heidelberg, 2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2020.
   ************************************************************************************************/
  void dgetrs_
    (char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
    /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, University of Heidelberg, 2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2020.
   ************************************************************************************************/
  void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
} // end of extern "C"

/*!*************************************************************************************************
 * \brief   Solve local system of equations.
 *
 * Solve linear (dense) system of equations \f$Ax=b\f$, where \$A\$ is an \f$n \times n\f$ square
 * matrix, which enters as a pointer to \c double, \f$n\f$ is provided via \c system_size, and the
 * \c double pointer \c rhs_b is both, input (i.e., \f$b\f$) and output (i.e., \f$x\f$) of the
 * function.
 *
 * Independent of \c const expressions of functions using \c lapack_solve one should not use
 * \c mat_a after calling this function. The input that has been within \c rhs_b will have been
 * replaced by the solution of the system of equations.
 *
 * \param  system_size  Size of the system of equations.
 * \param  mat_a        Pointer to the matrix describing the linear system of equations.
 * \param  rhs_b        Pointer to the right-hand side of the system.
 * \retval rhs_b        Pointer to the solution of the system of equations.
 **************************************************************************************************/
inline void lapack_solve(int system_size, double *mat_a, double *rhs_b)
{
  int one = 1, info = -1;
  int ipiv[system_size];
  dgesv_(&system_size, &one, mat_a, &system_size, ipiv, rhs_b, &system_size, &info);
  hy_assert( info == 0 ,
             "LAPACK's solve failed and the solution of the local problem might be inaccurate." );
}

/*!*************************************************************************************************
 * \brief   Solve local system of equations.
 *
 * Solve linear (dense) system of equations \f$Ax=b\f$, where \$A\$ is an \f$n \times n\f$ square
 * matrix, which enters as a \c std::array of \c double, \f$n\f$ is provided via \c system_size, and
 * the \c std::array of \c double \c rhs_b is both, input (i.e., \f$b\f$) and output (i.e., \f$x\f$)
 * of the function.
 *
 * Independent of \c const expressions of functions using \c lapack_solve one should not use
 * \c mat_a after calling this function. The input that has been within \c rhs_b will have been
 * replaced by the solution of the system of equations.
 *
 * \tparam system_size  Size of the system of equations.
 * \param  mat_a        Array comprising the matrix describing the linear system of equations.
 * \param  rhs_b        Array comprising the right-hand side of the system.
 * \retval rhs_b        Array comprising the solution of the system of equations.
 **************************************************************************************************/
template<unsigned int system_size> std::array<double, system_size> lapack_solve
(std::array<double, system_size * system_size>& dense_mat, std::array<double, system_size>& rhs)
{
  double *mat_a = dense_mat.data();
  double *rhs_b = rhs.data();
  lapack_solve(system_size, mat_a, rhs_b);
  return rhs;
}

#endif // end of ifndef LAPACKWRAPPER_HXX