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
 * \authors   Guido Kanschat, Heidelberg University, 2019.
 * \authors   Andreas Rupp, Heidelberg University, 2019.
 **************************************************************************************************/

#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/DenseLA.hxx>
#include <array>

// #include <exception>

/*!*************************************************************************************************
 * \brief   Exception to be thrown if LAPACK's solve fails.
 *
 * \todo    Is this the way to do it intended by Guido? If so, should we include <exception>? It
 *          works perfecly without the include. I suspect that array (or another by that point
 *          included package includes exception?!
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019.
 * \authors   Andreas Rupp, Heidelberg University, 2019.
 **************************************************************************************************/
struct LAPACKexception : public std::exception
{ const char * what () const throw ()  { return "LAPACK's function failed!"; } };

extern "C"
{
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void daxpy_(int* n, double* alpha, double* dx, int* incx, double* dy, int* incy);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  double dnrm2_(int* n, double* x, int* incx);

  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void dgetrf_(int* M, int *N, double* A, int* lda, int* IPIV, int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void dgetrs_
    (char* C, int* N, int* NRHS, double* A, int* LDA, int* IPIV, double* B, int* LDB, int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void dgesv_(int *n, int *nrhs, double *a, int *lda, int *ipiv, double *b, int *ldb, int *info);
  
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void dgeqr2_(int *m, int *n, double *a, int *lda, double *tau, double *work, int *info);   
  
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void saxpy_(int* n, float* alpha, float* dx, int* incx, float* dy, int* incy);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  double snrm2_(int* n, float* x, int* incx);

  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void sgetrf_(int* M, int *N, float* A, int* lda, int* IPIV, int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void sgetrs_
    (char* C, int* N, int* NRHS, float* A, int* LDA, int* IPIV, float* B, int* LDB, int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the 
   * use of functions \c lapack_solve that will be implemented below.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  void sgesv_(int *n, int *nrhs, float *a, int *lda, int *ipiv, float *b, int *ldb, int *info);
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
inline void lapack_solve(int system_size, int n_rhs_cols, double *mat_a, double *rhs_b)
{
  int info = -1;
  int *ipiv = new int[system_size];
  dgesv_(&system_size, &n_rhs_cols, mat_a, &system_size, ipiv, rhs_b, &system_size, &info);
  delete[] ipiv;
  if (info != 0)  throw LAPACKexception();
}

/*!*************************************************************************************************
 * \brief   QR decomposition.
 *
 * \todo
 *
 * \param  system_size  Size of the system of equations.
 * \param  mat_a        Pointer to the matrix describing the linear system of equations.
 * \param  rhs_b        Pointer to the right-hand side of the system.
 * \retval rhs_b        Pointer to the solution of the system of equations.
 **************************************************************************************************/
inline void lapack_qrDecomposition(int n_rows, int n_cols, double *mat_a, double *tau)
{
  int info = -1;
  double *work = new double[n_cols];
  dgeqr2_(&n_rows, &n_cols, mat_a, &n_rows, tau, work, &info);
  delete[] work;
  if (info != 0)  throw LAPACKexception();
}

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
inline void lapack_solve(int system_size, int n_rhs_cols, float *mat_a, float *rhs_b)
{
  int info = -1;
  int *ipiv = new int[system_size];
  sgesv_(&system_size, &n_rhs_cols, mat_a, &system_size, ipiv, rhs_b, &system_size, &info);
  delete[] ipiv;
  if (info != 0)  throw LAPACKexception();
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
template < unsigned int system_size, unsigned int n_rhs_cols = 1, typename lapack_float_t >
std::array<lapack_float_t, system_size * n_rhs_cols> lapack_solve
( 
  std::array<lapack_float_t, system_size * system_size>& dense_mat,
  std::array<lapack_float_t, system_size * n_rhs_cols>& rhs
)
{
  lapack_solve(system_size, n_rhs_cols, dense_mat.data(), rhs.data());
  return rhs;
}

/*!*************************************************************************************************
 * \brief   Solve local system of equations.
 *
 * \todo    All
 *
 * \tparam system_size  Size of the system of equations.
 * \param  mat_a        Array comprising the matrix describing the linear system of equations.
 * \param  rhs_b        Array comprising the right-hand side of the system.
 * \retval rhs_b        Array comprising the solution of the system of equations.
 **************************************************************************************************/
template < unsigned int n_rows, unsigned int n_cols, typename lapack_float_t >
std::array<lapack_float_t, n_rows * n_rows> lapack_qrQ
( std::array<lapack_float_t, n_rows * n_cols>& dense_mat )
{
  constexpr unsigned int rank = std::min(n_rows, n_cols);
  std::array<lapack_float_t, rank> tau;
  SmallMat<n_rows, n_rows, lapack_float_t> matQ = diagonal<n_rows, n_rows, lapack_float_t>(1.);
  SmallVec<n_rows, lapack_float_t> vec;

  lapack_qrDecomposition(n_rows, n_cols, dense_mat.data(), tau.data());
  
  for (unsigned int i = 0; i < rank; ++i)
  {
    for (unsigned int j = 0; j < n_rows; ++j)
      if (j < i)        vec[j] = 0.;
      else if (j == i)  vec[j] = 1.;
      else              vec[j] = dense_mat(j,i);
    matQ = matQ * ( diagonal<n_rows, n_rows, lapack_float_t>(1.) - tau[i] * dyadic_product<n_rows, n_rows, lapack_float_t>( vec, vec ) );
  }

  return matQ.data();
}