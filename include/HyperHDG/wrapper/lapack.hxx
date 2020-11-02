/*!*************************************************************************************************
 * \file    lapack.hxx
 * \brief   This file provides the function lapack_solve.
 *
 * This is a wrapper file to provide LAPACK based functions that have the ability to solve (dense)
 * local systems of linear equations in an efficient way. The functions \c daxpy_ and \c dnrm2_ are
 * needed to provide the three functions that solve local systems of equations. From these three
 * functions, one is chosen to be used in the remainder of the code (i.e., \c dgesv_ cf. LAPACK
 * manual for further details).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/

#pragma once  // Ensure that file is included only once in a single compilation.

#include <array>
// #include <exception>

namespace Wrapper
{
/*!*************************************************************************************************
 * \brief   Exception to be thrown if LAPACK's solve fails.
 *
 * \todo    Is this the way to do it intended by Guido? If so, should we include exception? It
 *          works perfecly without the include. I suspect that array (or another by that point
 *          included package includes exception?!
 **************************************************************************************************/
struct LAPACKexception : public std::exception
{
  const char* what() const throw() { return "LAPACK's function failed!"; }
};

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// MAIN FUNCTIONS THAT IMPLEMENT MATRIX OPERATIONS:
// Only the functions within this section are supposed to be used outside of this file!
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

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
 * \tparam  system_size Size of the system of equations.
 * \tparam  n_rhs_cols  Number of columns of the right hand side matrix. Defaults to 1.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the matrix describing the linear system of equations.
 * \param   rhs         Array comprising the right-hand side of the system.
 * \retval  rhs_b       Array comprising the solution of the system of equations.
 **************************************************************************************************/
template <unsigned int system_size, unsigned int n_rhs_cols = 1, typename lapack_float_t>
std::array<lapack_float_t, system_size * n_rhs_cols> lapack_solve(
  std::array<lapack_float_t, system_size * system_size>& dense_mat,
  std::array<lapack_float_t, system_size * n_rhs_cols>& rhs);
/*!*************************************************************************************************
 * \brief   Determinant of a rectangular system.
 *
 * Calculate the generalized determinant of a rectangular matrix. If the matrix is square, this is
 * the standard determinant. The determinant is determined by doing a QR decomposition based on
 * Householder transformations. Thus det(Q) = +1, if the number of rows if even and det(Q) = -1, if
 * the number of rows is odd. This number is multiplied by the diagonal entries of R.
 *
 * The matrix is destroyed using this function. Its entries might be used to generate matrices Q and
 * R of the QR descomposition.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the matrix describing the linear system of equations.
 * \retval  determinant Generalized determinant of the matrix.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
lapack_float_t lapack_det(std::array<lapack_float_t, n_rows * n_cols>& dense_mat);
/*!*************************************************************************************************
 * \brief   Matrix Q of QR decomposition.
 *
 * Return matrix Q of the Householder QR decomposition of the matrix. The matrix is destroyed using
 * this function. Its entries might be used to generate matrices Q and R of the QR descomposition.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the matrix describing the linear system of equations.
 * \retval  mat_q       Matrix Q of Householder QR decomposition.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
std::array<lapack_float_t, n_rows * n_rows> lapack_qr_decomp_q(
  std::array<lapack_float_t, n_rows * n_cols>& dense_mat);
/*!*************************************************************************************************
 * \brief   Matrix R of QR decomposition.
 *
 * Return matrix R of the Householder QR decomposition of the matrix. The matrix is destroyed using
 * this function. Its entries are the same as the ones of R of the QR descomposition.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the matrix describing the linear system of equations.
 * \retval  dense_mat   Matrix R of Householder QR decomposition.
 * \retval  mat_r       Matrix R of Householder QR decomposition.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
std::array<lapack_float_t, n_rows * n_cols>& lapack_qr_decomp_r(
  std::array<lapack_float_t, n_rows * n_cols>& dense_mat);
/*!*************************************************************************************************
 * \brief   Matrices Q and R of QR decomposition.
 *
 * Return matriices Q and R of the Householder QR decomposition of the matrix.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the matrix describing the linear system of equations.
 * \param   mat_q       Reference to empy matrix to be filled.
 * \retval  dense_mat   Matrix R of Householder QR decomposition.
 * \retval  mat_q       Matrix Q of Householder QR decomposition.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
void lapack_qr_decomp(std::array<lapack_float_t, n_rows * n_cols>& dense_mat,
                      std::array<lapack_float_t, n_rows * n_rows>& mat_q);
/*!*************************************************************************************************
 * \brief   Matrices Q and R of QR decomposition.
 *
 * Return matriices Q and R of the Householder QR decomposition of the matrix.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the matrix describing the linear system of equations.
 * \param   mat_q       Reference to empy matrix to be filled.
 * \param   mat_r       Reference to empy matrix to be filled.
 * \retval  dense_mat   Matrix R of Householder QR decomposition.
 * \retval  mat_q       Matrix Q of Householder QR decomposition.
 * \retval  mat_r       Square system of size n_cols that contains respective part of R.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
void lapack_qr_decomp(std::array<lapack_float_t, n_rows * n_cols>& dense_mat,
                      std::array<lapack_float_t, n_rows * n_rows>& mat_q,
                      std::array<lapack_float_t, n_cols * n_cols>& mat_r);

// -------------------------------------------------------------------------------------------------
/// \cond EXCLUDE_CODE
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// AUXILIARY FUNCTIONS THAT DIRECTLY USE LAPACK ROUTINES: Not to be used outside of this file!
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

extern "C"
{
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void daxpy_(int* n, double* alpha, double* dx, int* incx, double* dy, int* incy);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  double dnrm2_(int* n, double* x, int* incx);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void dgetrf_(int* M, int* N, double* A, int* lda, int* IPIV, int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void dgetrs_(char* C,
               int* N,
               int* NRHS,
               double* A,
               int* LDA,
               int* IPIV,
               double* B,
               int* LDB,
               int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void dgeqr2_(int* m, int* n, double* a, int* lda, double* tau, double* work, int* info);

  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void saxpy_(int* n, float* alpha, float* dx, int* incx, float* dy, int* incy);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  double snrm2_(int* n, float* x, int* incx);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void sgetrf_(int* M, int* N, float* A, int* lda, int* IPIV, int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void
  sgetrs_(char* C, int* N, int* NRHS, float* A, int* LDA, int* IPIV, float* B, int* LDB, int* INFO);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void sgesv_(int* n, int* nrhs, float* a, int* lda, int* ipiv, float* b, int* ldb, int* info);
  /*!***********************************************************************************************
   * \brief   This function is not (never) to be used.
   *
   * This function is \b not to be used in regular code. It only / solely is defined to allow the
   * use of functions \c lapack_solve that will be implemented below.
   ************************************************************************************************/
  void sgeqr2_(int* m, int* n, float* a, int* lda, float* tau, float* work, int* info);
}  // end of extern "C"

// -------------------------------------------------------------------------------------------------
/// \endcond
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Solve local system of equations with \c double floating point numbers --- DO NOT USE.
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
 * \param  n_rhs_cols   Number of right-hand sides the system is solved for.
 * \param  mat_a        Pointer to the matrix describing the linear system of equations.
 * \param  rhs_b        Pointer to the right-hand side of the system.
 * \retval rhs_b        Pointer to the solution of the system of equations.
 **************************************************************************************************/
inline void lapack_solve(int system_size, int n_rhs_cols, double* mat_a, double* rhs_b)
{
  int info = -1;
  int* ipiv = new int[system_size];
  dgesv_(&system_size, &n_rhs_cols, mat_a, &system_size, ipiv, rhs_b, &system_size, &info);
  delete[] ipiv;
  if (info != 0)
    throw LAPACKexception();
}
/*!*************************************************************************************************
 * \brief   Solve local system of equations with \c float floating point numbers --- DO NOT USE.
 *
 * Solve linear (dense) system of equations \f$Ax=b\f$, where \$A\$ is an \f$n \times n\f$ square
 * matrix, which enters as a pointer to \c float, \f$n\f$ is provided via \c system_size, and the
 * \c float pointer \c rhs_b is both, input (i.e., \f$b\f$) and output (i.e., \f$x\f$) of the
 * function.
 *
 * Independent of \c const expressions of functions using \c lapack_solve one should not use
 * \c mat_a after calling this function. The input that has been within \c rhs_b will have been
 * replaced by the solution of the system of equations.
 *
 * \param  system_size  Size of the system of equations.
 * \param  n_rhs_cols   Number of right-hand sides the system is solved for.
 * \param  mat_a        Pointer to the matrix describing the linear system of equations.
 * \param  rhs_b        Pointer to the right-hand side of the system.
 * \retval rhs_b        Pointer to the solution of the system of equations.
 **************************************************************************************************/
inline void lapack_solve(int system_size, int n_rhs_cols, float* mat_a, float* rhs_b)
{
  int info = -1;
  int* ipiv = new int[system_size];
  sgesv_(&system_size, &n_rhs_cols, mat_a, &system_size, ipiv, rhs_b, &system_size, &info);
  delete[] ipiv;
  if (info != 0)
    throw LAPACKexception();
}
/*!*************************************************************************************************
 * \brief   QR decomposition in \c double floating point arithmetic --- DO NOT USE.
 *
 * Perform QR decomposition of a given matrix using Householder transformations. The matrices Q and
 * R are stored according to the LAPACK function in an encoded way in matrix mat_a and vector tau.
 *
 * \param  n_rows       Number of rows of matrix.
 * \param  n_cols       Number of columns of matrix.
 * \param  mat_a        Matrix that is to be decomposed.
 * \param  tau          Space for return value, i.e. vector with auxiliary darta.
 * \retval mat_a        Encoded QR decomposition of the matrix.
 * \retval tau          Auxiliary parameters needed to reconstruct matrix Q.
 **************************************************************************************************/
inline void lapack_qr(int n_rows, int n_cols, double* mat_a, double* tau)
{
  int info = -1;
  double* work = new double[n_cols];
  dgeqr2_(&n_rows, &n_cols, mat_a, &n_rows, tau, work, &info);
  delete[] work;
  if (info != 0)
    throw LAPACKexception();
}
/*!*************************************************************************************************
 * \brief   QR decomposition in \c float floating point arithmetic --- DO NOT USE.
 *
 * Perform QR decomposition of a given matrix using Householder transformations. The matrices Q and
 * R are stored according to the LAPACK function in an encoded way in matrix mat_a and vector tau.
 *
 * \param  n_rows       Number of rows of matrix.
 * \param  n_cols       Number of columns of matrix.
 * \param  mat_a        Matrix that is to be decomposed.
 * \param  tau          Space for return value, i.e. vector with auxiliary darta.
 * \retval mat_a        Encoded QR decomposition of the matrix.
 * \retval tau          Auxiliary parameters needed to reconstruct matrix Q.
 **************************************************************************************************/
inline void lapack_qr(int n_rows, int n_cols, float* mat_a, float* tau)
{
  int info = -1;
  float* work = new float[n_cols];
  sgeqr2_(&n_rows, &n_cols, mat_a, &n_rows, tau, work, &info);
  delete[] work;
  if (info != 0)
    throw LAPACKexception();
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MAIN FUNCTIONS THAT DO NOT NEED A DENSE LINEAR ALGEBRA IMPLEMENTATION:
// Functions have been declared above and are only implemented here!
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// lapack_solve
// -------------------------------------------------------------------------------------------------

template <unsigned int system_size, unsigned int n_rhs_cols = 1, typename lapack_float_t>
std::array<lapack_float_t, system_size * n_rhs_cols> lapack_solve(
  std::array<lapack_float_t, system_size * system_size>& dense_mat,
  std::array<lapack_float_t, system_size * n_rhs_cols>& rhs)
{
  lapack_solve(system_size, n_rhs_cols, dense_mat.data(), rhs.data());
  return rhs;
}

// -------------------------------------------------------------------------------------------------
// lapack_det
// -------------------------------------------------------------------------------------------------

template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
lapack_float_t lapack_det(std::array<lapack_float_t, n_rows * n_cols>& dense_mat)
{
  constexpr unsigned int rank = std::min(n_rows, n_cols);
  std::array<lapack_float_t, rank> tau;
  lapack_qr(n_rows, n_cols, dense_mat.data(), tau.data());

  lapack_float_t determinant = 1.;
  for (unsigned int i = 0; i < rank; ++i)
    determinant *= -dense_mat[i * (n_rows + 1)];
  return determinant;
}

}  // end of namespace Wrapper

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF AUXILIARY FUNCTIONS THAT REQUIRE A DENSE LA IMPLEMENTATION:
// Do not use these functions outside of this file.
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

#include <HyperHDG/dense_la.hxx>  // Dense linear algebra that is utilized in the following.

namespace Wrapper
{
/*!*************************************************************************************************
 * \brief   Matrix Q of QR decomposition --- DO NOT USE.
 *
 * Return matrix Q of the Householder QR decomposition of the matrix. The matrix is supposed to
 * contain the Householder QR decomposition which is encoded in the style of LAPACK. The function
 * returns matrix Q in standard encoding.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the QR decomposed matrix.
 * \param   tau         Auxiliary data from LAPACK to encode matrix Q.
 * \retval  mat_q       Matrix Q of Householder QR decomposition.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, unsigned int rank, typename lapack_float_t>
inline std::array<lapack_float_t, n_rows * n_rows> get_q_from_lapack_qr_result(
  const std::array<lapack_float_t, n_rows * n_cols>& dense_mat,
  const std::array<lapack_float_t, rank>& tau)
{
  SmallMat unity = diagonal<n_rows, n_rows, lapack_float_t>(1.), matQ = unity;
  SmallVec<n_rows, lapack_float_t> vec;

  for (unsigned int i = 0; i < rank; ++i)  // i is column index column index, here!
  {
    for (unsigned int j = 0; j < n_rows; ++j)  // j is row index, here!
      if (j < i)
        vec[j] = 0.;
      else if (j == i)
        vec[j] = 1.;
      else
        vec[j] = dense_mat[i * n_rows + j];
    matQ = matQ * (unity - tau[i] * dyadic_product(vec, vec));
  }

  return matQ.data();
}
/*!*************************************************************************************************
 * \brief   Matrix Q of QR decomposition --- DO NOT USE.
 *
 * \todo    Check whether the implementation is correct, since I am not perfectly sure about the two
 *          \c std::move statements to do exactly, what I want.
 *
 * Return matrix Q of the Householder QR decomposition of the matrix. The matrix is supposed to
 * contain the Householder QR decomposition which is encoded in the style of LAPACK. The function
 * returns matrix Q in standard encoding.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the QR decomposed matrix.
 * \param   tau         Auxiliary data from LAPACK to encode matrix Q.
 * \param   mat_q       Matrix to be filled with entries of the return value.
 * \retval  mat_q       Matrix Q of Householder QR decomposition.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, unsigned int rank, typename lapack_float_t>
inline void get_q_from_lapack_qr_result(
  const std::array<lapack_float_t, n_rows * n_cols>& dense_mat,
  const std::array<lapack_float_t, rank>& tau,
  std::array<lapack_float_t, n_rows * n_rows>& mat_q)
{
  SmallMat unity = diagonal<n_rows, n_rows, lapack_float_t>(1.);
  SmallMat<n_rows, n_rows, lapack_float_t> matQ(std::move(mat_q));
  SmallVec<n_rows, lapack_float_t> vec;
  matQ = unity;

  for (unsigned int i = 0; i < rank; ++i)  // i is column index column index, here!
  {
    for (unsigned int j = 0; j < n_rows; ++j)  // j is row index, here!
      if (j < i)
        vec[j] = 0.;
      else if (j == i)
        vec[j] = 1.;
      else
        vec[j] = dense_mat[i * n_rows + j];
    matQ = matQ * (unity - tau[i] * dyadic_product(vec, vec));
  }

  mat_q = std::move(matQ.data());
}
/*!*************************************************************************************************
 * \brief   Matrix R of QR decomposition --- DO NOT USE.
 *
 * Return matrix R of the Householder QR decomposition of the matrix. The matrix is supposed to
 * contain the Householder QR decomposition which is encoded in the style of LAPACK. The function
 * returns matrix R in standard encoding.
 *
 * The function also manipulates the input matrix to be R.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   dense_mat   Array comprising the QR decomposed matrix.
 * \retval  mat_r       Matrix R of Householder QR decomposition.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
inline std::array<lapack_float_t, n_rows * n_cols>& get_r_from_lapack_qr_result(
  std::array<lapack_float_t, n_rows * n_cols>& dense_mat)
{
  for (unsigned int i = 0; i < n_cols; ++i)    // i is column index, here!
    for (unsigned int j = 0; j < n_rows; ++j)  // j is row index, here!
      if (j > i)
        dense_mat[i * n_rows + j] = 0.;
  return dense_mat;
}
/*!*************************************************************************************************
 * \brief   Matrix R of QR decomposition --- DO NOT USE.
 *
 * Return matrix R of the Householder QR decomposition of the matrix. The matrix is supposed to
 * contain the Householder QR decomposition which is encoded in the style of LAPACK. The function
 * returns matrix R in standard encoding, but with dimensions \c n_cols times \c n_cols.
 *
 * \tparam  n_rows      Number of rows of the matrix whose determinant should be calculated.
 * \tparam  n_cols      Number of columns of the matrix whose determinant should be calculated.
 * \tparam  float_t     Floating type which this function should be executed with. Only \c float and
 *                      \c double are supported.
 * \param   lapack_mat  Array comprising the QR decomposed matrix.
 * \param   mat_r       Matrix to be filled with the entries of the return value.
 * \retval  mat_r       Matrix R of Householder QR decomposition.
 **************************************************************************************************/
template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
inline void get_r_from_lapack_qr_result(
  const std::array<lapack_float_t, n_rows * n_cols>& lapack_mat,
  std::array<lapack_float_t, n_cols * n_cols>& mat_r)
{
  static_assert(n_rows >= n_cols, "Function only defined for these matrices!");
  for (unsigned int i = 0; i < n_cols; ++i)    // i is column index, here!
    for (unsigned int j = 0; j < n_cols; ++j)  // j is row index, here!
      if (j > i)
        mat_r[i * n_rows + j] = 0.;
      else
        mat_r[i * n_cols + j] = lapack_mat[i * n_rows + j];
}

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MAIN FUNCTIONS THAT NEED A DENSE LINEAR ALGEBRA IMPLEMENTATION:
// Functions have been declared above and are only implemented here!
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

// -------------------------------------------------------------------------------------------------
// lapack_qr_decomp_q
// -------------------------------------------------------------------------------------------------

template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
std::array<lapack_float_t, n_rows * n_rows> lapack_qr_decomp_q(
  std::array<lapack_float_t, n_rows * n_cols>& dense_mat)
{
  constexpr unsigned int rank = std::min(n_rows, n_cols);
  std::array<lapack_float_t, rank> tau;
  lapack_qr(n_rows, n_cols, dense_mat.data(), tau.data());
  return get_q_from_lapack_qr_result<n_rows, n_cols, rank, lapack_float_t>(dense_mat, tau);
}

// -------------------------------------------------------------------------------------------------
// lapack_qr_decomp_r
// -------------------------------------------------------------------------------------------------

template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
std::array<lapack_float_t, n_rows * n_cols>& lapack_qr_decomp_r(
  std::array<lapack_float_t, n_rows * n_cols>& dense_mat)
{
  constexpr unsigned int rank = std::min(n_rows, n_cols);
  std::array<lapack_float_t, rank> tau;
  lapack_qr(n_rows, n_cols, dense_mat.data(), tau.data());
  return get_r_from_lapack_qr_result<n_rows, n_cols, rank, lapack_float_t>(dense_mat);
}

// -------------------------------------------------------------------------------------------------
// lapack_qr_decomp
// -------------------------------------------------------------------------------------------------

template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
void lapack_qr_decomp(std::array<lapack_float_t, n_rows * n_cols>& dense_mat,
                      std::array<lapack_float_t, n_rows * n_rows>& mat_q)
{
  constexpr unsigned int rank = std::min(n_rows, n_cols);
  std::array<lapack_float_t, rank> tau;
  lapack_qr(n_rows, n_cols, dense_mat.data(), tau.data());

  get_q_from_lapack_qr_result<n_rows, n_cols, rank, lapack_float_t>(dense_mat, tau, mat_q);
  get_r_from_lapack_qr_result<n_rows, n_cols, lapack_float_t>(dense_mat);
}

// -------------------------------------------------------------------------------------------------
// lapack_qr_decomp
// -------------------------------------------------------------------------------------------------

template <unsigned int n_rows, unsigned int n_cols, typename lapack_float_t>
void lapack_qr_decomp(std::array<lapack_float_t, n_rows * n_cols>& dense_mat,
                      std::array<lapack_float_t, n_rows * n_rows>& mat_q,
                      std::array<lapack_float_t, n_cols * n_cols>& mat_r)
{
  constexpr unsigned int rank = std::min(n_rows, n_cols);
  std::array<lapack_float_t, rank> tau;
  lapack_qr(n_rows, n_cols, dense_mat.data(), tau.data());

  get_q_from_lapack_qr_result<n_rows, n_cols, rank, lapack_float_t>(dense_mat, tau, mat_q);
  get_r_from_lapack_qr_result<n_rows, n_cols, lapack_float_t>(dense_mat, mat_r);
  get_r_from_lapack_qr_result<n_rows, n_cols, lapack_float_t>(dense_mat);
}

}  // end of namespace Wrapper
