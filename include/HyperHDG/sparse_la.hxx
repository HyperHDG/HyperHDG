#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

// #include <cmath>
// #include <exception>

/*!*************************************************************************************************
 * \brief   A namespace containing different functions that implement basic linear algebra
 *          operations using large vectors.
 *
 * This namespace provides several functions to implement basic linear algebra operations of (long)
 * vector type in combination with a class providing a function \c matrix_vector_multiply. This
 * is mainly used for C++ examples and test cases that do not use the Python interface and its
 * version of an CG method, for example.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
namespace SparseLA
{
/*!*************************************************************************************************
 * \brief   Exception to be thrown if conjugate gradient fails.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019.
 * \authors   Andreas Rupp, Heidelberg University, 2019.
 **************************************************************************************************/
struct SolveException : public std::exception
{
  const char* what() const throw() { return "The sparse CG method did not converge."; }
};
/*!*************************************************************************************************
 * \brief   Evaluate the inner product of two vectors.
 *
 * Naive implementation of an Euclidean inner product of two large vectors which are supposed to be
 * of the same size. This function is needed to calculate a vector's 2 norm or to implement a CG
 * scheme.
 *
 * \tparam  vectorT             The class name of the large vector.
 * \param   left                Left argument of the inner product.
 * \param   right               Right argument of the inner product.
 * \retval  product             Inner product of the two arguments.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename vectorT>
typename vectorT::value_type inner_product(const vectorT& left, const vectorT& right)
{
  hy_assert(left.size() == right.size(), "Both vectors of inner product must be of same size!");

  typename vectorT::value_type product = 0.;
  for (decltype(left.size()) i = 0; i < left.size(); ++i)
    product += left[i] * right[i];
  return product;
}
/*!*************************************************************************************************
 * \brief   Evaluate 2 norm of a vector.
 *
 * Naive implementation of an 2 norm of a vector. This is the square root of the \c inner_product of
 * a vector paired with itself.
 *
 * \tparam  vectorT             The class name of the large vector.
 * \param   vec                 Vector whose 2 norm is to be calculates.
 * \retval  norm                2 norm of given vector.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename vectorT>
typename vectorT::value_type norm_2(const vectorT& vec)
{
  return std::sqrt(inner_product(vec, vec));
}
/*!*************************************************************************************************
 * \brief   Evaluate linear combination of vectors and return the result.
 *
 * This functions takes two large vectors and two floating points, and returns their linear
 * combination "leftFac * leftVec + rightFac * rightVec" as a new vector (in contrast to just a
 * reference to a vector).
 *
 * \tparam  vectorT             The class name of the large vector.
 * \param   leftFac             Scaling factor of left vector.
 * \param   leftVec             Left vector in linear combination.
 * \param   rightFac            Scaling factor of right vector.
 * \param   rightVec            Right vector in linear combination.
 * \retval  lin_comb            Linear combination of vectors with respective coefficients.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename vectorT>
vectorT linear_combination(const typename vectorT::value_type leftFac,
                           const vectorT& leftVec,
                           const typename vectorT::value_type rightFac,
                           const vectorT& rightVec)
{
  hy_assert(leftVec.size() == rightVec.size(), "Linear combined vectors must be of same size!");

  vectorT lin_comb(leftVec.size(), 0.);
  for (decltype(leftVec.size()) i = 0; i < leftVec.size(); ++i)
    lin_comb[i] = leftFac * leftVec[i] + rightFac * rightVec[i];
  return lin_comb;
}
/*!*************************************************************************************************
 * \brief   Evaluate linear combination of vectors and return reference to result.
 *
 * This functions takes two large vectors and two floating points, and returns their linear
 * combination "leftFac * leftVec + rightFac * rightVec" as a reference to a vector. This vector
 * needs to be passed to the function
 *
 * \tparam  vectorT             The class name of the large vector.
 * \param   leftFac             Scaling factor of left vector.
 * \param   leftV               Left vector in linear combination.
 * \param   rightFac            Scaling factor of right vector.
 * \param   rightV              Right vector in linear combination.
 * \param   result              Reference to vector whicb is supposed to contain the result.
 * \retval  result              Linear combination of vectors with respective coefficients.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename vectorT>
void linear_combination(const typename vectorT::value_type leftFac,
                        const vectorT& leftV,
                        const typename vectorT::value_type rightFac,
                        const vectorT& rightV,
                        vectorT& result)
{
  hy_assert(leftV.size() == rightV.size() && leftV.size() == result.size(),
            "All three vectors of linear combination must be of same size!");

  for (decltype(leftV.size()) i = 0; i < result.size(); ++i)
    result[i] = leftFac * leftV[i] + rightFac * rightV[i];
}
/*!*************************************************************************************************
 * \brief   Execute conjugate gradient algorithm to find solution to system of equations.
 *
 * Execute conjugate gradient algorithm where the matrix is not explicitly given, but the template
 * class \c ProblemT is supposed to implement a function \c matrix_vector_multiply which only takes
 * a \c large vector and generates the matrix vector product from that. The associated matrix is
 * assumed to be square and symmetric positive definite.
 *
 * \tparam  ProblemT            Class to implement matrix vector multiplication.
 * \tparam  vectorT             The class name of the large vector.
 * \param   b                   Right-hand side of linear system of equations.
 * \param   problem             Class instantiation to implement matrix vector multiplication.
 * \param   n_iterations        Maximum number of iterations. 0 is default and the size of b.
 * \param   tolerance           Absolute tolerance value in 2 norm. Default is 1e-9.
 * \retval  solution            Vector sufficing Ax = b up to given tolerance if converged.
 * \retval  n_iterations        Number of needed iterations. -1 indicates no convergence.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <class ProblemT, typename vectorT>
vectorT conjugate_gradient(const vectorT& b,
                           ProblemT& problem,
                           unsigned int n_iterations = 0,
                           const typename vectorT::value_type tolerance = 1e-9)
{
  vectorT x(b.size(), (typename vectorT::value_type)0.);
  vectorT r = b;  // b - A * x (with x = 0)
  vectorT d = r;

  typename vectorT::value_type r_square_new = inner_product(r, r);
  typename vectorT::value_type r_square_old = r_square_new;

  if (r_square_new < tolerance * tolerance)
  {
    n_iterations = 0;
    return x;
  }

  if (n_iterations == 0)
    n_iterations = b.size();
  hy_assert(n_iterations > 0, "Number of allowed iterations of CG solver must be positive.");

  for (unsigned int k = 0; k < n_iterations; ++k)
  {
    vectorT z = problem.trace_to_flux(d);
    r_square_old = r_square_new;

    typename vectorT::value_type alpha = r_square_old / inner_product(d, z);
    linear_combination((typename vectorT::value_type)1., x, alpha, d, x);
    linear_combination((typename vectorT::value_type)1., r, -alpha, z, r);

    r_square_new = inner_product(r, r);

    typename vectorT::value_type beta = r_square_new / r_square_old;
    linear_combination((typename vectorT::value_type)1., r, beta, d, d);

    if (r_square_new < tolerance * tolerance)
    {
      n_iterations = k + 1;
      return x;
    }
  }

  throw SparseLA::SolveException();

  return x;
}

template <class ProblemT, typename vectorT, typename paramT>
vectorT conjugate_gradient(const vectorT& b,
                           ProblemT& problem,
                           paramT time,
                           unsigned int n_iterations = 0,
                           const typename vectorT::value_type tolerance = 1e-9)
{
  vectorT x(b.size(), (typename vectorT::value_type)0.);
  vectorT r = b;  // b - A * x (with x = 0)
  vectorT d = r;

  typename vectorT::value_type r_square_new = inner_product(r, r);
  typename vectorT::value_type r_square_old = r_square_new;

  if (r_square_new < tolerance * tolerance)
  {
    n_iterations = 0;
    return x;
  }

  if (n_iterations == 0)
    n_iterations = b.size();
  hy_assert(n_iterations > 0, "Number of allowed iterations of CG solver must be positive.");

  for (unsigned int k = 0; k < n_iterations; ++k)
  {
    vectorT z = problem.trace_to_flux(d, time);
    r_square_old = r_square_new;

    typename vectorT::value_type alpha = r_square_old / inner_product(d, z);
    linear_combination((typename vectorT::value_type)1., x, alpha, d, x);
    linear_combination((typename vectorT::value_type)1., r, -alpha, z, r);

    r_square_new = inner_product(r, r);

    typename vectorT::value_type beta = r_square_new / r_square_old;
    linear_combination((typename vectorT::value_type)1., r, beta, d, d);

    if (r_square_new < tolerance * tolerance)
    {
      n_iterations = k + 1;
      return x;
    }
  }

  throw SparseLA::SolveException();

  return x;
}

}  // end of namespace SparseLA
