#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

#include <algorithm>  // std::max
#include <cmath>
#include <exception>
#include <functional>
#include <limits>
#include <type_traits>
#include <vector>

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
 * \brief   Why an iterative solver stopped (cf. PETSc's KSPConvergedReason).
 **************************************************************************************************/
enum class StopReason
{
  converged_rtol,           //!< ||r|| <= rtol * ||b|| (and ||r|| > atol)
  converged_atol,           //!< ||r|| <= atol (checked first, as in PETSc)
  diverged_dtol,            //!< ||r|| > dtol * ||r_0||
  diverged_max_iterations,  //!< max_iterations done without convergence
  diverged_breakdown,       //!< d^T A d is zero: A is singular or indefinite on the Krylov space
  diverged_nan              //!< a residual norm or d^T A d is not finite
};
/*!*************************************************************************************************
 * \brief   The name of a StopReason, e.g. for log output.
 **************************************************************************************************/
inline const char* to_string(const StopReason reason)
{
  switch (reason)
  {
    case StopReason::converged_rtol:
      return "converged_rtol";
    case StopReason::converged_atol:
      return "converged_atol";
    case StopReason::diverged_dtol:
      return "diverged_dtol";
    case StopReason::diverged_max_iterations:
      return "diverged_max_iterations";
    case StopReason::diverged_breakdown:
      return "diverged_breakdown";
    case StopReason::diverged_nan:
      return "diverged_nan";
  }
  return "unknown";
}
/*!*************************************************************************************************
 * \brief   Stopping criteria and monitoring of an iterative solver; the defaults are PETSc's.
 *
 * The solver stops as converged if ||r|| <= max(rtol * ||b||, atol), and as diverged if ||r|| >
 * dtol * ||r_0|| or after max_iterations iterations. ||.|| is the Euclidean norm, r the residual of
 * the solver's recursion, r_0 the initial residual b - A x_0.
 **************************************************************************************************/
struct SolverOptions
{
  double rtol = 1e-5;
  double atol = 1e-50;
  double dtol = 1e5;
  unsigned int max_iterations = 10000;
  /*!
   * Called with (iteration, ||r||) before every iteration and for the final residual, if set.
   */
  std::function<void(unsigned int, double)> monitor;
};
/*!*************************************************************************************************
 * \brief   What an iterative solver did.
 **************************************************************************************************/
struct SolverResult
{
  StopReason reason;
  unsigned int iterations;       //!< matrix-vector products in the iteration (not counting r_0)
  double residual_norm;          //!< ||r|| at the end
  double initial_residual_norm;  //!< ||r_0||, the reference of dtol
  double rhs_norm;               //!< ||b||, the reference of rtol
  bool converged() const
  {
    return reason == StopReason::converged_rtol || reason == StopReason::converged_atol;
  }
};
/*!*************************************************************************************************
 * \brief   Conjugate gradients for A x = b with A symmetric and definite (either sign).
 *
 * The operator is any callable apply(in, out) that writes out = A in; in and out are std::vector
 * of the value type of b, e.g. [&](const auto& in, auto& out) { loop.apply(in, out); } for
 * GlobalLoop::Generic. b and x are vector-like (size(), operator[]); x is the initial guess on
 * entry (r_0 = b is used without a product if x is zero) and the last iterate on return, also if
 * the solver did not converge.
 *
 * \param   apply     The operator: apply(in, out) sets out = A in.
 * \param   b         Right-hand side.
 * \param   x         Initial guess on entry, solution on return.
 * \param   options   Stopping criteria and an optional monitor.
 * \retval  result    Stop reason, iterations and residual norms.
 **************************************************************************************************/
template <class OperatorT, class RhsVectorT, class VectorT>
SolverResult cg(OperatorT&& apply,
                const RhsVectorT& b,
                VectorT& x,
                const SolverOptions& options = SolverOptions())
{
  using value_t = std::remove_cvref_t<decltype(b[0])>;
  const std::size_t n = b.size();
  hy_assert(x.size() == n, "Solution and right-hand side must be of the same size.");
  auto dot = [n](const std::vector<value_t>& u, const std::vector<value_t>& v)
  {
    double sum = 0.;
    for (std::size_t i = 0; i < n; ++i)
      sum += (double)u[i] * (double)v[i];
    return sum;
  };

  std::vector<value_t> r(n), d(n), z(n);
  for (std::size_t i = 0; i < n; ++i)
    r[i] = b[i];
  const double rhs_norm = std::sqrt(dot(r, r));
  bool zero_guess = true;
  for (std::size_t i = 0; i < n; ++i)
    zero_guess = zero_guess && x[i] == (value_t)0.;
  if (!zero_guess)
  {
    for (std::size_t i = 0; i < n; ++i)
      d[i] = x[i];
    apply(d, z);
    for (std::size_t i = 0; i < n; ++i)
      r[i] -= z[i];
  }
  for (std::size_t i = 0; i < n; ++i)
    d[i] = r[i];
  double r_square = dot(r, r);
  const double initial_residual_norm = std::sqrt(r_square);
  const double target = std::max(options.rtol * rhs_norm, options.atol);

  SolverResult result{StopReason::diverged_max_iterations, 0, initial_residual_norm,
                      initial_residual_norm, rhs_norm};
  auto stop = [&](const StopReason reason)
  {
    result.reason = reason;
    if (options.monitor)
      options.monitor(result.iterations, result.residual_norm);
    return result;
  };
  for (;; ++result.iterations)
  {
    result.residual_norm = std::sqrt(r_square);
    if (!std::isfinite(result.residual_norm))
      return stop(StopReason::diverged_nan);
    if (result.residual_norm <= target)
      return stop(result.residual_norm <= options.atol ? StopReason::converged_atol
                                                       : StopReason::converged_rtol);
    if (result.residual_norm > options.dtol * initial_residual_norm)
      return stop(StopReason::diverged_dtol);
    if (result.iterations == options.max_iterations)
      return stop(StopReason::diverged_max_iterations);
    if (options.monitor)
      options.monitor(result.iterations, result.residual_norm);

    apply(d, z);
    const double curvature = dot(d, z);
    if (!std::isfinite(curvature))
      return stop(StopReason::diverged_nan);
    if (curvature == 0.)
      return stop(StopReason::diverged_breakdown);
    const value_t alpha = r_square / curvature;
    for (std::size_t i = 0; i < n; ++i)
    {
      x[i] += alpha * d[i];
      r[i] -= alpha * z[i];
    }
    const double r_square_old = r_square;
    r_square = dot(r, r);
    const value_t beta = r_square / r_square_old;
    for (std::size_t i = 0; i < n; ++i)
      d[i] = r[i] + beta * d[i];
  }
}
/*!*************************************************************************************************
 * \brief   Execute conjugate gradient algorithm to find solution to system of equations (legacy).
 *
 * The interface of the older C++ tests and examples, on top of \c cg: the template class
 * \c ProblemT implements \c trace_to_flux, which takes a large vector and returns the matrix vector
 * product. The associated matrix is assumed to be square and symmetric positive definite.
 *
 * \tparam  ProblemT            Class to implement matrix vector multiplication.
 * \tparam  vectorT             The class name of the large vector.
 * \param   b                   Right-hand side of linear system of equations.
 * \param   problem             Class instantiation to implement matrix vector multiplication.
 * \param   n_iterations        Maximum number of iterations. 0 is default and the size of b.
 * \param   tolerance           Absolute tolerance value in 2 norm. Default is 1e-9.
 * \retval  solution            Vector sufficing Ax = b up to given tolerance if converged; throws
 *                              SolveException otherwise.
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
  SolverOptions options;
  options.rtol = 0.;
  options.atol = tolerance;
  options.dtol = std::numeric_limits<double>::infinity();
  options.max_iterations = n_iterations == 0 ? b.size() : n_iterations;
  vectorT x(b.size(), (typename vectorT::value_type)0.);
  const SolverResult result =
    cg([&](const auto& in, auto& out) { out = problem.trace_to_flux(in); }, b, x, options);
  if (!result.converged())
    throw SparseLA::SolveException();
  return x;
}

}  // end of namespace SparseLA
