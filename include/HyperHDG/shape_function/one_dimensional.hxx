#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

#include <cmath>

/*!*************************************************************************************************
 * \brief   Namespace for auxiliary functions and classes needed for struct shap_function.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
 **************************************************************************************************/
namespace ShapeType
{
/*!*************************************************************************************************
 * \brief   Struct that handles the evaluation of one-dimensional Legendre polynomials.
 *
 * We use Clenshaw's algorithm (cf. https://en.wikipedia.org/wiki/Clenshaw_algorithm; date: Jan, 09,
 * 2021) to evaluate the three-term recusion formula defining Legendre polynomials as defined on
 * https://en.wikipedia.org/wiki/Legendre_polynomials; date Jan, 09, 2021).
 *
 * The shape functions in this struct, however, are defined with respect to the unit interval [0,1]
 * and normalized to be L^2 orthonormal (not just orthogonal).
 *
 * \tparam  poly_deg  The maximum degree of the polynomial.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
 * \authors   Guido Kanschat, Heidelberg University, 2021.
 **************************************************************************************************/
template <unsigned int poly_deg>
struct Legendre
{
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_fun() { return degree() + 1; }
  /*!***********************************************************************************************
   * \brief   Dimension of abscissas of the polynomials.
   ************************************************************************************************/
  static constexpr unsigned int dim() { return 1U; }
  /*!***********************************************************************************************
   * \brief   Maximum degree of all polynomials.
   ************************************************************************************************/
  static constexpr unsigned int degree() { return poly_deg; }

  /*!***********************************************************************************************
   * \brief   Evaluate variable alpha of the Clenshaw algorithm.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t clenshaw_alpha(const unsigned int k, const input_t& x)
  {
    return (((return_t)(2 * k + 1)) / ((return_t)(k + 1))) * ((return_t)x);
  }
  /*!***********************************************************************************************
   * \brief   Evaluate variable beta of the Clenshaw algorithm.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t clenshaw_beta(const unsigned int k, const input_t& UNUSED(x))
  {
    return -((return_t)k) / ((return_t)(k + 1));
  }
  /*!***********************************************************************************************
   * \brief   Evaluate variable b of the Clenshaw algorithm.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static constexpr return_t clenshaw_b(const unsigned int k, const unsigned int n, const input_t& x)
  {
    if (k > n)
      return (return_t)0.;
    else if (k == n)
      return (return_t)1.;
    else
      return clenshaw_alpha<return_t>(k, x) * clenshaw_b<return_t>(k + 1, n, x) +
             clenshaw_beta<return_t>(k + 1, x) * clenshaw_b<return_t>(k + 2, n, x);
  }

  /*!***********************************************************************************************
   * \brief   Evaluate value of one shape function.
   *
   * Evaluates value of the \c index one-dimensional shape function on the reference interval
   * \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  input_t       Floating type specification for input value.
   * \param   index         Index of evaluated shape function.
   * \param   x_val         Abscissa of evaluated shape function.
   * \param   normalized    Decide whether L^2 normalization is conducted. Defaults to true.
   * \retval  fct_value     Evaluated value of shape function.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static inline return_t fct_val(const unsigned int index,
                                 const input_t& x_val,
                                 const bool normalized = true)
  {
    // Transform x value from reference interval [0,1] -> [-1,1].
    const return_t x = 2. * (return_t)x_val - 1.;

    // Evaluate the value of the Legendre polynomial at x.
    return_t legendre_val;
    switch (index)
    {
      case 0:
        legendre_val = 1.;
        break;
      case 1:
        legendre_val = x;
        break;
      default:
        legendre_val =
          x * clenshaw_b<return_t>(1, index, x) - 0.5 * clenshaw_b<return_t>(2, index, x);
    }

    // Return L^2 normalized value.
    if (normalized)
      legendre_val *= std::sqrt((return_t)(2. * index + 1.));

    return legendre_val;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate value of the derivative of orthonormal shape function.
   *
   * Evaluates the value of the derivative of the \c index orthonormal, one-dimensional shape
   * function on the reference interval \f$[0,1]\f$ at abscissa \c x_val.
   *
   * \tparam  return_t      Floating type specification for return value.
   * \tparam  input_t       Floating type specification for input value.
   * \param   index         Index of evaluated shape function.
   * \param   x_val         Abscissa of evaluated shape function.
   * \param   normalized    Decide whether L^2 normalization is conducted. Defaults to true.
   * \param   unit_interval Derivative is evaluated with respect to unit interval. Defauts to true.
   * \retval  fct_value     Evaluated value of shape function's derivative.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static inline return_t der_val(const unsigned int index,
                                 const input_t& x_val,
                                 const bool normalized = true,
                                 const bool unit_interval = true)
  {
    if (index == 0)
      return (return_t)0.;

    // Transform x value from reference interval [0,1] -> [-1,1].
    const return_t x = 2. * (return_t)x_val - 1.;

    // Non-recusrive version of:
    // return_t legendre_val = ((return_t)index) * fct_val<return_t>(index - 1, x_val, false) +
    //                         x * der_val<return_t>(index - 1, x_val, false, false);
    return_t legendre_val = 0., x_pot = 1.;
    for (unsigned int p = index; p > 0; --p)
    {
      legendre_val += ((return_t)p) * x_pot * fct_val<return_t>(p - 1, x_val, false);
      x_pot *= x;
    }

    // Return L^2 normalized value.
    if (normalized)
      legendre_val *= std::sqrt((return_t)(2. * index + 1.));

    // Derivatives with respect to unit interval must be multiplied by 2.
    if (unit_interval)
      legendre_val *= 2.;

    return legendre_val;
  }
};  // end of struct Legendre

}  // end of namespace ShapeType
