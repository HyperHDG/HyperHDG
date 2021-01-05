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
 * \tparam  poly_deg  The maximum degree of the polynomial.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2021.
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
   * \brief   Evaluate value of one shape function.
   *
   * Evaluates value of the \c index one-dimensional shape function on the reference interval
   * \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t    Floating type specification for return value.
   * \tparam  input_t     Floating type specification for input value.
   * \param   index       Index of evaluated shape function.
   * \param   x_val       Abscissa of evaluated shape function.
   * \retval  fct_value   Evaluated value of shape function.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static inline return_t fct_val(const unsigned int index, const input_t x_val)
  {
    hy_assert(index <= 5 && index < poly_deg + 1,
              "The index of a shape function must be non-negative and smaller than or equal to 5 "
                << "at the moment. Your choice has been " << index << ".");
    hy_assert(0. <= x_val && x_val <= 1.,
              "The abscissa / x for which the shape function is evaluated has been set to be in "
                << "the closed interval [0,1]. Your choice has been " << x_val << ".");

    const return_t x = (return_t)x_val;

    switch (index)
    {
      case 0:
        return 1.;
      case 1:
        return std::sqrt(3) * (1. - 2. * x);
      case 2:
        return std::sqrt(5) * ((6. * x - 6.) * x + 1.);
      case 3:
        return std::sqrt(7) * (((20. * x - 30.) * x + 12.) * x - 1.);
      case 4:
        return std::sqrt(9) * ((((70. * x - 140.) * x + 90.) * x - 20.) * x + 1.);
      case 5:
        return std::sqrt(11) * (((((252. * x - 630.) * x + 560.) * x - 210.) * x + 30.) * x - 1.);
      default:
        hy_assert(0 == 1, "This shape function has not yet been implemented.");
    }

    hy_assert(0 == 1, "Something went wrong when evaluating a shape function. This message, "
                        << "however, is never supposed to appear. Something went wrong in the core "
                        << "program.");
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate value of the derivative of orthonormal shape function.
   *
   * Evaluates the value of the derivative of the \c index orthonormal, one-dimensional shape
   * function on the reference interval \f$[0,1]\f$ at abscissa \c x_val.
   *
   * \tparam  return_t    Floating type specification for return value.
   * \tparam  input_t     Floating type specification for input value.
   * \param   index       Index of evaluated shape function.
   * \param   x_val       Abscissa of evaluated shape function.
   * \retval  fct_value   Evaluated value of shape function's derivative.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static inline return_t der_val(const unsigned int index, const input_t x_val)
  {
    hy_assert(index <= 5 && index < poly_deg + 1,
              "The index of a shape function must be non-negative and smaller than or equal to 5 "
                << "at the moment. Your choice has been " << index << ".");
    hy_assert(0. <= x_val && x_val <= 1.,
              "The abscissa / x for which the shape function is evaluated has been set to be in "
                << "the closed interval [0,1]. Your choice has been " << x_val << ".");

    const return_t x = (return_t)x_val;

    switch (index)
    {
      case 0:
        return 0.;
      case 1:
        return -std::sqrt(12);
      case 2:
        return std::sqrt(5) * (12. * x - 6.);
      case 3:
        return std::sqrt(7) * ((60. * x - 60.) * x + 12.);
      case 4:
        return std::sqrt(9) * (((280. * x - 420.) * x + 180.) * x - 20.);
      case 5:
        return std::sqrt(11) * ((((1260. * x - 2520.) * x + 1680.) * x - 420.) * x + 30.);
      default:
        hy_assert(0 == 1, "This shape function has not yet been implemented.");
    }

    hy_assert(0 == 1, "Something went wrong when evaluating a shape function. This message, "
                        << "however, is never supposed to appear. Something went wrong in the core "
                        << "program.");
    return 0.;
  }
};  // end of struct Legendre

}  // end of namespace ShapeType
