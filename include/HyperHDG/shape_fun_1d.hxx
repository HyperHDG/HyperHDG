#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

#include <array>
#include <cmath>

/*!*************************************************************************************************
 * \todo  Introduce enum of structs that can be plugged into ShapeFun1D's template parameters.
 **************************************************************************************************/


// Definition of Legendre type shape functions

struct Legendre
{
  /*!***********************************************************************************************
   * \brief   Evaluate value of orthonormal shape function.
   *
   * Evaluates the value of the \c index orthonormal, one-dimensional shape function on the
   * reference interval \f$[0,1]\f$ at abscissa \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  input_t             Floating type specification for input value.
   * \param   index               Index of evaluated shape function.
   * \param   x_val               Abscissa of evaluated shape function.
   * \retval  fct_value           Evaluated value of shape function.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static inline return_t shape_fct_eval(const unsigned int index, const input_t x_val)
  {
    hy_assert(index <= 5,
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
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  input_t             Floating type specification for input value.
   * \param   index               Index of evaluated shape function.
   * \param   x_val               Abscissa of evaluated shape function.
   * \retval  fct_value           Evaluated value of shape function's derivative.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static inline return_t shape_der_eval(const unsigned int index, const input_t x_val)
  {
    hy_assert(index <= 5,
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

// Shape functions & their derivatives (evaluation):

template <typename shape_t>
struct ShapeFun1D
{
  /*!***********************************************************************************************
   * \brief   Evaluate value of one shape function.
   *
   * Evaluates value of the \c index one-dimensional shape function on the reference interval
   * \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  input_t             Floating type specification for input value.
   * \param   index               Index of evaluated shape function.
   * \param   x_val               Abscissa of evaluated shape function.
   * \retval  fct_value           Evaluated value of shape function.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static inline return_t shape_fct_eval(const unsigned int index, const input_t x_val)
  {
    return shape_t::template shape_fct_eval<return_t>(index, x_val);
  }
  /*!***********************************************************************************************
   * \brief   Evaluate values of the derivative of one shape function.
   *
   * Evaluates value of the derivative of the \c index one-dimensional shape function on the
   * reference interval \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  input_t             Floating type specification for input value.
   * \param   index               Index of evaluated shape function.
   * \param   x_val               Abscissa of evaluated shape function.
   * \retval  fct_value           Evaluated value of shape function's derivative.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t, typename input_t>
  static inline return_t shape_der_eval(const unsigned int index, const input_t x_val)
  {
    return shape_t::template shape_der_eval<return_t>(index, x_val);
  }

  /*!***********************************************************************************************
   * \brief   Evaluate several values of one shape function.
   *
   * Evaluates several values of the \c index one-dimensional shape function on the reference
   * interval \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  input_t             Floating type specification for input value.
   * \tparam  sizeX               Size of array of x values.
   * \param   index               Index of evaluated shape function.
   * \param   x_val               Abscissas of evaluated shape function.
   * \retval  fct_value           Evaluated value of shape function.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t, typename input_t, std::size_t sizeX>
  static inline std::array<return_t, sizeX> shape_fct_eval(const unsigned int index,
                                                    const std::array<input_t, sizeX>& x_val)
  {
    std::array<return_t, sizeX> result;
    for (unsigned int k = 0; k < sizeX; ++k)
      result[k] = shape_fct_eval<return_t, shape_t>(index, x_val[k]);
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate several values of the derivative of one shape function.
   *
   * Evaluates several values of the derivative of the \c index one-dimensional shape function in
   * the reference interval \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  input_t             Floating type specification for input value.
   * \tparam  sizeX               Size of array of x values.
   * \param   index               Index of evaluated shape function.
   * \param   x_val               Abscissas of evaluated shape function.
   * \retval  fct_value           Evaluated value of shape function's derivative.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t, typename input_t, std::size_t sizeX>
  static inline std::array<return_t, sizeX> shape_der_eval(const unsigned int index,
                                                    const std::array<input_t, sizeX>& x_val)
  {
    std::array<return_t, sizeX> result;
    for (unsigned int k = 0; k < sizeX; ++k)
      result[k] = shape_der_eval<return_t, shape_t>(index, x_val[k]);
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate one value of several shape functions.
   *
   * Evaluate the value of the \c index one-dimensional shape functions on the reference interval
   * \f$[0,1]\f$ at abscissa \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  input_t             Floating type specification for input value.
   * \tparam  sizeInd             Size of array of inidces of polynomial degrees.
   * \param   index               Indices of evaluated shape functions.
   * \param   x_val               Abscissa of evaluated shape function.
   * \retval  fct_value           Evaluated value of shape function.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t, typename input_t, std::size_t sizeInd>
  static inline std::array<return_t, sizeInd> shape_fct_eval(
    const std::array<unsigned int, sizeInd>& index,
    const input_t x_val)
  {
    std::array<return_t, sizeInd> result;
    for (unsigned int k = 0; k < sizeInd; ++k)
      result[k] = shape_fct_eval<return_t, shape_t>(index[k], x_val);
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate one value of several shape functions' derivatives.
   *
   * Evaluates the value of the derivatives of \c index one-dimensional shape functions on the
   * reference interval \f$[0,1]\f$ at abscissa \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  input_t             Floating type specification for input value.
   * \tparam  sizeInd             Size of array of inidces of polynomial degrees.
   * \param   index               Indices of evaluated shape functions.
   * \param   x_val               Abscissa of evaluated shape function.
   * \retval  fct_value           Evaluated value of shape functions' derivatives.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t, typename input_t, std::size_t sizeInd>
  static inline std::array<return_t, sizeInd> shape_der_eval(
    const std::array<unsigned int, sizeInd>& index,
    const input_t x_val)
  {
    std::array<return_t, sizeInd> result;
    for (unsigned int k = 0; k < sizeInd; ++k)
      result[k] = shape_der_eval<return_t, shape_t>(index[k], x_val);
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate several values of several shape functions.
   *
   * Evaluates the values of the \c index one-dimensional shape functions on the reference interval
   * \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  input_t             Floating type specification for input value.
   * \tparam  sizeInd             Size of array of inidces of polynomial degrees.
   * \tparam  sizeX               Size of array of x values.
   * \param   index               Indices of evaluated shape function.
   * \param   x_val               Abscissas of evaluated shape function.
   * \retval  fct_value           Evaluated values of shape functions.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t,
            typename input_t,
            std::size_t sizeInd,
            std::size_t sizeX>
  static inline std::array<std::array<return_t, sizeX>, sizeInd> shape_fct_eval(
    const std::array<unsigned int, sizeInd>& index,
    const std::array<input_t, sizeX>& x_val)
  {
    std::array<std::array<return_t, sizeX>, sizeInd> result;
    for (unsigned int k = 0; k < sizeInd; ++k)
      result[k] = shape_fct_eval<return_t, shape_t>(index[k], x_val);
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate several values of several shape functions' derivatives.
   *
   * Evaluates the values of the \c index one-dimensional shape functions' derivatives on
   * the reference interval \f$[0,1]\f$ at abscissas \c x_val.
   *
   * \tparam  return_t            Floating type specification for return value.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  input_t             Floating type specification for input value.
   * \tparam  sizeInd             Size of array of inidces of polynomial degrees.
   * \tparam  sizeX               Size of array of x values.
   * \param   index               Indices of evaluated shape function.
   * \param   x_val               Abscissas of evaluated shape function.
   * \retval  fct_value           Evaluated values of shape functions' derivatives.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t,
            typename input_t,
            std::size_t sizeInd,
            std::size_t sizeX>
  static inline std::array<std::array<return_t, sizeX>, sizeInd> shape_der_eval(
    const std::array<unsigned int, sizeInd>& index,
    const std::array<input_t, sizeX>& x_val)
  {
    std::array<std::array<return_t, sizeX>, sizeInd> result;
    for (unsigned int k = 0; k < sizeInd; ++k)
      result[k] = shape_der_eval<return_t, shape_t>(index[k], x_val);
    return result;
  }

  /*!***********************************************************************************************
   * \brief   Shape functions evaluated at end points of unit interval.
   *
   * Returns the values of the shape functions on \f$[0,1]\f$ of degree at most
   * \c max_poly_degree at the value \f$0\f$ (at index 0) and at \f$1\f$ (at index 1).
   *
   * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  return_t            Floating type specification. Default is double.
   * \retval  corner_vals         \c std::array of polynomial degrees containing \c std::array of
   *                              corner indices (the shape functions are evaluated at).
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <unsigned int max_poly_degree, typename return_t = double>
  static inline std::array<std::array<return_t, 2>, max_poly_degree + 1> shape_fcts_at_bdrs()
  {
    std::array<return_t, 2> bdrs = {0., 1.};
    std::array<unsigned int, max_poly_degree + 1> poly_deg_index;
    for (unsigned int i = 0; i < poly_deg_index.size(); ++i)
      poly_deg_index[i] = i;

    return shape_fct_eval<return_t, shape_t>(poly_deg_index, bdrs);
  }
  /*!***********************************************************************************************
   * \brief   Derivatives of shape functions evaluated at end points of unit interval.
   *
   * Returns the values of the shape functions' derivatives on \f$[0,1]\f$ of degree at most
   * \c max_poly_degree at the value \f$0\f$ (at index 0) and at \f$1\f$ (at index 1).
   *
   * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
   * \tparam  shape_t             Type of shape functions, e.g. Legendre.
   * \tparam  return_t            Floating type specification. Default is double.
   * \retval  corner_vals         \c std::array of polynomial degrees containing \c std::array of
   *                              corner indices (the shape functions' derivatives are evaluated
   *                              at).
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <unsigned int max_poly_degree, typename return_t = double>
  static inline std::array<std::array<return_t, 2>, max_poly_degree + 1> shape_ders_at_bdrs()
  {
    std::array<return_t, 2> bdrs = {0., 1.};
    std::array<unsigned int, max_poly_degree + 1> poly_deg_index;
    for (unsigned int i = 0; i < poly_deg_index.size(); ++i)
      poly_deg_index[i] = i;

    return shape_der_eval<return_t, shape_t>(poly_deg_index, bdrs);
  }
};  // end of struct ShapeFun1D
