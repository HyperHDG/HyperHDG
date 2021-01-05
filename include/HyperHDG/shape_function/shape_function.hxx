#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/shape_function/tensorial.hxx>

/*!*************************************************************************************************
 * \brief   Class that handles the evaluation of tensorial shape functions in tensor product points.
 *
 * \tparam  dimT              Local dimension of the shape function's domain or the point.
 * \tparam  return_t          The floating type that is returned.
 * \tparam  polynomial_degree Maximum polynomial degree that can be evaluated.
 * \tparam  abscissas_sizeT   Size of abszissas array.
 * \tparam  abscissa_float_t  Floating type for the abscissa values.
 * \tparam  shape_t           Type of 1D shapefunctions.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 * \authors   Simon Schmidt, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename shape_t>
struct ShapeFunction
{
  typedef shape_t shape_fun_t;
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_fun() { return shape_t::n_fun(); }

  static constexpr unsigned int dim() { return shape_t::dim(); }

  static constexpr unsigned int degree() { return shape_t::degree(); }
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
  template <typename return_t, typename point_t>
  static inline return_t fct_val(const unsigned int index, const point_t& point)
  {
    return shape_t::template fct_val<return_t>(index, point);
  }
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
  template <typename return_t, typename point_t>
  static inline return_t der_val(const unsigned int index,
                                 const point_t& point,
                                 const unsigned int der_dim)
  {
    return shape_t::template der_val<return_t>(index, point, der_dim);
  }

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
  template <typename return_t, typename point_t>
  static inline return_t grad_val(const unsigned int index, const point_t& point)
  {
    static_assert(return_t::size() == point_t::size(), "Gradient and point have same dimension.");
    static_assert(point_t::size() == dim(), "Gradient and dimension need to be equal.");

    return_t value;
    for (unsigned int der_dim = 0; der_dim < dim(); ++der_dim)
      value[der_dim] = der_val<return_t::value_type>(index, point, der_dim);
    return value;
  }
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
  template <typename return_t, typename point_t>
  static inline return_t div_val(const unsigned int index, const point_t& point)
  {
    static_assert(point_t::size() == dim(), "Gradient and dimension need to be equal.");

    return_t value = 0.;
    for (unsigned int der_dim = 0; der_dim < dim(); ++der_dim)
      value += der_val<return_t>(index, point, der_dim);
    return value;
  }

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
  template <typename return_t, typename coeffs_t, typename point_t>
  static inline return_t lin_comb_fct_val(const coeffs_t& coeffs, const point_t& point)
  {
    static_assert(point_t::size() == dim(), "Point needs to have correct dimension.");
    static_assert(coeffs_t::size() == n_fun(), "Coeffs size must coincide with poly space!");

    return_t value = 0.;
    for (unsigned int index = 0; index < coeffs.size(); ++index)
      value += coeffs[index] * fct_val<return_t>(index, point);
    return value;
  }
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
  template <typename return_t, typename coeffs_t, typename point_t>
  static inline return_t lin_comb_der_val(const coeffs_t& coeffs,
                                          const point_t& point,
                                          const unsigned int der_dim)
  {
    static_assert(point_t::size() == dim(), "Point needs to have correct dimension.");
    static_assert(coeffs_t::size() == n_fun(), "Coeffs size must coincide with poly space!");

    return_t value = 0.;
    for (unsigned int index = 0; index < coeffs.size(); ++index)
      value += coeffs[index] * der_val<return_t>(index, point, der_dim);
    return value;
  }
};