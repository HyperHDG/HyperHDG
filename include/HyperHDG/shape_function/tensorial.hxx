#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/shape_function/one_dimensional.hxx>

#include <array>

namespace ShapeType
{
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
template <template <unsigned int> typename shape_t, unsigned int dimT, unsigned int poly_deg>
struct Tensorial
{
  using shape_fun_1d = shape_t<poly_deg>;
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_shape_fun() { return Hypercube<dimT>::pow(shape_fun_1d::n_shape_fun()); }

  static constexpr unsigned int dim() { return dimT; }

  static constexpr unsigned int degree() { return poly_deg; }

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
    static_assert(point_t::size() == dimT, "Point needs to have correct dimension.");

    if constexpr (dimT == 0)
      return (return_t)1.;

    return_t value = 0.;
    std::array<unsigned int, std::max(dimT, 1U)> index_dim = index_decompose(index,poly_deg+1);
    for (unsigned int dim = 0; dim < dimT; ++dim)
      value += shape_fun_1d::template fct_val<return_t>(index_dim[dim], point[dim]);

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
  static inline return_t der_val(const unsigned int index,
                                 const point_t& point,
                                 const unsigned int der_dim)
  {
    static_assert(point_t::size() == dimT, "Point needs to have correct dimension.");
    hy_assert(der_dim < dimT, "The derivative needs to be with respect to a valid dimension.");

    return_t value = 0.;
    std::array<unsigned int, std::max(dimT, 1U)> index_dim = index_decompose(index,poly_deg+1);
    for (unsigned int dim = 0; dim < dimT; ++dim)
      if (der_dim == dim)
        value += shape_fun_1d::template der_val<return_t>(index_dim[dim], point[dim]);
      else
        value += shape_fun_1d::template fct_val<return_t>(index_dim[dim], point[dim]);

    return value;
  }
};

}  // end of namespace ShapeType