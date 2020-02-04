/*!*************************************************************************************************
 * \file    FuncAndQuad.hxx
 * \brief   Contains functions to evaluate shape functions (e.g. at quadrature points).
 *
 * This file provides several functions to evaluate shape functions (and their derivatives) that are
 * L^2 orthonormal with respect to the unit interval \f$[0,1]\f$. (Obviously, this does not hold for
 * the derivatives.)
 * 
 * From these one-dimensional shape frunctions, multi-dimensional shape functions are constructed by
 * evaluating tensor / dyadic products. Thus, also the multi-dimensional shape functions are L^2
 * orthonormal with respect to the unit hypercube  \f$[0,1]^d\f$, where \f$d\f$ denotes the spatial
 * dimension.
 * 
 * Moreover, Gaussian quadrature rules are provided that contain at most nine quadrature points in
 * one dimension (and therefore are exact for polynomials of degree at most 17). Exploiting the
 * tensor structure of the reference hypercube allows to also extend these rules to be applicable
 * in several spatial dimensions.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#ifndef FUNC_AND_QUAD_HXX
#define FUNC_AND_QUAD_HXX

#include <array>

/*!*************************************************************************************************
 * \brief   Contains functions to evaluate shape functions (e.g. at quadrature points).
 *
 * This namespace provides several functions to evaluate shape functions (and their derivatives)
 * that are L^2 orthonormal with respect to the unit interval \f$[0,1]\f$. (Obviously, this does not
 * hold for the derivatives.)
 * 
 * From these one-dimensional shape frunctions, multi-dimensional shape functions are constructed by
 * evaluating tensor / dyadic products. Thus, also the multi-dimensional shape functions are L^2
 * orthonormal with respect to the unit hypercube  \f$[0,1]^d\f$, where \f$d\f$ denotes the spatial
 * dimension.
 * 
 * Moreover, Gaussian quadrature rules are provided that contain at most nine quadrature points in
 * one dimension (and therefore are exact for polynomials of degree at most 17). Exploiting the
 * tensor structure of the reference hypercube allows to also extend these rules to be applicable
 * in several spatial dimensions.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
namespace FuncQuad
{

/*!*************************************************************************************************
 * \brief   Calculate the amount of quadrature points at compile time.
 * 
 * Naive implementation to calculate the amount of needed quadrature points when a rule of accuracy
 * \c max_quad_degree is desired in \c local_dimensions dimensions using an orthogonal product of
 * Gaussian quadrature rules.
 * 
 * \param   max_quad_degree     Desired degree of accuracy.
 * \param   local_dimensions    Dimension of the underlying domain.
 * \retval  n_quad_points       Amount of needed quadrature points.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
constexpr const unsigned int compute_n_quad_points
(const unsigned int max_quad_degree, const unsigned int local_dimensions = 1)
{
  unsigned int amount = 1, amount1D = 1;
  for ( ; 2 * amount1D - 1 < max_quad_degree; ++amount1D ) ;
  for ( unsigned int dim = 0; dim < local_dimensions; ++dim )  amount *= amount1D;
  return amount;
}

/*!*************************************************************************************************
 * \brief   Evaluate value of orthonormal shape function.
 * 
 * Evaluates the value of the \c index orthonormal, one-dimensional shape function on the reference
 * interval \f$[0,1]\f$ at abscissa \c x_value.
 * 
 * \param   index               Index of evaluated shape function.
 * \param   x_value             Abscissa of evaluated shape function.
 * \retval  fct_value           Evaluated value of shape function.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
double shape_fct_eval(const unsigned int index, const double x_value);
/*!*************************************************************************************************
 * \brief   Evaluate value of the derivative of orthonormal shape function.
 * 
 * Evaluates the value of the derivative of the \c index orthonormal, one-dimensional shape function
 * on the reference interval \f$[0,1]\f$ at abscissa \c x_value.
 * 
 * \param   index               Index of evaluated shape function.
 * \param   x_value             Abscissa of evaluated shape function.
 * \retval  fct_value           Evaluated value of shape function's derivative.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
double shape_der_eval(const unsigned int index, const double x_value);
/*!*************************************************************************************************
 * \brief   Gaussian quadrature points on one-dimensional unit interval.
 * 
 * Returns the quadrature points of the quadrature rule with accuracy order \c max_quad_degree on a
 * one-dimensional unit interval \f$[0,1]\f$.
 * 
 * \tparam  max_quad_degree     Desired degree of accuracy.
 * \retval  quad_points         \c std::array containing the quadrature points.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_quad_degree>
std::array<double, compute_n_quad_points(max_quad_degree)> quad_points();
/*!*************************************************************************************************
 * \brief   Gaussian quadrature weights on one-dimensional unit interval.
 * 
 * Returns the quadrature weights of the quadrature rule with accuracy order \c max_quad_degree on a
 * one-dimensional unit interval \f$[0,1]\f$.
 * 
 * \tparam  max_quad_degree     Desired degree of accuracy.
 * \retval  quad_weights        \c std::array containing the quadrature weights.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_quad_degree>
std::array<double, compute_n_quad_points(max_quad_degree)> quad_weights();
/*!*************************************************************************************************
 * \brief   Orthonormal shape functions evaluated at Gaussian quadrature points.
 * 
 * Returns the values of the orthonormal shape functions on \f$[0,1]\f$ of degree at most
 * \c max_poly_degree at the quadrature rule with accuracy order \c max_quad_degree on a
 * one-dimensional unit interval \f$[0,1]\f$.
 * 
 * \tparam  max_quad_degree     Desired degree of accuracy.
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \retval  quad_vals           \c std::array of polynomial degrees containing \c std::array of 
 *                              quadrature points (the shape functions are evaluated at).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_poly_degree, unsigned int max_quad_degree>
std::array< std::array<double, compute_n_quad_points(max_quad_degree)> , max_poly_degree + 1 >
shape_fcts_at_quad_points();
/*!*************************************************************************************************
 * \brief   Derivatives of orthonormal shape functions evaluated at Gaussian quadrature points.
 * 
 * Returns the values of the derivatives of orthonormal shape functions on \f$[0,1]\f$ of degree at
 * most \c max_poly_degree at the quadrature rule with accuracy order \c max_quad_degree on a
 * one-dimensional unit interval \f$[0,1]\f$.
 * 
 * \tparam  max_quad_degree     Desired degree of accuracy.
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \retval  quad_vals           \c std::array of polynomial degrees containing \c std::array of 
 *                              quadrature points (the shape functions' derivatives are evaluated).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_poly_degree, unsigned int max_quad_degree>
std::array< std::array<double, compute_n_quad_points(max_quad_degree)> , max_poly_degree + 1 >
shape_ders_at_quad_points();
/*!*************************************************************************************************
 * \brief   Orthonormal shape functions evaluated at end points of unit interval.
 * 
 * Returns the values of the orthonormal shape functions on \f$[0,1]\f$ of degree at most
 * \c max_poly_degree at the value \f$0\f$ (at index 0) and at \f$1\f$ (at index 1).
 * 
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \retval  corner_vals         \c std::array of polynomial degrees containing \c std::array of 
 *                              corner indices (the shape functions are evaluated at).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_poly_degree>
std::array< std::array<double, 2> , max_poly_degree + 1 > shape_fcts_at_bdrs();
/*!*************************************************************************************************
 * \brief   Derivatives of orthonormal shape functions evaluated at end points of unit interval.
 * 
 * Returns the values of the orthonormal shape functions' derivatives on \f$[0,1]\f$ of degree at
 * most \c max_poly_degree at the value \f$0\f$ (at index 0) and at \f$1\f$ (at index 1).
 * 
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \retval  corner_vals         \c std::array of polynomial degrees containing \c std::array of 
 *                              corner indices (the shape functions' derivatives are evaluated at).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_poly_degree>
std::array< std::array<double, 2> , max_poly_degree + 1 > shape_ders_at_bdrs();

} // end of namespace FuncQuad

#endif // end of ifndef FUNC_AND_QUAD_HXX
