/*!*************************************************************************************************
 * \file    FuncAndQuad.hxx
 * \brief   A namespace containing different functions that evaluate the trial functions at
 *          quadrature points.
 *
 * This namespace aims to ultimately provide the opportunity to choose from several types of trial
 * functions and different quadrature rules by choosing the appropriate \c namesapce. However, at
 * the moment only orthonormal trial functions and Gaussian quadrature are implemented resultin in
 * the currently very general name FuncQuad.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#ifndef FUNC_AND_QUAD_HXX
#define FUNC_AND_QUAD_HXX

#include <array>

/*!*************************************************************************************************
 * \brief   A namespace containing different functions that evaluate the trial functions at
 *          quadrature points.
 * 
 * \todo    Check, whether this construction makes sense.
 *
 * This namespace aims to ultimately provide the opportunity to choose from several types of trial
 * functions and different quadrature rules by choosing the appropriate \c namesapce. However, at
 * the moment only orthonormal trial functions and Gaussian quadrature are implemented resultin in
 * the currently very general name FuncQuad.
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
constexpr const unsigned int compute_n_quad_points(const unsigned int max_quad_degree,
                                                   const unsigned int local_dimensions = 1)
{
  unsigned int amount = 1, amount1D = 1;
  for ( ; 2 * amount1D - 1 < max_quad_degree; ++amount1D ) ;
  for ( unsigned int dim = 0; dim < local_dimensions; ++dim )  amount *= amount1D;
  return amount;
}

/*!*************************************************************************************************
 * \brief   Evaluate value of orthonormal trial function.
 * 
 * Evaluates the value of the \c index orthonormal, one-dimensional trial function on the reference
 * interval \f$[0,1]\f$ at abscissa \c x_value.
 * 
 * \param   index               Index of evaluated trial function.
 * \param   x_value             Abscissa of evaluated trial function.
 * \retval  fct_value           Evaluated value of trial function.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
double trial_function_eval(const unsigned int index, const double x_value);
/*!*************************************************************************************************
 * \brief   Evaluate value of the derivatibe of orthonormal trial function.
 * 
 * Evaluates the value of the derivative of the \c index orthonormal, one-dimensional trial function
 * on the reference interval \f$[0,1]\f$ at abscissa \c x_value.
 * 
 * \param   index               Index of evaluated trial function.
 * \param   x_value             Abscissa of evaluated trial function.
 * \retval  fct_value           Evaluated value of trial function's derivative.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
double deriv_of_trial_eval(const unsigned int index, const double x_value);
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
std::array<double, compute_n_quad_points(max_quad_degree)> quadrature_points();
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
std::array<double, compute_n_quad_points(max_quad_degree)> quadrature_weights();
/*!*************************************************************************************************
 * \brief   Orthonormal trial functions evaluated at Gaussian quadrature points.
 * 
 * Returns the values of the orthonormal trial functions on \f$[0,1]\f$ of degree at most
 * \c max_poly_degree at the quadrature rule with accuracy order \c max_quad_degree on a
 * one-dimensional unit interval \f$[0,1]\f$.
 * 
 * \tparam  max_quad_degree     Desired degree of accuracy.
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \retval  quad_vals           \c std::array of polynomial degrees containing \c std::array of 
 *                              quadrature points (the trial functions are evaluated at).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_poly_degree, unsigned int max_quad_degree>
std::array< std::array<double, compute_n_quad_points(max_quad_degree)> , max_poly_degree + 1 >
trial_functions_at_quadrature_points();
/*!*************************************************************************************************
 * \brief   Derivatives of orthonormal trial functions evaluated at Gaussian quadrature points.
 * 
 * Returns the values of the derivatives of orthonormal trial functions on \f$[0,1]\f$ of degree at
 * most \c max_poly_degree at the quadrature rule with accuracy order \c max_quad_degree on a
 * one-dimensional unit interval \f$[0,1]\f$.
 * 
 * \tparam  max_quad_degree     Desired degree of accuracy.
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \retval  quad_vals           \c std::array of polynomial degrees containing \c std::array of 
 *                              quadrature points (the trial functions' derivatives are evaluated).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_poly_degree, unsigned int max_quad_degree>
std::array< std::array<double, compute_n_quad_points(max_quad_degree)> , max_poly_degree + 1 >
derivs_of_trial_at_quadrature_points();
/*!*************************************************************************************************
 * \brief   Orthonormal trial functions evaluated at end points of unit interval.
 * 
 * Returns the values of the orthonormal trial functions on \f$[0,1]\f$ of degree at most
 * \c max_poly_degree at the value \f$0\f$ (at index 0) and at \f$1\f$ (at index 1).
 * 
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \retval  corner_vals         \c std::array of polynomial degrees containing \c std::array of 
 *                              corner indices (the trial functions are evaluated at).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_poly_degree>
std::array< std::array<double, 2> , max_poly_degree + 1 > trial_functions_at_boundaries();
/*!*************************************************************************************************
 * \brief   Derivatives of orthonormal trial functions evaluated at end points of unit interval.
 * 
 * Returns the values of the orthonormal trial functions' derivatives on \f$[0,1]\f$ of degree at
 * most \c max_poly_degree at the value \f$0\f$ (at index 0) and at \f$1\f$ (at index 1).
 * 
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \retval  corner_vals         \c std::array of polynomial degrees containing \c std::array of 
 *                              corner indices (the trial functions' derivatives are evaluated at).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<unsigned int max_poly_degree>
std::array< std::array<double, 2> , max_poly_degree + 1 > derivs_of_trial_at_boundaries();

} // end of namespace FuncQuad

#endif // end of ifndef FUNC_AND_QUAD_HXX
