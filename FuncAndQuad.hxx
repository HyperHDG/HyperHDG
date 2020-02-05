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

#include <HyAssert.hxx>

#include <array>
#include <cmath>

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
 * \todo    This function needs to be inline or the linker does not work. Why?
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
inline double shape_fct_eval(const unsigned int index, const double x_val)
{
  hy_assert( 0 <= index && index <= 5 ,
             "The index of a shape function must be non-negative and smaller than or equal to 5 at "
             << "the moment. Your choice has been " << index << "." );
  hy_assert( 0. <= x_val && x_val <= 1. ,
             "The abscissa / x for which the shape function is evaluated has been set to be in the "
             << "closed interval [0,1]. Your choice has been " << x_val << "." );

  switch (index)
  {
    case 0: return 1.;
    case 1: return std::sqrt(3)*(1.-2.*x_val);
    case 2: return std::sqrt(5)*((6.*x_val-6.)*x_val+1.);
    case 3: return std::sqrt(7)*(((20.*x_val-30.)*x_val+12.)*x_val-1.);
    case 4: return std::sqrt(9)*((((70.*x_val-140.)*x_val+90.)*x_val-20.)*x_val+1.);
    case 5: return std::sqrt(11)*(((((252.*x_val-630.)*x_val+560.)*x_val-210.)*x_val+30.)*x_val-1.);
    default: hy_assert( 0 == 1 , "This shape function has not yet been implemented." );
  }

  hy_assert( 0 == 1 , "Something went wrong when evaluating a shape function. This message however "
                      << "is never supposed to appear. Something went wrong in the core program." );
  return 0.;
}
/*!*************************************************************************************************
 * \brief   Evaluate value of the derivative of orthonormal shape function.
 *
 * \todo    This function needs to be inline or the linker does not work. Why?
 * 
 * Evaluates the value of the derivative of the \c index orthonormal, one-dimensional shape function
 * on the reference interval \f$[0,1]\f$ at abscissa \c x_val.
 * 
 * \param   index               Index of evaluated shape function.
 * \param   x_val               Abscissa of evaluated shape function.
 * \retval  fct_value           Evaluated value of shape function's derivative.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
inline double shape_der_eval(const unsigned int index, const double x_val)
{
  hy_assert( 0 <= index && index <= 5 ,
             "The index of a shape function must be non-negative and smaller than or equal to 5 at "
             << "the moment. Your choice has been " << index << "." );
  hy_assert( 0. <= x_val && x_val <= 1. ,
             "The abscissa / x for which the shape function is evaluated has been set to be in the "
             << "closed interval [0,1]. Your choice has been " << x_val << "." );

  switch (index)
  {
    case 0: return 0.;
    case 1: return -std::sqrt(12);
    case 2: return std::sqrt(5)*(12.*x_val-6.);
    case 3: return std::sqrt(7)*((60.*x_val-60.)*x_val+12.);
    case 4: return std::sqrt(9)*(((280.*x_val-420.)*x_val+180.)*x_val-20.);
    case 5: return std::sqrt(11)*((((1260.*x_val-2520.)*x_val+1680.)*x_val-420.)*x_val+30.);
    default: hy_assert( 0 == 1 , "This shape function has not yet been implemented." );
  }
  
  hy_assert( 0 == 1 , "Something went wrong when evaluating a shape function. This message however "
                      << "is never supposed to appear. Something went wrong in the core program." );
  return 0.;
}
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
std::array<double, compute_n_quad_points(max_quad_degree)> quad_points()
{
  constexpr unsigned int n_points = compute_n_quad_points(max_quad_degree);
  static_assert( 1 <= n_points && n_points <= 9 , "Amount of points needs to be smaller than 10!");
  std::array<double, n_points> quad_points;
  
  if constexpr (n_points == 1)
    quad_points = { 0. };
  if constexpr (n_points == 2)
    quad_points = { -std::sqrt(1./3.) , std::sqrt(1./3.) };
  if constexpr (n_points == 3)
    quad_points = { -std::sqrt(3./5.) , 0. , std::sqrt(3./5.) };
  if constexpr (n_points == 4)
    quad_points = { -std::sqrt(3./7.+2./7.*std::sqrt(6./5.)) , 
                    -std::sqrt(3./7.-2./7.*std::sqrt(6./5.)) ,
                     std::sqrt(3./7.-2./7.*std::sqrt(6./5.)) ,
                     std::sqrt(3./7.+2./7.*std::sqrt(6./5.)) };
  if constexpr (n_points == 5)
    quad_points = { -std::sqrt(5.+2.*std::sqrt(10./7.))/3. ,
                    -std::sqrt(5.-2.*std::sqrt(10./7.))/3. , 0. ,
                     std::sqrt(5.-2.*std::sqrt(10./7.))/3. ,
                     std::sqrt(5.+2.*std::sqrt(10./7.))/3. };
  if constexpr (n_points == 6)
    quad_points = { 0.6612093864662645, -0.6612093864662645, -0.2386191860831969,
                    0.2386191860831969, -0.9324695142031521,  0.9324695142031521 };
  if constexpr (n_points == 7)
    quad_points = { 0.0000000000000000,  0.4058451513773972, -0.4058451513773972,
                   -0.7415311855993945,  0.7415311855993945, -0.9491079123427585,
                    0.9491079123427585 };
  if constexpr (n_points == 8)
    quad_points = {-0.1834346424956498,  0.1834346424956498, -0.5255324099163290,
                    0.5255324099163290, -0.7966664774136267,  0.7966664774136267,
                   -0.9602898564975363,  0.9602898564975363 };
  if constexpr (n_points == 9)
    quad_points = { 0.0000000000000000, -0.8360311073266358,  0.8360311073266358,
                   -0.9681602395076261,  0.9681602395076261, -0.3242534234038089,
                    0.3123470770400029,  0.2606106964029354,  0.2606106964029354 };
  
  hy_assert( n_points == quad_points.size() ,
             "The number of points should equal the size of the array to be returned. In this case "
             << "the number of points is " << n_points << " and the size of the array is " <<
             quad_points.size() );

  // Transform quadrature points from [-1,1] -> [0,1]
  for (unsigned int index = 0; index < quad_points.size(); ++index)
    quad_points[index] = 0.5 * ( quad_points[index] + 1. );

  return quad_points;
}
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
std::array<double, compute_n_quad_points(max_quad_degree)> quad_weights()
{
  constexpr unsigned int n_points = compute_n_quad_points(max_quad_degree);
  static_assert( 1 <= n_points && n_points <= 9 , "Amount of points needs to be smaller than 10!");
  std::array<double, n_points> quad_weights;
  
  if constexpr (n_points == 1)
    quad_weights = { 2. };
  if constexpr (n_points == 2)
    quad_weights = { 1. , 1. };
  if constexpr (n_points == 3)
    quad_weights = { 5./9. , 8./9. , 5./9. };
  if constexpr (n_points == 4)
    quad_weights = { 1./36.*(18. - std::sqrt(30.)) , 1./36.*(18. + std::sqrt(30.)) ,
                     1./36.*(18. + std::sqrt(30.)) , 1./36.*(18. - std::sqrt(30.)) };
  if constexpr (n_points == 5)
    quad_weights = { 1./900.*(322.-13.*std::sqrt(70.)) , 1./900.*(322.+13.*std::sqrt(70.)) ,
                     1./900.*(322.+190.) ,
                     1./900.*(322.+13.*std::sqrt(70.)) , 1./900.*(322.-13.*std::sqrt(70.)) };
  if constexpr (n_points == 6)
    quad_weights = { 0.3607615730481386,  0.3607615730481386,  0.4679139345726910,
                     0.4679139345726910,  0.1713244923791704,  0.1713244923791700 };
  if constexpr (n_points == 7)
    quad_weights = { 0.4179591836734694,  0.3818300505051189,  0.3818300505051189,
                     0.2797053914892766,  0.2797053914892766,  0.1294849661688697,
                     0.1294849661688697 };
  if constexpr (n_points == 8)
    quad_weights = { 0.3626837833783620,  0.3626837833783620,  0.3137066458778873,
                     0.3137066458778873,  0.2223810344533745,  0.2223810344533745,
                     0.1012285362903763,  0.1012285362903763 };
  if constexpr (n_points == 9)
    quad_weights = { 0.3302393550012598,  0.1806481606948574,  0.1806481606948574,
                     0.0812743883615744,  0.0812743883615744,  0.3123470770400029,
                     0.3123470770400029,  0.2606106964029354,  0.2606106964029354 };

  hy_assert( n_points == quad_weights.size() ,
             "The number of points should equal the size of the array to be returned. In this case "
             << "the number of points is " << n_points << " and the size of the array is " <<
             quad_weights.size() );

  // Transform quadrature points from [-1,1] -> [0,1]
  for (unsigned int index = 0; index < quad_weights.size(); ++index)  quad_weights[index] *= 0.5;

  return quad_weights;
}
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
shape_fcts_at_quad_points()
{
  constexpr unsigned int n_points = compute_n_quad_points(max_quad_degree);
  static_assert( 1 <= max_poly_degree && max_poly_degree <= 5 , "Maximum polynomial degree is 5!");
  static_assert( 1 <= n_points && n_points <= 9 , "Amount of points needs to be smaller than 10!");
  
  std::array<double, n_points> quad_points = FuncQuad::quad_points<max_quad_degree>();
  std::array< std::array<double, n_points> , max_poly_degree + 1 > fct_val;
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
    for (unsigned int j = 0; j < n_points; ++j)
      fct_val[i][j] = FuncQuad::shape_fct_eval(i, quad_points[j]);
  
  return fct_val;
}
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
shape_ders_at_quad_points()
{
  constexpr unsigned int n_points = compute_n_quad_points(max_quad_degree);
  static_assert( 1 <= max_poly_degree && max_poly_degree <= 5 , "Maximum polynomial degree is 5!");
  static_assert( 1 <= n_points && n_points <= 9 , "Amount of points needs to be smaller than 10!");
  
  std::array<double, n_points> quad_points = FuncQuad::quad_points<max_quad_degree>();
  std::array< std::array<double, n_points> , max_poly_degree + 1 > der_val;
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
    for (unsigned int j = 0; j < n_points; ++j)
      der_val[i][j] = FuncQuad::shape_der_eval(i, quad_points[j]);
  
  return der_val;
}
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
std::array< std::array<double, 2> , max_poly_degree + 1 > shape_fcts_at_bdrs()
{
  static_assert( 1 <= max_poly_degree && max_poly_degree <= 5 , "Maximum polynomial degree is 5!");
  std::array< std::array<double, 2> , max_poly_degree + 1 > fct_val;
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
  {
    fct_val[i][0] = FuncQuad::shape_fct_eval(i, 0.);
    fct_val[i][1] = FuncQuad::shape_fct_eval(i, 1.);
  }
  
  return fct_val;
}
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
std::array< std::array<double, 2> , max_poly_degree + 1 > shape_ders_at_bdrs()
{
  static_assert( 1 <= max_poly_degree && max_poly_degree <= 5 , "Maximum polynomial degree is 5!");
  std::array< std::array<double, 2> , max_poly_degree + 1 > der_val;
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
  {
    der_val[i][0] = FuncQuad::shape_der_eval(i, 0.);
    der_val[i][1] = FuncQuad::shape_der_eval(i, 1.);
  }
  
  return der_val;
}

} // end of namespace FuncQuad

#endif // end of ifndef FUNC_AND_QUAD_HXX
