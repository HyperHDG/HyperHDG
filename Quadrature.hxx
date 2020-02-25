#ifndef QUADRATURE_HXX
#define QUADRATURE_HXX

#include <HyAssert.hxx>

#include <array>
#include <cmath>


struct Gaussian
{
  /*!***********************************************************************************************
   * \brief   Calculate the amount of quadrature points at compile time.
   * 
   * Naive implementation to calculate the amount of needed quadrature points if a rule of accuracy
   * \c max_quad_degree is desired in \c local_dimensions dimensions using an orthogonal product of
   * Gaussian quadrature rules.
   * 
   * \param   max_quad_degree     Desired degree of accuracy.
   * \param   local_dimensions    Dimension of the underlying domain.
   * \retval  n_quad_points       Amount of needed quadrature points.
   * 
   * \authors   Guido Kanschat, University of Heidelberg, 2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2020.
   ************************************************************************************************/
  static constexpr unsigned int compute_n_quad_points
  (const unsigned int max_quad_degree, const unsigned int local_dimensions = 1)
  {
    unsigned int amount = 1, amount1D = 1;
    for ( ; 2 * amount1D - 1 < max_quad_degree; ++amount1D ) ;
    for ( unsigned int dim = 0; dim < local_dimensions; ++dim )  amount *= amount1D;
    return amount;
  }
  /*!***********************************************************************************************
   * \brief   Gaussian quadrature points on one-dimensional unit interval.
   * 
   * Returns the quadrature points of the quadrature rule with accuracy order \c max_quad_degree on
   * a one-dimensional unit interval \f$[0,1]\f$.
   * 
   * \tparam  max_quad_degree     Desired degree of accuracy.
   * \tparam  return_t            Floating type specification. Default is double.
   * \retval  quad_points         \c std::array containing the quadrature points.
   * 
   * \authors   Guido Kanschat, University of Heidelberg, 2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2020.
   ************************************************************************************************/
  template< unsigned int max_quad_degree, typename return_t = double >
  static inline std::array<return_t, compute_n_quad_points(max_quad_degree)> quad_points()
  {
    constexpr unsigned int n_points = compute_n_quad_points(max_quad_degree);
    static_assert( 1 <= n_points && n_points <= 9, "Amount of points needs to be smaller than 10!");
    std::array<return_t, n_points> quad_points;
    
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
                      -std::sqrt(5.-2.*std::sqrt(10./7.))/3. ,
                       0. ,
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
               "The number of points should equal the size of the array to be returned. In this "
               << "case the number of points is " << n_points << " and the size of the array is " <<
               quad_points.size() );

    // Transform quadrature points from [-1,1] -> [0,1]
    for (unsigned int index = 0; index < quad_points.size(); ++index)
      quad_points[index] = 0.5 * ( quad_points[index] + 1. );
    
    return quad_points;
  }
  /*!***********************************************************************************************
   * \brief   Gaussian quadrature weights on one-dimensional unit interval.
   * 
   * Returns the quadrature weights of the quadrature rule with accuracy order \c max_quad_degree on
   * a one-dimensional unit interval \f$[0,1]\f$.
   * 
   * \tparam  max_quad_degree     Desired degree of accuracy.
   * \tparam  return_t        Floating type specification. Default is double.
   * \retval  quad_weights        \c std::array containing the quadrature weights.
   * 
   * \authors   Guido Kanschat, University of Heidelberg, 2020.
   * \authors   Andreas Rupp, University of Heidelberg, 2020.
   ************************************************************************************************/
  template< unsigned int max_quad_degree, typename return_t = double >
  static inline std::array<return_t, compute_n_quad_points(max_quad_degree)> quad_weights()
  {
    constexpr unsigned int n_points = compute_n_quad_points(max_quad_degree);
    static_assert( 1 <= n_points && n_points <= 9, "Amount of points needs to be smaller than 10!");
    std::array<return_t, n_points> quad_weights;
    
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
               "The number of points should equal the size of the array to be returned. In this "
               << "case the number of points is " << n_points << " and the size of the array is " <<
               quad_weights.size() );

    // Transform quadrature points from [-1,1] -> [0,1]
    for (unsigned int index = 0; index < quad_weights.size(); ++index)  quad_weights[index] *= 0.5;
    
    return quad_weights;
  }
}; // end of struct Gaussian

// Shape functions & their derivatives evaluated at quadrature's points:

/*!*************************************************************************************************
 * \brief   Orthonormal shape functions evaluated at Gaussian quadrature points.
 * 
 * Returns the values of the orthonormal shape functions on \f$[0,1]\f$ of degree at most
 * \c max_poly_degree at the quadrature rule with accuracy order \c max_quad_degree on a
 * one-dimensional unit interval \f$[0,1]\f$.
 * 
 * \tparam  max_quad_degree     Desired degree of accuracy.
 * \tparam  max_poly_degree     Maximum degree of evaluated polynomials.
 * \tparam  return_t        Floating type specification. Default is double.
 * \retval  quad_vals           \c std::array of polynomial degrees containing \c std::array of 
 *                              quadrature points (the shape functions are evaluated at).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template
< 
  unsigned int max_poly_degree, unsigned int max_quad_degree,
  typename quadrature_t, typename shape_t, typename return_t = double
>
std::array
<
  std::array < return_t, quadrature_t::compute_n_quad_points(max_quad_degree) > , 
  max_poly_degree + 1
>
shape_fcts_at_quad_points()
{
  constexpr unsigned int n_points = quadrature_t::compute_n_quad_points(max_quad_degree);
  
  std::array<return_t, n_points> quad_points
    = quadrature_t::template quad_points<max_quad_degree, return_t>();
  std::array<unsigned int, max_poly_degree + 1> poly_deg_index;
  for (unsigned int i = 0; i < poly_deg_index.size(); ++i)  poly_deg_index[i] = i;
  
  return shape_fct_eval<return_t, shape_t>(poly_deg_index,quad_points);
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
 * \tparam  return_t        Floating type specification. Default is double.
 * \retval  quad_vals           \c std::array of polynomial degrees containing \c std::array of 
 *                              quadrature points (the shape functions' derivatives are evaluated).
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template
< 
  unsigned int max_poly_degree, unsigned int max_quad_degree,
  typename quadrature_t, typename shape_t, typename return_t = double
>
std::array
<
  std::array < return_t, quadrature_t::template compute_n_quad_points(max_quad_degree) > , 
  max_poly_degree + 1
>
shape_ders_at_quad_points()
{
  constexpr unsigned int n_points = quadrature_t::compute_n_quad_points(max_quad_degree);
  
  std::array<return_t, n_points> quad_points
    = quadrature_t::template quad_points<max_quad_degree,return_t>();
  std::array<unsigned int, max_poly_degree + 1> poly_deg_index;
  for (unsigned int i = 0; i < poly_deg_index.size(); ++i)  poly_deg_index[i] = i;
  
  return shape_der_eval<return_t, shape_t>(poly_deg_index,quad_points);
}

// General integrator class

template
< 
  unsigned int max_poly_degree, unsigned int max_quad_degree,
  typename quadrature_t, typename shape_t, typename return_t = double
>
class IntegratorTensorial
{
  private:
    const std::array < return_t, quadrature_t::compute_n_quad_points(max_quad_degree)>
      quad_points_, quad_weights_;
    const std::array
    <
      std::array < return_t, quadrature_t::compute_n_quad_points(max_quad_degree) > , 
      max_poly_degree + 1
    >
    shape_fcts_at_quad_, shape_ders_at_quad_;
  public:
    IntegratorTensorial()
    : quad_points_(quadrature_t::template quad_points<max_quad_degree,return_t>()),
      quad_weights_(quadrature_t::template quad_weights<max_quad_degree,return_t>()),
      shape_fcts_at_quad_(shape_fcts_at_quad_points
        <max_poly_degree, max_quad_degree, quadrature_t, shape_t, return_t>()),
      shape_ders_at_quad_(shape_ders_at_quad_points
        <max_poly_degree, max_quad_degree, quadrature_t, shape_t, return_t>())
    { }
    
    template < std::size_t size_fct, std::size_t size_der = 0 >
    return_t integrate1D
    ( std::array<unsigned int, size_fct> fct_ind,
      std::array<unsigned int, size_der> der_ind = { } ) const
    {
      return_t result = 0., helper;
      
      for (unsigned int q = 0; q < quad_weights_.size(); ++q)
      {
        helper = 1.;
        for (unsigned int i = 0; i < fct_ind.size(); ++i)
          helper *= shape_fcts_at_quad_[fct_ind[i]][q];
        for (unsigned int i = 0; i < der_ind.size(); ++i)
          helper *= shape_ders_at_quad_[der_ind[i]][q];
        result += quad_weights_[q] * helper;
      }
      
      return result;
    }
}; // end of class Integrator

#endif // end of ifndef QUADRATURE_HXX
