#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/tensorial_shape_fun.hxx>

#include <array>
#include <cmath>
//#include <numeric>

namespace Quadrature
{
template <unsigned int quad_deg>
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
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  static constexpr unsigned int n_quad_points()
  {
    unsigned int amount = 1;
    for (; 2 * amount - 1 < quad_deg; ++amount)
      ;
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
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t = double>
  static inline std::array<return_t, n_quad_points()> quad_points()
  {
    constexpr unsigned int n_points = n_quad_points();
    static_assert(1 <= n_points && n_points <= 9, "Amount of points needs to be smaller than 10!");
    std::array<return_t, n_points> quad_points;

    if constexpr (n_points == 1)
      quad_points = {static_cast<return_t>(0.)};
    if constexpr (n_points == 2)
      quad_points = {static_cast<return_t>(-std::sqrt(1. / 3.)),
                     static_cast<return_t>(std::sqrt(1. / 3.))};
    if constexpr (n_points == 3)
      quad_points = {static_cast<return_t>(-std::sqrt(3. / 5.)), static_cast<return_t>(0.),
                     static_cast<return_t>(std::sqrt(3. / 5.))};
    if constexpr (n_points == 4)
      quad_points = {static_cast<return_t>(-std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.))),
                     static_cast<return_t>(-std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.))),
                     static_cast<return_t>(std::sqrt(3. / 7. - 2. / 7. * std::sqrt(6. / 5.))),
                     static_cast<return_t>(std::sqrt(3. / 7. + 2. / 7. * std::sqrt(6. / 5.)))};
    if constexpr (n_points == 5)
      quad_points = {static_cast<return_t>(-std::sqrt(5. + 2. * std::sqrt(10. / 7.)) / 3.),
                     static_cast<return_t>(-std::sqrt(5. - 2. * std::sqrt(10. / 7.)) / 3.),
                     static_cast<return_t>(0.),
                     static_cast<return_t>(std::sqrt(5. - 2. * std::sqrt(10. / 7.)) / 3.),
                     static_cast<return_t>(std::sqrt(5. + 2. * std::sqrt(10. / 7.)) / 3.)};
    if constexpr (n_points == 6)
      quad_points = {
        static_cast<return_t>(0.6612093864662645),  static_cast<return_t>(-0.6612093864662645),
        static_cast<return_t>(-0.2386191860831969), static_cast<return_t>(0.2386191860831969),
        static_cast<return_t>(-0.9324695142031521), static_cast<return_t>(0.9324695142031521)};
    if constexpr (n_points == 7)
      quad_points = {
        static_cast<return_t>(0.0000000000000000),  static_cast<return_t>(0.4058451513773972),
        static_cast<return_t>(-0.4058451513773972), static_cast<return_t>(-0.7415311855993945),
        static_cast<return_t>(0.7415311855993945),  static_cast<return_t>(-0.9491079123427585),
        static_cast<return_t>(0.9491079123427585)};
    if constexpr (n_points == 8)
      quad_points = {
        static_cast<return_t>(-0.1834346424956498), static_cast<return_t>(0.1834346424956498),
        static_cast<return_t>(-0.5255324099163290), static_cast<return_t>(0.5255324099163290),
        static_cast<return_t>(-0.7966664774136267), static_cast<return_t>(0.7966664774136267),
        static_cast<return_t>(-0.9602898564975363), static_cast<return_t>(0.9602898564975363)};
    if constexpr (n_points == 9)
      quad_points = {
        static_cast<return_t>(0.0000000000000000), static_cast<return_t>(-0.8360311073266358),
        static_cast<return_t>(0.8360311073266358), static_cast<return_t>(-0.9681602395076261),
        static_cast<return_t>(0.9681602395076261), static_cast<return_t>(-0.3242534234038089),
        static_cast<return_t>(0.3123470770400029), static_cast<return_t>(0.2606106964029354),
        static_cast<return_t>(0.2606106964029354)};

    hy_assert(n_points == quad_points.size(),
              "The number of points should equal the size of the array to be returned. In this "
                << "case the number of points is " << n_points << " and the size of the array is "
                << quad_points.size());

    // Transform quadrature points from [-1,1] -> [0,1]
    for (unsigned int index = 0; index < quad_points.size(); ++index)
      quad_points[index] = 0.5 * (quad_points[index] + 1.);

    return quad_points;
  }
  /*!***********************************************************************************************
   * \brief   Gaussian quadrature weights on one-dimensional unit interval.
   *
   * Returns the quadrature weights of the quadrature rule with accuracy order \c max_quad_degree on
   * a one-dimensional unit interval \f$[0,1]\f$.
   *
   * \tparam  max_quad_degree     Desired degree of accuracy.
   * \tparam  return_t            Floating type specification. Default is double.
   * \retval  quad_weights        \c std::array containing the quadrature weights.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2020.
   * \authors   Andreas Rupp, Heidelberg University, 2020.
   ************************************************************************************************/
  template <typename return_t = double>
  static inline std::array<return_t, n_quad_points()> quad_weights()
  {
    constexpr unsigned int n_points = n_quad_points();
    static_assert(1 <= n_points && n_points <= 9, "Amount of points needs to be smaller than 10!");
    std::array<return_t, n_points> quad_weights;

    if constexpr (n_points == 1)
      quad_weights = {static_cast<return_t>(2.)};
    if constexpr (n_points == 2)
      quad_weights = {static_cast<return_t>(1.), static_cast<return_t>(1.)};
    if constexpr (n_points == 3)
      quad_weights = {static_cast<return_t>(5. / 9.), static_cast<return_t>(8. / 9.),
                      static_cast<return_t>(5. / 9.)};
    if constexpr (n_points == 4)
      quad_weights = {static_cast<return_t>(1. / 36. * (18. - std::sqrt(30.))),
                      static_cast<return_t>(1. / 36. * (18. + std::sqrt(30.))),
                      static_cast<return_t>(1. / 36. * (18. + std::sqrt(30.))),
                      static_cast<return_t>(1. / 36. * (18. - std::sqrt(30.)))};
    if constexpr (n_points == 5)
      quad_weights = {static_cast<return_t>(1. / 900. * (322. - 13. * std::sqrt(70.))),
                      static_cast<return_t>(1. / 900. * (322. + 13. * std::sqrt(70.))),
                      static_cast<return_t>(1. / 900. * (322. + 190.)),
                      static_cast<return_t>(1. / 900. * (322. + 13. * std::sqrt(70.))),
                      static_cast<return_t>(1. / 900. * (322. - 13. * std::sqrt(70.)))};
    if constexpr (n_points == 6)
      quad_weights = {
        static_cast<return_t>(0.3607615730481386), static_cast<return_t>(0.3607615730481386),
        static_cast<return_t>(0.4679139345726910), static_cast<return_t>(0.4679139345726910),
        static_cast<return_t>(0.1713244923791704), static_cast<return_t>(0.1713244923791700)};
    if constexpr (n_points == 7)
      quad_weights = {
        static_cast<return_t>(0.4179591836734694), static_cast<return_t>(0.3818300505051189),
        static_cast<return_t>(0.3818300505051189), static_cast<return_t>(0.2797053914892766),
        static_cast<return_t>(0.2797053914892766), static_cast<return_t>(0.1294849661688697),
        static_cast<return_t>(0.1294849661688697)};
    if constexpr (n_points == 8)
      quad_weights = {
        static_cast<return_t>(0.3626837833783620), static_cast<return_t>(0.3626837833783620),
        static_cast<return_t>(0.3137066458778873), static_cast<return_t>(0.3137066458778873),
        static_cast<return_t>(0.2223810344533745), static_cast<return_t>(0.2223810344533745),
        static_cast<return_t>(0.1012285362903763), static_cast<return_t>(0.1012285362903763)};
    if constexpr (n_points == 9)
      quad_weights = {
        static_cast<return_t>(0.3302393550012598), static_cast<return_t>(0.1806481606948574),
        static_cast<return_t>(0.1806481606948574), static_cast<return_t>(0.0812743883615744),
        static_cast<return_t>(0.0812743883615744), static_cast<return_t>(0.3123470770400029),
        static_cast<return_t>(0.3123470770400029), static_cast<return_t>(0.2606106964029354),
        static_cast<return_t>(0.2606106964029354)};

    hy_assert(n_points == quad_weights.size(),
              "The number of points should equal the size of the array to be returned. In this "
                << "case the number of points is " << n_points << " and the size of the array is "
                << quad_weights.size());

    // Transform quadrature points from [-1,1] -> [0,1]
    for (unsigned int index = 0; index < quad_weights.size(); ++index)
      quad_weights[index] *= 0.5;

    return quad_weights;
  }
};  // end of struct Gaussian

}  // end of namespace Quadrature