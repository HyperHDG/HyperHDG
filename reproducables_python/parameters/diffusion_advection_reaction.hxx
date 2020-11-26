#pragma once  // Ensure that file is included only once in a single compilation.

#include <cmath>

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TestParametersSinEllipt
{
  static constexpr double pi = acos(-1);
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                               const param_float_t time = 0.)
  {
    return 1. / cos(0.25 * pi * point[0]);
  }
  static SmallVec<space_dimT, param_float_t> advective_velocity(
    const Point<space_dimT, param_float_t>& point,
    const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> advection;
    advection[0] = -0.25 * pi * sin(0.25 * pi * point[0]) + 1.;
    return inverse_diffusion_coeff(point, time) * advection;
  }
  static param_float_t reaction_rate(const Point<space_dimT, param_float_t>& point,
                                     const param_float_t time = 0.)
  {
    return 1.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return sin(0.25 * pi * point[0]);
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return sin(0.25 * pi * point[0]) + 0.25 * pi * cos(0.25 * pi * point[0]);
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return analytic_result(point);
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return analytic_result(point);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& point,
                                     const param_float_t time = 0.)
  {
    return 0.;
  }
};

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT>
struct HG
{
  template <unsigned int space_dimT, typename param_float_t = double>
  struct TestParametersQuadEllipt
  {
    static constexpr double pi = acos(-1);
    static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
      1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
      14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    static constexpr std::array<unsigned int, 0U> neumann_nodes{};
    static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                                 const param_float_t time = 0.)
    {
      return 1.;
    }
    static SmallVec<space_dimT, param_float_t> advective_velocity(
      const Point<space_dimT, param_float_t>& point,
      const param_float_t time = 0.)
    {
      SmallVec<space_dimT, param_float_t> advection;
      for (unsigned int dim = 0; dim < space_dimT; ++dim)
        advection[dim] = (point[dim] - 2.) / analytic_result(point, time);
      return inverse_diffusion_coeff(point, time) * advection;
    }
    static param_float_t reaction_rate(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
    {
      return 1.;
    }
    static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                         const param_float_t time = 0.)
    {
      param_float_t result = 0;
      for (unsigned int dim = 0; dim < space_dimT; ++dim)
        result -= (point[dim] - 2.) * (point[dim] - 2.);
      return result;
    }

    static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                         const param_float_t time = 0.)
    {
      return 3. * (param_float_t)hyEdge_dimT + analytic_result(point, time);
    }

    static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                         const param_float_t time = 0.)
    {
      return analytic_result(point);
    }
    static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                                 const param_float_t time = 0.)
    {
      return analytic_result(point);
    }
    static param_float_t neumann_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
    {
      return 0.;
    }
  };
};
