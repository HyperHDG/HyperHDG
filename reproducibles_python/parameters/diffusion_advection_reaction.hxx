#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <cmath>

template <unsigned int space_dimT, typename param_float_t = double>
struct TestParametersSinEllipt
{
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                               const param_float_t = 0.)
  {
    return 1. / cos(0.25 * M_PI * point[0]);
  }
  static SmallVec<space_dimT, param_float_t> advective_velocity(
    const Point<space_dimT, param_float_t>& point,
    const param_float_t = 0.)
  {
    SmallVec<space_dimT, param_float_t> advection;
    advection[0] = -0.25 * M_PI * sin(0.25 * M_PI * point[0]) + 1.;
    return inverse_diffusion_coeff(point) * advection;
  }
  static param_float_t reaction_rate(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 1.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t = 0.)
  {
    return sin(0.25 * M_PI * point[0]);
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t = 0.)
  {
    return sin(0.25 * M_PI * point[0]) + 0.25 * M_PI * cos(0.25 * M_PI * point[0]);
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t = 0.)
  {
    return analytic_result(point);
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t = 0.)
  {
    return analytic_result(point);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <unsigned int hyEdge_dimT>
struct HG
{
  template <unsigned int space_dimT, typename param_float_t = double>
  struct TestParametersQuadEllipt
  {
    static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
      1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
      14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    static constexpr std::array<unsigned int, 0U> neumann_nodes{};
    static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                                 const param_float_t = 0.)
    {
      return 1.;
    }
    static SmallVec<space_dimT, param_float_t> advective_velocity(
      const Point<space_dimT, param_float_t>& point,
      const param_float_t = 0.)
    {
      SmallVec<space_dimT, param_float_t> advection;
      for (unsigned int dim = 0; dim < space_dimT; ++dim)
        advection[dim] = (point[dim] - 2.) / analytic_result(point);
      return inverse_diffusion_coeff(point) * advection;
    }
    static param_float_t reaction_rate(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
    {
      return 1.;
    }
    static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                         const param_float_t = 0.)
    {
      param_float_t result = 0;
      for (unsigned int dim = 0; dim < space_dimT; ++dim)
        result -= (point[dim] - 2.) * (point[dim] - 2.);
      return result;
    }
    static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                         const param_float_t = 0.)
    {
      return 3. * (param_float_t)hyEdge_dimT + analytic_result(point);
    }
    static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                         const param_float_t = 0.)
    {
      return analytic_result(point);
    }
    static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                                 const param_float_t = 0.)
    {
      return analytic_result(point);
    }
    static param_float_t neumann_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t = 0.)
    {
      return 0.;
    }
  };
};
