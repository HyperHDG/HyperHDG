#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <cmath>

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TestParametersSinParab
{
  static constexpr std::array<unsigned int, 0U> boundary_nodes{};
  static constexpr std::array<unsigned int, 9U> dirichlet_nodes{1, 4, 7, 10, 13, 16, 19, 22, 25};
  static constexpr std::array<unsigned int, 9U> outflow_nodes{2, 5, 8, 11, 14, 17, 20, 23, 26};
  static SmallVec<space_dimT, param_float_t> velocity(const Point<space_dimT, param_float_t>&,
                                                      const param_float_t = 0.)
  {
    SmallVec<space_dimT, param_float_t> velocity;
    velocity[0] = 1. / M_PI;
    return velocity;
  }

  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return sin(M_PI * (point[0] + time));
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return (1. + M_PI) * cos(M_PI * (point[0] + time));
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
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
template <unsigned int space_dimT, typename param_float_t = double>
struct LeVeque
{
  static constexpr std::array<unsigned int, 26U> boundary_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> dirichlet_nodes{};
  static constexpr std::array<unsigned int, 0U> outflow_nodes{};
  static SmallVec<space_dimT, param_float_t> velocity(const Point<space_dimT, param_float_t>& point,
                                                      const param_float_t = 0.)
  {
    static_assert(space_dimT == 2, "Example is defined in two spatial dimensions.");
    SmallVec<space_dimT, param_float_t> velocity;
    velocity[0] = 0.5 - point[1];
    velocity[1] = point[0] - 0.5;
    return velocity;
  }

  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t = 0.)
  {
    static_assert(space_dimT == 2, "Example is defined in two spatial dimensions.");
    const param_float_t x = point[0], y = point[1], r = 0.15;

    Point<space_dimT, param_float_t> center;
    center[0] = 0.5;
    center[1] = 0.75;
    if (norm_2(point - center) <= r && (x <= 0.475 || x >= 0.525 || y >= 0.85))
      return 1.;
    center[0] = 0.5;
    center[1] = 0.25;
    if (norm_2(point - center) <= r)
      return 1. - norm_2(point - center) / r;
    center[0] = 0.25;
    center[1] = 0.5;
    if (norm_2(point - center) <= r)
      return 0.25 * (1. + cos(M_PI * norm_2(point - center) / r));
    return 0.;
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
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
template <unsigned int space_dimT, typename param_float_t = double>
struct injection
{
  static constexpr std::array<unsigned int, 0U> boundary_nodes{};
  static constexpr std::array<unsigned int, 1U> dirichlet_nodes{1};
  static constexpr std::array<unsigned int, 1U> outflow_nodes{2};
  static SmallVec<space_dimT, param_float_t> velocity(const auto&, const param_float_t = 0.)
  {
    SmallVec<space_dimT, param_float_t> velocity;
    velocity[0] = 1.;
    return velocity;
  }

  static param_float_t analytic_result(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 1.;
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};
