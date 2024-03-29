#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <cmath>

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
  struct DiffusionElliptic
  {
    static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
      1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
      14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
    static constexpr std::array<unsigned int, 0U> neumann_nodes{};
    static param_float_t diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                         const param_float_t = 0.)
    {
      return 1.;
    }
    static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                                 const param_float_t = 0.)
    {
      return 1.;
    }

    static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                         const param_float_t = 0.)
    {
      param_float_t result = 0;
      for (unsigned int dim = 0; dim < space_dimT; ++dim)
        result -= point[dim] * point[dim];
      return result;
    }

    static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                         const param_float_t = 0.)
    {
      return 2. * (param_float_t)hyEdge_dimT;
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
};
