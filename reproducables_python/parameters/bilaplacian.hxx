#pragma once // Ensure that file is included only once in a single compilation.

#include <cmath>

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int space_dimT, typename param_float_t = double >
struct TestParametersSin
{
  static constexpr double pi = acos(-1);
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 26U> dirichlet_laplacian_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 0U> neumann_nodes {};
  static constexpr std::array<unsigned int, 0U> neumann_laplcian_nodes {};
  static param_float_t inverse_bilaplacian_coefficient( const Point<space_dimT,param_float_t>& pt )
  { return pi*pi; }
  static param_float_t analytic_result( const Point<space_dimT,param_float_t>& pt )
  { return cos(pi * pt[0]); }
  static param_float_t right_hand_side( const Point<space_dimT,param_float_t>& pt )
  { return pi * pi * cos(pi * pt[0]); }
  static param_float_t dirichlet_value( const Point<space_dimT,param_float_t>& pt )
  { return analytic_result(pt); }
  static param_float_t dirichlet_laplace_value( const Point<space_dimT,param_float_t>& pt )
  { return analytic_result(pt); }
  static param_float_t neumann_value( const Point<space_dimT,param_float_t>& pt )
  { return 0.; }
  static param_float_t neumann_laplace_value( const Point<space_dimT,param_float_t>& pt )
  { return 0.; }
};