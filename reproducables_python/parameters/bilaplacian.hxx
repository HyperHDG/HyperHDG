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
  
  static param_float_t inverse_bilaplacian_coefficient
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return pi; }
  
  static param_float_t analytic_result
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return sin(pi * point[0]); }
  
  static param_float_t right_hand_side
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return pi * pi * pi * sin(pi * point[0]); }
  
  static param_float_t dirichlet_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point); }
  static param_float_t dirichlet_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return pi * analytic_result(point); }
  
  static param_float_t initial
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  static param_float_t initial_laplace
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  
  static param_float_t neumann_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  static param_float_t neumann_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
};

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int space_dimT, typename param_float_t = double >
struct TestParametersSinParab
{
  static constexpr double pi = acos(-1);
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 26U> dirichlet_laplacian_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 0U> neumann_nodes {};
  static constexpr std::array<unsigned int, 0U> neumann_laplcian_nodes {};
  
  static param_float_t inverse_bilaplacian_coefficient
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return pi; }
  
  static param_float_t analytic_result
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return sin(pi * (point[0] + time)); }
  
  static param_float_t right_hand_side
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return pi * pi * pi * sin(pi * (point[0] + time)) + pi * cos(pi * (point[0] + time)); }
  
  static param_float_t dirichlet_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  static param_float_t dirichlet_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  
  static param_float_t initial
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  static param_float_t initial_laplace
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  
  static param_float_t neumann_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  static param_float_t neumann_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
};

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int space_dimT, typename param_float_t = double >
struct TestParametersTrivParab
{
  static constexpr double pi = acos(-1);
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 26U> dirichlet_laplacian_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 0U> neumann_nodes {};
  static constexpr std::array<unsigned int, 0U> neumann_laplcian_nodes {};
  
  static param_float_t inverse_bilaplacian_coefficient
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 1.; }
  
  static param_float_t analytic_result
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 1. - time; }
  
  static param_float_t right_hand_side
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return -1.; }
  
  static param_float_t dirichlet_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  static param_float_t dirichlet_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  
  static param_float_t initial
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  static param_float_t initial_laplace
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  
  static param_float_t neumann_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  static param_float_t neumann_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
};

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int space_dimT, typename param_float_t = double >
struct TestParametersHomo
{
  static constexpr double pi = acos(-1);
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 26U> dirichlet_laplacian_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 0U> neumann_nodes {};
  static constexpr std::array<unsigned int, 0U> neumann_laplcian_nodes {};
  
  static param_float_t inverse_bilaplacian_coefficient
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 1.; }
  
  static param_float_t analytic_result
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  
  static param_float_t right_hand_side
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  
  static param_float_t dirichlet_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  static param_float_t dirichlet_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return analytic_result(point, time); }
  
  static param_float_t initial
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  {
    param_float_t approx = 1.;
    for (unsigned int dim = 0; dim < space_dimT; ++dim)
      approx *= sin(pi * point[dim]) + 0*1e-6 * ((std::rand() % 201) - 100); // Random noise!
    return approx;
  }
  static param_float_t initial_laplace
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { 
    param_float_t approx = 1.;
    for (unsigned int dim = 0; dim < space_dimT; ++dim)
      approx *= pi * pi * sin(pi * point[dim]) + 0*1e-6 * ((std::rand() % 201) - 100); // Random noise!
    return approx;
  }
  
  static param_float_t neumann_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
  static param_float_t neumann_laplace_value
  ( const Point<space_dimT,param_float_t>& point , const param_float_t time = 0. )
  { return 0.; }
};
