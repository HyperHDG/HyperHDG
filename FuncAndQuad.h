/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 * 
 * Definition of test and trial functions together with Gauss quadrature rules.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef FUNC_AND_QUAD_H
#define FUNC_AND_QUAD_H

#include <array>

// Naive implementation finding the minimal amount of quadrature points to exactly integrate a polynomial of degree at most
// max_quad_degree in hyperedge_dim dimensions
constexpr const unsigned int compute_n_quad_points(const unsigned int max_quad_degree, const unsigned int local_dimensions = 1)
{
  unsigned int amount = 1, amount1D = 1;
  for ( ; 2 * amount1D - 1 < max_quad_degree; ++amount1D ) ;
  for ( unsigned int dim = 0; dim < local_dimensions; ++dim )  amount *= amount1D;
  return amount;
}


// Functions returning the value of the (derivative of the) one-dimensional trial function number "index" at point "x_value".
double trial_function_eval(const unsigned int index, const double x_value);
double deriv_of_trial_eval(const unsigned int index, const double x_value);


template<unsigned int max_quad_degree>
std::array<double, compute_n_quad_points(max_quad_degree)>
quadrature_points();

template<unsigned int max_quad_degree>
std::array<double, compute_n_quad_points(max_quad_degree)>
quadrature_weights();

template<unsigned int max_poly_degree, unsigned int max_quad_degree>
std::array< std::array<double, compute_n_quad_points(max_quad_degree)> , max_poly_degree + 1 >
trial_functions_at_quadrature_points();

template<unsigned int max_poly_degree, unsigned int max_quad_degree>
std::array< std::array<double, compute_n_quad_points(max_quad_degree)> , max_poly_degree + 1 >
derivs_of_trial_at_quadrature_points();

template<unsigned int max_poly_degree>
std::array< std::array<double, 2> , max_poly_degree + 1 >
trial_functions_at_boundaries();

template<unsigned int max_poly_degree>
std::array< std::array<double, 2> , max_poly_degree + 1 >
derivs_of_trial_at_boundaries();

#endif
