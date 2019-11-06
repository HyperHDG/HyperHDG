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

#include <vector>

double trial_function_eval(const unsigned int index, const double x_value);
double deriv_of_trial_eval(const unsigned int index, const double x_value);
std::vector<double> quadrature_points(const unsigned int num_of_points);
std::vector<double> quadrature_weights(const unsigned int num_of_points);
std::vector< std::vector<double> >  trial_functions_at_quadrature_points
  (const unsigned int max_poly_degree, const unsigned int num_of_points);
std::vector< std::vector<double> > derivs_of_trial_at_quadrature_points
  (const unsigned int max_poly_degree, const unsigned int num_of_points);
std::vector< std::vector<double> >  trial_functions_at_boundaries(const unsigned int max_poly_degree);
std::vector< std::vector<double> > derivs_of_trial_at_boundaries(const unsigned int max_poly_degree);

#endif
