/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "FuncAndQuad.h"
#include <algorithm>
#include <cmath>
#include <cassert>

using namespace std;


double trial_function_eval(const unsigned int index, const double x_value)
{
  assert(0 <= index && index <= 5);
  assert(0. <= x_value && x_value <= 1.);
  switch (index)
  {
    case 0:  return 1.;
    case 1:  return sqrt(3) * (1. - 2. * x_value);
    case 2:  return sqrt(5) * ( (6. * x_value - 6.) * x_value + 1. );
    case 3:  return sqrt(7) * ( ( (20. * x_value - 30.) * x_value + 12. ) * x_value - 1. );
    case 4:  return sqrt(9) * ( ( ( (70. * x_value - 140.) * x_value + 90. ) * x_value - 20. ) * x_value + 1. );
    case 5:  return sqrt(11)* ( ( ( ( (252. * x_value - 630.) * x_value + 560. ) * x_value - 210. ) * x_value + 30. ) * x_value - 1. );
    default: assert(0 == 1);
  }
  assert(0 == 1);
  return 0.;
}


double deriv_of_trial_eval(const unsigned int index, const double x_value)
{
  assert(0 <= index && index <= 5);
  assert(0. <= x_value && x_value <= 1.);
  switch (index)
  {
    case 0:  return 0.;
    case 1:  return -sqrt(12);
    case 2:  return sqrt(5) * ( 12. * x_value - 6. );
    case 3:  return sqrt(7) * ( (60. * x_value - 60.) * x_value + 12. );
    case 4:  return sqrt(9) * ( ( (280. * x_value - 420.) * x_value + 180. ) * x_value - 20. );
    case 5:  return sqrt(11)* ( ( ( (1260. * x_value - 2520.) * x_value + 1680. ) * x_value - 420. ) * x_value + 30. );
    default: assert(0 == 1);
  }
  assert(0 == 1);
  return 0.;
}


vector<double> quadrature_points(const unsigned int num_of_points)
{
  assert(1 <= num_of_points && num_of_points <= 9);
  vector<double> quad_points;
  
  switch (num_of_points)
  {
    case 1: quad_points = { 0. };
            break;
    case 2: quad_points = { -sqrt(1./3.) , sqrt(1./3.) };
            break;
    case 3: quad_points = { -sqrt(3./5.) , 0. , sqrt(3./5.) };
            break;
    case 4: quad_points = { -sqrt(3./7.+2./7.*sqrt(6./5.)) , -sqrt(3./7.-2./7.*sqrt(6./5.)) ,
                             sqrt(3./7.-2./7.*sqrt(6./5.)) ,  sqrt(3./7.+2./7.*sqrt(6./5.)) };
            break;
    case 5: quad_points = { -sqrt(5.+2.*sqrt(10./7.))/3. , -sqrt(5.-2.*sqrt(10./7.))/3., 0. ,
                             sqrt(5.-2.*sqrt(10./7.))/3. ,  sqrt(5.+2.*sqrt(10./7.))/3. };
            break;
    case 6: quad_points = { 0.6612093864662645, -0.6612093864662645, -0.2386191860831969,
                            0.2386191860831969, -0.9324695142031521,  0.9324695142031521 };
            break;
    case 7: quad_points = { 0.0000000000000000,  0.4058451513773972, -0.4058451513773972,
                           -0.7415311855993945,  0.7415311855993945, -0.9491079123427585,
                            0.9491079123427585 };
            break;
    case 8: quad_points = {-0.1834346424956498,  0.1834346424956498, -0.5255324099163290,
                            0.5255324099163290, -0.7966664774136267,  0.7966664774136267,
                           -0.9602898564975363,  0.9602898564975363 };
            break;
    case 9: quad_points = { 0.0000000000000000, -0.8360311073266358,  0.8360311073266358,
                           -0.9681602395076261,  0.9681602395076261, -0.3242534234038089,
                            0.3123470770400029,  0.2606106964029354,  0.2606106964029354 };
            break;
    default: assert( 0 == 1 );
  }
  
  assert (num_of_points == quad_points.size() );
  // Transform quadrature points from [-1,1] -> [0,1]
  for_each(quad_points.begin(), quad_points.end(), [](double& pt){ pt = 0.5 * ( pt + 1. ); });
  return quad_points;
}


vector<double> quadrature_weights(const unsigned int num_of_points)
{
  assert(1 <= num_of_points && num_of_points <= 9);
  vector<double> quad_weights;
  
  switch (num_of_points)
  {
    case 1: quad_weights = { 2. };
            break;
    case 2: quad_weights = { 1. , 1. };
            break;
    case 3: quad_weights = { 5./9. , 8./9. , 5./9. };
            break;
    case 4: quad_weights = { 1./36.*(18. - sqrt(30.)) , 1./36.*(18. + sqrt(30.)) ,
                             1./36.*(18. + sqrt(30.)) , 1./36.*(18. - sqrt(30.)) };
            break;
    case 5: quad_weights = { 1./900.*(322.-13.*sqrt(70.)) , 1./900.*(322.+13.*sqrt(70.)) , 1./900.*(322.+190.) ,
                             1./900.*(322.+13.*sqrt(70.)) , 1./900.*(322.-13.*sqrt(70.)) };
            break;
    case 6: quad_weights = { 0.3607615730481386,  0.3607615730481386,  0.4679139345726910,
                             0.4679139345726910,  0.1713244923791704,  0.1713244923791700 };
            break;
    case 7: quad_weights = { 0.4179591836734694,  0.3818300505051189,  0.3818300505051189,
                             0.2797053914892766,  0.2797053914892766,  0.1294849661688697,
                             0.1294849661688697 };
            break;
    case 8: quad_weights = { 0.3626837833783620,  0.3626837833783620,  0.3137066458778873,
                             0.3137066458778873,  0.2223810344533745,  0.2223810344533745,
                             0.1012285362903763,  0.1012285362903763 };
            break;
    case 9: quad_weights = { 0.3302393550012598,  0.1806481606948574,  0.1806481606948574,
                             0.0812743883615744,  0.0812743883615744,  0.3123470770400029,
                             0.3123470770400029,  0.2606106964029354,  0.2606106964029354 };
            break;
    default: assert( 0 == 1 );
  }
  
  assert (num_of_points == quad_weights.size() );
  // Transform quadrature points from [-1,1] -> [0,1]
  for_each(quad_weights.begin(), quad_weights.end(), [](double& wt){ wt *= 0.5; });
  return quad_weights;
}


vector< vector<double> >  trial_functions_at_quadrature_points
  (const unsigned int max_poly_degree, const unsigned int num_of_points)
{
  assert( 1 <= max_poly_degree && max_poly_degree <= 5 );
  assert( 1 <= num_of_points && num_of_points <= 9);
  
  vector<double> quad_points = quadrature_points(num_of_points);
  vector< vector<double> > fct_val(max_poly_degree+1);
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
  {
    fct_val[i].clear();
    for (unsigned int j = 0; j < num_of_points; ++j)
      fct_val[i].push_back(trial_function_eval(i, quad_points[j]));
  }
  
  return fct_val;
}


vector< vector<double> > derivs_of_trial_at_quadrature_points
  (const unsigned int max_poly_degree, const unsigned int num_of_points)
{
  assert( 1 <= max_poly_degree && max_poly_degree <= 5 );
  assert( 1 <= num_of_points && num_of_points <= 9);
  
  vector<double> quad_points = quadrature_points(num_of_points);
  vector< vector<double> > der_val(max_poly_degree+1);
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
  {
    der_val[i].clear();
    for (unsigned int j = 0; j < num_of_points; ++j)
      der_val[i].push_back(deriv_of_trial_eval(i, quad_points[j]));
  }
  
  return der_val;
}


vector< vector<double> >  trial_functions_at_boundaries(const unsigned int max_poly_degree)
{
  assert( 1 <= max_poly_degree && max_poly_degree <= 5 );
  
  vector< vector<double> > fct_val(max_poly_degree+1);
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
  {
    fct_val[i].clear();
    fct_val[i].push_back(trial_function_eval(i, 0.));
    fct_val[i].push_back(trial_function_eval(i, 1.));
  }
  
  return fct_val;
}


vector< vector<double> > derivs_of_trial_at_boundaries(const unsigned int max_poly_degree)
{
  assert( 1 <= max_poly_degree && max_poly_degree <= 5 );
  
  vector< vector<double> > der_val(max_poly_degree+1);
  
  for (unsigned int i = 0; i < max_poly_degree+1; ++i)
  {
    der_val[i].clear();
    der_val[i].push_back(deriv_of_trial_eval(i, 0.));
    der_val[i].push_back(deriv_of_trial_eval(i, 1.));
  }
  
  return der_val;
}
