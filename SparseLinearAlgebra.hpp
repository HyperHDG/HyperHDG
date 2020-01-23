#ifndef SPARSELINEARALGEBRA_HPP
#define SPARSELINEARALGEBRA_HPP

#include "HyAssert.h"
#include "TypeDefs.h"
#include <vector>
#include <cmath>


dof_value_type inner_product ( const std::vector<dof_value_type>& left,
                               const std::vector<dof_value_type>& right )
{
  hy_assert( left.size() == right.size() ,
             "Both vectors of inner product must be of same size!" );
  
  dof_value_type product = 0.;
  for (dof_index_type i = 0; i < left.size(); ++i)  product += left[i] * right[i];
  return product;
};


dof_value_type norm_2 ( const std::vector<dof_value_type>& vec )
{
  return std::sqrt( inner_product(vec,vec) );
};


std::vector<dof_value_type> linear_combination ( const dof_value_type leftFac,
                                                 const std::vector<dof_value_type>& leftVec,
                                                 const dof_value_type rightFac,
                                                 const std::vector<dof_value_type>& rightVec )
{
  hy_assert( leftVec.size() == rightVec.size() ,
             "Both vectors of linear combination must be of same size!" );
  
  std::vector<dof_value_type> lin_comb ( leftVec.size() , 0. );
  for (dof_index_type i = 0; i < leftVec.size(); ++i)
    lin_comb[i] = leftFac * leftVec[i] + rightFac * rightVec[i];
  return lin_comb;
};


void linear_combination ( const dof_value_type leftFac,  const std::vector<dof_value_type>& leftV,
                          const dof_value_type rightFac, const std::vector<dof_value_type>& rightV,
                          std::vector<dof_value_type>& result )
{
  hy_assert( leftV.size() == rightV.size() && leftV.size() == result.size() ,
             "All three vectors of linear combination must be of same size!" );
  
  for (dof_index_type i = 0; i < result.size(); ++i)
    result[i] = leftFac * leftV[i] + rightFac * rightV[i];
};


template<class ProblemT>
std::vector<dof_value_type> conjugate_gradient
( const std::vector<dof_value_type>& b, const ProblemT& problem,
  int& number_of_iterations = 0, const dof_value_type tolerance = 1e-9 )
{
  std::vector<dof_value_type> x (b.size(), 0.);
  std::vector<dof_value_type> r = b; // Wiki: b - A x (with x = 0)
  std::vector<dof_value_type> d = r;
  
  dof_value_type r_square_old;
  dof_value_type r_square_new = inner_product(r,r);
  
  if (number_of_iterations == 0)  number_of_iterations = b.size();
  
  for (unsigned int k = 0; k < number_of_iterations; ++k)
  {
    std::vector<dof_value_type> z = problem.matrix_vector_multiply(d);
    r_square_old = r_square_new;
    
    dof_value_type alpha = r_square_old / inner_product(d,z);
    linear_combination(1.,x, alpha,d,  x);
    linear_combination(1.,r, -alpha,z, r);
    
    r_square_new = inner_product(r,r);
    
    dof_value_type beta = r_square_new / r_square_old;
    linear_combination(1., r, beta,d,  d);
    
    if ( std::sqrt(r_square_new) < tolerance )  
    {
      number_of_iterations = k;
      return x;
    }
  }

  hy_assert( 0 == 1 ,
             "CG method did not converge! The final residual after " << number_of_iterations <<
             " iterations turned out to be " << std::sqrt(r_square_new) << ", while the needed "
             << "tolerance is " << tolerance << "." );
             
  number_of_iterations = -1;
  return x;
}

#endif // end of ifndef SPARSELINEARALGEBRA_HPP
