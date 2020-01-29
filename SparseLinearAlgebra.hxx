/*!*************************************************************************************************
 * \file    SparseLinearAlgebra.hxx
 * \brief   A file containing different functions that implement basic linear algebra operations
 *          using large \c std::vector.
 *
 * This namespace provides several functions to implement basic linear algebra operations of (long)
 * \c std::vector in combination with a class providing a function \c matrix_vector_multiply. This
 * is mainly used for C++ examples and test cases that do not use the Python interface and its
 * version of an CG method, for example.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#ifndef SPARSELINEARALGEBRA_HXX
#define SPARSELINEARALGEBRA_HXX

#include "HyAssert.hxx"
#include "TypeDefs.hxx"
#include <vector>
#include <cmath>

/*!*************************************************************************************************
 * \brief   A namespace containing different functions that implement basic linear algebra
 *          operations using large \c std::vector.
 * 
 * \todo    Check, whether this construction makes sense.
 *
 * This namespace provides several functions to implement basic linear algebra operations of (long)
 * \c std::vector in combination with a class providing a function \c matrix_vector_multiply. This
 * is mainly used for C++ examples and test cases that do not use the Python interface and its
 * version of an CG method, for example.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
namespace SparseLA
{

/*!*************************************************************************************************
 * \brief   Evaluate the inner product of two \c std::vector.
 * 
 * Naive implementation of an Euclidean inner product of two \c std::vector which are supposed to be
 * of the same size. This function is needed to calculate a vector's 2 norm or to implement a CG
 * scheme.
 * 
 * \param   left                Left argument of the inner product.
 * \param   right               Right argument of the inner product.
 * \retval  product             Inner product of the two arguments.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
dof_value_t inner_product ( const std::vector<dof_value_t>& left,
                               const std::vector<dof_value_t>& right )
{
  hy_assert( left.size() == right.size() ,
             "Both vectors of inner product must be of same size!" );
  
  dof_value_t product = 0.;
  for (dof_index_type i = 0; i < left.size(); ++i)  product += left[i] * right[i];
  return product;
};

/*!*************************************************************************************************
 * \brief   Evaluate 2 norm of a \c std::vector.
 * 
 * Naive implementation of an 2 norm of a vector. This is the square root of the \c inner_product of
 * a vector paired with itself.
 * 
 * \param   vec                 Vector whose 2 norm is to be calculates.
 * \retval  norm                2 norm of given vector.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
dof_value_t norm_2 ( const std::vector<dof_value_t>& vec )
{
  return std::sqrt( inner_product(vec,vec) );
};

/*!*************************************************************************************************
 * \brief   Evaluate linear combination of vectors and return the result.
 * 
 * This functions takes two \c std::vector and two \c dof_value_t and returns their linear
 * combination "leftFac * leftVec + rightFac * rightVec" as a new vector (in contrast to just a
 * reference to a vector).
 * 
 * \param   leftFac             Scaling factor of left vector.
 * \param   leftVec             Left vector in linear combination.
 * \param   rightFac            Scaling factor of right vector.
 * \param   rightVec            Right vector in linear combination.
 * \retval  lin_comb            Linear combination of vectors with respective coefficients.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
std::vector<dof_value_t> linear_combination ( const dof_value_t leftFac,
                                                 const std::vector<dof_value_t>& leftVec,
                                                 const dof_value_t rightFac,
                                                 const std::vector<dof_value_t>& rightVec )
{
  hy_assert( leftVec.size() == rightVec.size() ,
             "Both vectors of linear combination must be of same size!" );
  
  std::vector<dof_value_t> lin_comb ( leftVec.size() , 0. );
  for (dof_index_type i = 0; i < leftVec.size(); ++i)
    lin_comb[i] = leftFac * leftVec[i] + rightFac * rightVec[i];
  return lin_comb;
};

/*!*************************************************************************************************
 * \brief   Evaluate linear combination of vectors and return reference to result.
 * 
 * This functions takes two \c std::vector and two \c dof_value_t and returns their linear
 * combination "leftFac * leftVec + rightFac * rightVec" as a reference to a vector. This vector
 * needs to be passed to the function
 * 
 * \param   leftFac             Scaling factor of left vector.
 * \param   leftVec             Left vector in linear combination.
 * \param   rightFac            Scaling factor of right vector.
 * \param   rightVec            Right vector in linear combination.
 * \param   result              Reference to vector whicb is supposed to contain the result.
 * \retval  result              Linear combination of vectors with respective coefficients.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
void linear_combination ( const dof_value_t leftFac,  const std::vector<dof_value_t>& leftV,
                          const dof_value_t rightFac, const std::vector<dof_value_t>& rightV,
                          std::vector<dof_value_t>& result )
{
  hy_assert( leftV.size() == rightV.size() && leftV.size() == result.size() ,
             "All three vectors of linear combination must be of same size!" );
  
  for (dof_index_type i = 0; i < result.size(); ++i)
    result[i] = leftFac * leftV[i] + rightFac * rightV[i];
};

/*!*************************************************************************************************
 * \brief   Execute conjugate gradient algorithm to find solution to system of equations.
 * 
 * Execute conjugate gradient algorithm where the matrix is not explicitly given, but the template
 * class \c ProblemT is supposed to implement a function \c matrix_vector_multiply which only takes
 * a \c std::vector and generates the matrix vector product from that. The associated matrix is
 * assumed to be square and symmetric positive definite.
 * 
 * \tparam  ProblemT            Class to implement matrix vector multiplication.
 * \param   b                   Right-hand side of linear system of equations.
 * \param   problme             Class instantiation to implement matrix vector multiplication.
 * \param   n_iterations        Maximum number of iterations. 0 is default and the size of b.
 * \param   tolerance           Absolute tolerance value in 2 norm. Default is 1e-9.
 * \retval  solution            Vector sufficing Ax = b up to given tolerance if converged.
 * \retval  n_iterations        Number of needed iterations. -1 indicates no convergence.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<class ProblemT>
std::vector<dof_value_t> conjugate_gradient
( const std::vector<dof_value_t>& b, const ProblemT& problem,
  int& n_iterations = 0, const dof_value_t tolerance = 1e-9 )
{
  std::vector<dof_value_t> x (b.size(), 0.);
  std::vector<dof_value_t> r = b; // Wiki: b - A x (with x = 0)
  std::vector<dof_value_t> d = r;
  
  dof_value_t r_square_old;
  dof_value_t r_square_new = inner_product(r,r);
  
  if (n_iterations == 0)  n_iterations = b.size();
  
  for (unsigned int k = 0; k < n_iterations; ++k)
  {
    std::vector<dof_value_t> z = problem.matrix_vector_multiply(d);
    r_square_old = r_square_new;
    
    dof_value_t alpha = r_square_old / inner_product(d,z);
    linear_combination(1.,x, alpha,d,  x);
    linear_combination(1.,r, -alpha,z, r);
    
    r_square_new = inner_product(r,r);
    
    dof_value_t beta = r_square_new / r_square_old;
    linear_combination(1., r, beta,d,  d);
    
    if ( std::sqrt(r_square_new) < tolerance )  
    {
      n_iterations = k;
      return x;
    }
  }

  hy_assert( 0 == 1 ,
             "CG method did not converge! The final residual after " << n_iterations <<
             " iterations turned out to be " << std::sqrt(r_square_new) << ", while the needed "
             << "tolerance is " << tolerance << "." );
             
  n_iterations = -1;
  return x;
};

} // end of namespace SparseLA

#endif // end of ifndef SPARSELINEARALGEBRA_HXX
