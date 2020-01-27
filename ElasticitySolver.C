/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "ElasticitySolver.h"
#include "HyAssert.h"
#include "LapackWrapper.h"
#include <cmath>
#include <vector>


using namespace std;
using namespace FuncQuad;
using namespace Geometry;
#include "ElasticitySolver.inst"


inline vector<double> dyadic_product(const vector<double>& left, const vector<double>& right)
{
  vector<double> result(left.size() * right.size());
  for (unsigned int i = 0; i < right.size(); ++i)
    for (unsigned int j = 0; j < left.size(); ++j)
      result[i * left.size() + j] = left[j] * right[i];
  return result;
}


inline vector< vector<double> > double_dyadic_product(const vector< vector<double> >& left, const vector< vector<double> >& right)
{
  vector< vector<double> > result(left.size() * right.size());
  for (unsigned int i = 0; i < right.size(); ++i)
    for (unsigned int j = 0; j < left.size(); ++j)
      result[i * left.size() + j] = dyadic_product(left[j], right[i]);
  return result;
}


inline vector<double> get_relevant_coeffs_indicator(const unsigned int hyperedge_dim, const unsigned int max_poly_degree, const unsigned int dimension, const unsigned int ansatz)
{
  vector<double> unity_vec(max_poly_degree+1, 1.);
  vector<double> ansatz_vec(max_poly_degree+1, 0.);
  hy_assert( ansatz < ansatz_vec.size() , "Ansatz function index in one dimension must be smaller than maximal degree + 1." );
  ansatz_vec[ansatz] = 1;
  vector<double> result;
  if( hyperedge_dim == 1 )  result = unity_vec;
  else                      result = ansatz_vec;
  
  for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
    if (dim == 1 && dimension == 0)  result = dyadic_product(unity_vec, result);
    else if (dim == dimension)       result = dyadic_product(result, unity_vec);
    else                             result = dyadic_product(result, ansatz_vec);
  
  return result;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
ElasticitySolver_RegularQuad(const constructor_value_type& tau)
: tau_(tau)
{ 
  static_assert( hyperedge_dim == 1 , "This has only been implemented for one dimensional hyperedges." );
//  hy_assert( 0 == 1 , "Not yet implemented!" );
  array<double, compute_n_quad_points(max_quad_degree)>
    quad_weights1D = quadrature_weights<max_quad_degree>();
  array< array<double, compute_n_quad_points(max_quad_degree)> , max_poly_degree + 1 >
    trials_at_quad1D = trial_functions_at_quadrature_points<max_poly_degree, max_quad_degree>();
  array< array<double, compute_n_quad_points(max_quad_degree)> , max_poly_degree + 1 >
    derivs_at_quad1D = derivs_of_trial_at_quadrature_points<max_poly_degree, max_quad_degree>();
  array< array<double, 2> , max_poly_degree + 1 >
    trials_at_bdr1D = trial_functions_at_boundaries<max_poly_degree>();
//  array< array<double, 2> , max_poly_degree + 1 >
//    derivs_at_bdr1D = derivs_of_trial_at_boundaries<max_poly_degree>();
    
  // In the one-dimensional case, we are done now.
  if constexpr (hyperedge_dim == 1)
  {
    quad_weights_ = quad_weights1D;
    quad_bdr_ = { 1. };
    trials_quad_ = trials_at_quad1D;
    bound_trials_quad_[0] = { 1. };
    derivs_quad_[0] = derivs_at_quad1D;
    trials_in_corners_ = trials_at_bdr1D;
    for (unsigned int side = 0; side < 2; ++side)
      for (unsigned int i = 0; i < max_poly_degree + 1; ++i)
        trials_bound_[side][i][0] = trials_at_bdr1D[i][side];
  }
  else // There is more to be done and we need to use growing containers.
  {
    // Deal with quadrature points that remain a vector.
    vector<double> quad_weights1D_vec(quad_weights1D.begin(), quad_weights1D.end());
    vector<double> quad_weights_vec(quad_weights1D.begin(), quad_weights1D.end());
    if (quad_bdr_.size() == quad_weights1D_vec.size())
      for (unsigned int i = 0; i < quad_bdr_.size(); ++i)  quad_bdr_[i] = quad_weights1D_vec[i];
    for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
    {
      if (dim == hyperedge_dim - 2)
      {
        hy_assert( quad_bdr_.size() == quad_weights_vec.size() ,
                   "Amount of quadrature points should be equal to amount of quadrature weights." );
        for (unsigned int i = 0; i < quad_bdr_.size(); ++i)  quad_bdr_[i] = quad_weights_vec[i];
      }
      quad_weights_vec = dyadic_product(quad_weights_vec, quad_weights1D_vec);
    }
    hy_assert( quad_weights_.size() == quad_weights_vec.size() ,
               "Size of array and vector that will become array must be equal." );
    for (unsigned int i = 0; i < quad_weights_.size(); ++i)  quad_weights_[i] = quad_weights_vec[i];
    
    // Deal with trials in corners that is a vector (func) of vectors (corners).
    vector< vector<double> > trials_in_corners1D_vec(trials_at_bdr1D.size());
    for (unsigned int i = 0; i < trials_at_bdr1D.size(); ++i)
    {
      trials_in_corners1D_vec[i].resize(trials_at_bdr1D[i].size());
      for (unsigned int j = 0; j < trials_at_quad1D[i].size(); ++j)
        trials_in_corners1D_vec[i][j] = trials_at_bdr1D[i][j];
    }
    vector< vector<double> > trials_in_corners_vec = trials_in_corners1D_vec;
    for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
      trials_in_corners_vec = double_dyadic_product(trials_in_corners_vec, trials_in_corners1D_vec);
    hy_assert( trials_in_corners_.size() == trials_in_corners_vec.size() ,
               "Size of array and vector that will become array must be equal." );
    for (unsigned int i = 0; i < trials_in_corners_.size(); ++i)
    {
      hy_assert( trials_in_corners_[i].size() == trials_in_corners_vec[i].size() ,
                 "Size of array and vector that will become array must be equal." );
      for (unsigned int j = 0; j < trials_quad_[i].size(); ++j)
        trials_in_corners_[i][j] = trials_in_corners_vec[i][j];
    }
    
    // Deal with trials at quadrature points that remain a vector (func) of vectors (quad).
    vector< vector<double> > trials_at_quad1D_vec(trials_at_quad1D.size());
    for (unsigned int i = 0; i < trials_at_quad1D.size(); ++i)
    {
      trials_at_quad1D_vec[i].resize(trials_at_quad1D[i].size());
      for (unsigned int j = 0; j < trials_at_quad1D[i].size(); ++j)
        trials_at_quad1D_vec[i][j] = trials_at_quad1D[i][j];
    }
    vector< vector<double> > trials_at_quad_vec = trials_at_quad1D_vec;
    for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
      trials_at_quad_vec = double_dyadic_product(trials_at_quad_vec, trials_at_quad1D_vec);
    hy_assert( trials_quad_.size() == trials_at_quad_vec.size() ,
               "Size of array and vector that will become array must be equal." );
    for (unsigned int i = 0; i < trials_quad_.size(); ++i)
    {
      hy_assert( trials_quad_[i].size() == trials_at_quad_vec[i].size() ,
                 "Size of array and vector that will become array must be equal." );
      for (unsigned int j = 0; j < trials_quad_[i].size(); ++j)
        trials_quad_[i][j] = trials_at_quad_vec[i][j];
    }
    
    // Deal with boundary trials at quadrature points that remain a vector (func) of vectors (quad).
    vector< vector<double> > bound_trials_at_quad_vec = trials_at_quad1D_vec;
    for (unsigned int dim = 1; dim < hyperedge_dim - 1; ++dim)
      bound_trials_at_quad_vec = double_dyadic_product(bound_trials_at_quad_vec, trials_at_quad1D_vec);
    hy_assert( bound_trials_quad_.size() == bound_trials_at_quad_vec.size() ,
               "Size of array and vector that will become array must be equal." );
    for (unsigned int i = 0; i < bound_trials_quad_.size(); ++i)
    {
      hy_assert( bound_trials_quad_[i].size() == bound_trials_at_quad_vec[i].size() ,
                 "Size of array and vector that will become array must be equal." );
      for (unsigned int j = 0; j < bound_trials_quad_[i].size(); ++j)
        bound_trials_quad_[i][j] = bound_trials_at_quad_vec[i][j];
    }
    
    // Deal with derivatives with respect to first_index of function of second_index at quadrature point
    // third_index. In 1D only first_index = 0 makes sense.
    vector< vector<double> > derivs_at_quad1D_vec(derivs_at_quad1D.size());
    for (unsigned int i = 0; i < derivs_at_quad1D.size(); ++i)
    {
      derivs_at_quad1D_vec[i].resize(derivs_at_quad1D[i].size());
      for (unsigned int j = 0; j < derivs_at_quad1D[i].size(); ++j)
        derivs_at_quad1D_vec[i][j] = derivs_at_quad1D[i][j];
    }
    vector< vector< vector<double> > > derivs_at_quad_vec(hyperedge_dim);
    for (unsigned int dim_deriv = 0; dim_deriv < hyperedge_dim; ++dim_deriv)
    {
      derivs_at_quad_vec[dim_deriv] = trials_at_quad1D_vec;
      for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
        if(dim_deriv == 0 && dim == 1)  derivs_at_quad_vec[dim_deriv] = double_dyadic_product(derivs_at_quad1D_vec, trials_at_quad1D_vec);
        else if (dim_deriv == dim)      derivs_at_quad_vec[dim_deriv] = double_dyadic_product(derivs_at_quad_vec[dim_deriv], derivs_at_quad1D_vec);
        else                            derivs_at_quad_vec[dim_deriv] = double_dyadic_product(derivs_at_quad_vec[dim_deriv], trials_at_quad1D_vec);
    }
    hy_assert( derivs_quad_.size() == derivs_at_quad_vec.size() ,
               "Size of array and vector that will become array must be equal." );
    for (unsigned int i = 0; i < derivs_quad_.size(); ++i)
    {
      hy_assert( derivs_quad_[i].size() == derivs_at_quad_vec[i].size() ,
                 "Size of array and vector that will become array must be equal." );
      for (unsigned int j = 0; j < derivs_quad_[i].size(); ++j)
      {
        hy_assert( derivs_quad_[i][j].size() == derivs_at_quad_vec[i][j].size() ,
                   "Size of array and vector that will become array must be equal." );
        for (unsigned int k = 0; k < derivs_quad_[i][j].size(); ++k)
          derivs_quad_[i][j][k] = derivs_at_quad_vec[i][j][k];
      }
    }
    
    // Deal with trials at boundary which is a vector (boundary) of vectors (function) of vectors (quad).
    vector< vector< vector<double> > > trials_bound_vec(2 * hyperedge_dim);
    for (unsigned int side = 0; side < 2; ++side)
    {
      vector< vector<double> > helperling(trials_at_bdr1D.size());
      for (unsigned int i = 0; i < helperling.size(); ++i)
        helperling[i] = vector<double>(1, trials_at_bdr1D[i][side]);
      for (unsigned int dim_face = 0; dim_face < hyperedge_dim; ++dim_face)
      {
        trials_bound_vec[2*dim_face+side] = trials_at_quad1D_vec;
        for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
          if (dim_face == 0 && dim == 1)  trials_bound_vec[2*dim_face+side] = double_dyadic_product(helperling, trials_at_quad1D_vec);
          else if (dim_face == dim)       trials_bound_vec[2*dim_face+side] = double_dyadic_product(trials_bound_vec[2*dim_face+side], helperling);
          else                            trials_bound_vec[2*dim_face+side] = double_dyadic_product(trials_bound_vec[2*dim_face+side], trials_at_quad1D_vec);
      }
    }
    hy_assert( trials_bound_.size() == trials_bound_vec.size() ,
               "Size of array and vector that will become array must be equal." );
    for (unsigned int i = 0; i < trials_bound_.size(); ++i)
    {
      hy_assert( trials_bound_[i].size() == trials_bound_vec[i].size() ,
                 "Size of array and vector that will become array must be equal." );
      for (unsigned int j = 0; j < trials_bound_[i].size(); ++j)
      {
        hy_assert( trials_bound_[i][j].size() == trials_bound_vec[i][j].size() ,
                   "Size of array and vector that will become array must be equal." );
        for (unsigned int k = 0; k < trials_bound_[i][j].size(); ++k)
          trials_bound_[i][j][k] = trials_bound_vec[i][j][k];
      }
    }
  }
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
inline unsigned int ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
loc_matrix_index(const unsigned int row, const unsigned int column) const
{
  hy_assert( 0 <= row ,
             "Row index should be larger than or equal to zero." );
  hy_assert( row < (hyperedge_dim + 1) * pow((max_poly_degree + 1), hyperedge_dim) ,
             "Row index should be smaller than total amount of rows." );
  hy_assert( 0 <= column ,
             "Column index should be larger than or equal to zero." );
  hy_assert( column < (hyperedge_dim + 1) * pow((max_poly_degree + 1), hyperedge_dim) ,
             "Column index should smaller than total amount of columns." );
  return column * (hyperedge_dim + 1) * pow((max_poly_degree + 1), hyperedge_dim) + row;  // Transposed for LAPACK
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
inline auto // array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_>
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
assemble_loc_mat() const
{
  double hyperedge_area = 1.;
  array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_> local_mat;
  local_mat.fill(0.);
  
  for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
    for (unsigned int i = 0; i < num_ansatz_fct_; ++i)
    {
      for (unsigned int j = 0; j < num_ansatz_fct_; ++j)
        for (unsigned int q = 0; q < num_of_quad_; ++q)
          local_mat[loc_matrix_index( dim * num_ansatz_fct_ + i , dim * num_ansatz_fct_ + j )] += 
            quad_weights_[q] * hyperedge_area * trials_quad_[i][q] * trials_quad_[j][q];
      for (unsigned int j = 0; j < num_ansatz_fct_; ++j)
        for (unsigned int q = 0; q < num_of_quad_; ++q)
          local_mat[loc_matrix_index(  dim * num_ansatz_fct_ + i , hyperedge_dim * num_ansatz_fct_ + j )] -=
            quad_weights_[q] * derivs_quad_[dim][i][q] * trials_quad_[j][q];
    }
  
  for (unsigned int i = 0; i < num_ansatz_fct_; ++i)
  {
    for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
    {
      for (unsigned int j = 0; j < num_ansatz_fct_; ++j)
      {
        for (unsigned int q = 0; q < num_of_quad_; ++q)
          local_mat[loc_matrix_index( hyperedge_dim * num_ansatz_fct_ + i , dim * num_ansatz_fct_ + j )] -=
            quad_weights_[q] * derivs_quad_[dim][i][q] * trials_quad_[j][q];
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          local_mat[loc_matrix_index( hyperedge_dim * num_ansatz_fct_ + i , dim * num_ansatz_fct_ + j )] +=
            + quad_bdr_[q] * trials_bound_[2*dim+1][i][q] * trials_bound_[2*dim+1][j][q]
            - quad_bdr_[q] * trials_bound_[2*dim+0][i][q] * trials_bound_[2*dim+0][j][q];
      }
      for (unsigned int j = 0; j < num_ansatz_fct_; ++j)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          local_mat[loc_matrix_index( hyperedge_dim * num_ansatz_fct_ + i , hyperedge_dim * num_ansatz_fct_ + j )] +=
            tau_ * quad_bdr_[q] * (trials_bound_[2*dim+0][i][q] * trials_bound_[2*dim+0][j][q]
                                   + trials_bound_[2*dim+1][i][q] * trials_bound_[2*dim+1][j][q]);
    }
  }
/*  
  cout << endl << endl;
  for (unsigned int i = 0; i < (hyperedge_dim + 1) * pow((max_poly_degree + 1), hyperedge_dim); ++i)
  {
    for(unsigned int j = 0; j < (hyperedge_dim + 1) * pow((max_poly_degree + 1), hyperedge_dim); ++j)
      cout << local_mat[loc_matrix_index(i,j)] << "  ";
    cout << endl;
  }
*/  
  return local_mat;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
inline auto // array<double, (hyperedge_dim+1) * num_ansatz_fct_>
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
assemble_rhs(const array< array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const
{
  array<double, (hyperedge_dim+1) * num_ansatz_fct_> right_hand_side;
  right_hand_side.fill(0.);
  
  hy_assert( lambda_values.size() == 2 * hyperedge_dim ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyperedge_dim; ++i)
    hy_assert( lambda_values[i].size() == num_ansatz_bdr_ ,
               "The size of the lambda values should be the amount of ansatz functions ar boundary." );
  
  for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
    for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
      for (unsigned int bdr_ans = 0; bdr_ans < num_ansatz_bdr_; ++bdr_ans)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          right_hand_side[dim * num_ansatz_fct_ + ansatz] += 
            + quad_bdr_[q] * lambda_values[2*dim+0][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+0][ansatz][q]
            - quad_bdr_[q] * lambda_values[2*dim+1][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+1][ansatz][q];
  
  for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
    for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
      for (unsigned int bdr_ans = 0; bdr_ans < num_ansatz_bdr_; ++bdr_ans)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          right_hand_side[hyperedge_dim * num_ansatz_fct_ + ansatz] += tau_ * quad_bdr_[q] *
            ( lambda_values[2*dim+0][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+0][ansatz][q]
            + lambda_values[2*dim+1][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+1][ansatz][q] );
/*  
  cout << endl;
  for (unsigned int i = 0; i < right_hand_side.size(); ++i)  cout << "  " << right_hand_side[i];
*/  
  return right_hand_side;
}


/*
template<unsigned int dim, unsigned int unknown_dim>
vector<double> DiffusionSolver<dim,unknown_dim>::solve_local_system_of_eq(const vector<double>& loc_matrix, const vector<double>& loc_rhs) const
{
  assert( loc_matrix.size() == loc_rhs.size() * loc_rhs.size() );
  const int system_size = loc_rhs.size();
  vector<double> solution(system_size, 0.);
  
  Eigen::MatrixXd eigen_matrix(system_size, system_size);
  Eigen::VectorXd eigen_rhs(system_size), eigen_sol(system_size);
  
  for(int i = 0; i < system_size; ++i)
    for(int j = 0; j < system_size; ++j)
      eigen_matrix(i,j) = loc_matrix[loc_matrix_index(i,j)];
      
  for(int i = 0; i < system_size; ++i)  eigen_rhs(i) = loc_rhs[i];
  
  eigen_sol = eigen_matrix.colPivHouseholderQr().solve(eigen_rhs);
  
  for(int i = 0; i < system_size; ++i)  solution[i] = eigen_sol(i);
  
  return solution;
}
*/


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array<double, (hyperedge_dim+1) * num_ansatz_fct_>
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
solve_local_system_of_eq(array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_>& loc_matrix,
                         array<double, (hyperedge_dim+1) * num_ansatz_fct_>& loc_rhs) const
{
  hy_assert( loc_matrix.size() == loc_rhs.size() * loc_rhs.size() ,
             "The size of a local matrix should be the size of the right-hand side squared." );
  const int system_size = loc_rhs.size();
  double *mat_a=loc_matrix.data(), *rhs_b = loc_rhs.data();
  int info = -1;  
  lapack_solve(system_size, mat_a, rhs_b, &info);
  hy_assert( info == 0 ,
             "LAPACK's solve failed and the solution of the local problem might be inaccurate." );
  return loc_rhs;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
inline auto // array<double, (hyperedge_dim+1) * num_ansatz_fct_>
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
solve_local_problem(const array< array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const
{
  array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_> local_matrix;
  array<double, (hyperedge_dim+1) * num_ansatz_fct_> right_hand_side;
  right_hand_side = assemble_rhs(lambda_values);
  local_matrix = assemble_loc_mat();
  
  array<double, (hyperedge_dim+1) * num_ansatz_fct_> solution = solve_local_system_of_eq(local_matrix, right_hand_side);
/*  
  cout << endl;
  for (unsigned int i = 0; i < solution.size(); ++i) cout << "  " << solution[i];
  cout << endl;
*/  
  return solution;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
inline auto // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
dual_at_boundary(const array<double, (hyperedge_dim+1) * num_ansatz_fct_>& coeffs) const
{
  array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim > bdr_values;
  
  for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
  {
    bdr_values[2*dim].fill(0.);
    bdr_values[2*dim+1].fill(0.);
    for(unsigned int bdr_ansatz = 0; bdr_ansatz < num_ansatz_bdr_; ++bdr_ansatz)
      for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
        {
          bdr_values[2*dim+0][bdr_ansatz] -= quad_bdr_[q] * coeffs[dim * num_ansatz_fct_ + ansatz] * trials_bound_[2*dim+0][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
          bdr_values[2*dim+1][bdr_ansatz] += quad_bdr_[q] * coeffs[dim * num_ansatz_fct_ + ansatz] * trials_bound_[2*dim+1][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
        }
  }
  
  return bdr_values;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
inline auto // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim > 
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
primal_at_boundary(const array<double, (hyperedge_dim+1) * num_ansatz_fct_>& coeffs) const
{
  array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim > bdr_values;
  
  for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
  {
    bdr_values[2*dim].fill(0.);
    bdr_values[2*dim+1].fill(0.);
    for (unsigned int bdr_ansatz = 0; bdr_ansatz < num_ansatz_bdr_; ++bdr_ansatz)
      for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
        {
          bdr_values[2*dim+0][bdr_ansatz] += quad_bdr_[q] * coeffs[hyperedge_dim * num_ansatz_fct_ + ansatz] * trials_bound_[2*dim+0][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
          bdr_values[2*dim+1][bdr_ansatz] += quad_bdr_[q] * coeffs[hyperedge_dim * num_ansatz_fct_ + ansatz] * trials_bound_[2*dim+1][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
        }
  }
  
  return bdr_values;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
numerical_flux_at_boundary(const array< array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values, const array<double, (hyperedge_dim+1) * num_ansatz_fct_>& coeffs) const
{
  array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim > bdr_values;
  array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim > primals = primal_at_boundary(coeffs);
  array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim > duals   = dual_at_boundary(coeffs);
  
  for (unsigned int i = 0; i < lambda_values.size(); ++i)
    for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
      bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];
/*  
  for(unsigned int i = 0; i < bdr_values.size(); ++i)
  {
    for(unsigned int j = 0; j < bdr_values[i].size(); ++j)
      cout << bdr_values[i][j] << "  ";
    cout << endl;
  }
*/  
  return bdr_values;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
array<double, compute_n_corners_of_cube(hyperedge_dim)>
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
primal_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const
{
  array<double, (hyperedge_dim+1) * num_ansatz_fct_> coefficients = solve_local_problem(lambda_values);
  array<double, compute_n_corners_of_cube(hyperedge_dim)> primal_in_corners;
  primal_in_corners.fill(0.);
  for (unsigned int corner = 0; corner < compute_n_corners_of_cube(hyperedge_dim); ++corner)
    for (unsigned int ansatz_fct = 0; ansatz_fct < num_ansatz_fct_; ++ansatz_fct)
      primal_in_corners[corner] += coefficients[hyperedge_dim * num_ansatz_fct_ + ansatz_fct] * trials_in_corners_[ansatz_fct][corner];
  return primal_in_corners;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
array< array<double, hyperedge_dim> , compute_n_corners_of_cube(hyperedge_dim) >
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
dual_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const
{
  array<double, (hyperedge_dim+1) * num_ansatz_fct_> coefficients = solve_local_problem(lambda_values);
  array< array<double, hyperedge_dim> , compute_n_corners_of_cube(hyperedge_dim) > dual_in_corners;
  for (unsigned int corner = 0; corner < compute_n_corners_of_cube(hyperedge_dim); ++corner)
  {
    dual_in_corners[corner].fill(0.);
    for (unsigned int ansatz_fct = 0; ansatz_fct < num_ansatz_fct_; ++ansatz_fct)
      for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
        dual_in_corners[corner][dim] += coefficients[dim * num_ansatz_fct_ + ansatz_fct] * trials_in_corners_[ansatz_fct][corner];
  }
  return dual_in_corners;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
array< array<double, compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim > // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
numerical_flux_from_lambda(const array< array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const
{
  return numerical_flux_at_boundary(lambda_values, solve_local_problem(lambda_values));
}

template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
array< array<double, compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim > // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
preprocess_data( array< array<double, space_dim * compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2*hyperedge_dim >& hyperedge_dofs,
                 HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim>& geometry ) const
{
  array< array<double, compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2*hyperedge_dim > result;
  hy_assert( result.size() == 2 , "Only implemented in one dimension!" );
  for (unsigned int i = 0; i < result.size(); ++i){
    hy_assert( result[i].size() == 1 , "Only implemented in one dimension!" );
    result[i].fill(0.);
  }
  
  for (unsigned int i = 0; i < 2 * hyperedge_dim; ++i)
  {
    Point<space_dim> normal_vector = geometry.normal(1);
    for (unsigned int dim = 0; dim < space_dim; ++dim)  result[i][0] += normal_vector[dim] * hyperedge_dofs[i][dim];
  }
  
  return result;
}


template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
array< array<double, space_dim * compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2*hyperedge_dim >
ElasticitySolver_RegularQuad<hyperedge_dim, space_dim, max_poly_degree, max_quad_degree>::
postprocess_data( array< array<double, compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2*hyperedge_dim >& hyperedge_dofs,
                  HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim>& geometry ) const
{
  std::array< std::array<double, space_dim * compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2*hyperedge_dim > result;
  for (unsigned int i = 0; i < result.size(); ++i)  result[i].fill(0.);
  
  for (unsigned int i = 0; i < 2 * hyperedge_dim; ++i)
  {
    Point<space_dim> normal_vector = geometry.normal(1);
    for (unsigned int dim = 0; dim < space_dim; ++dim)  result[i][dim] += normal_vector[dim] * hyperedge_dofs[i][0];
  }
  
  return result;
}

