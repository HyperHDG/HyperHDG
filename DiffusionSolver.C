/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "DiffusionSolver.h"
// #include "FuncAndQuad.h"
#include "LapackWrapper.h"
#include <cassert>
#include <cmath>
// #include <Eigen/Dense>


#include <iostream>


using namespace std;


template class DiffusionSolver_RegularQuad<1,1,1>;
template class DiffusionSolver_RegularQuad<1,1,2>;
template class DiffusionSolver_RegularQuad<1,1,3>;
template class DiffusionSolver_RegularQuad<1,2,2>;
template class DiffusionSolver_RegularQuad<1,2,4>;
template class DiffusionSolver_RegularQuad<1,2,6>;
template class DiffusionSolver_RegularQuad<1,3,3>;
template class DiffusionSolver_RegularQuad<1,3,6>;
template class DiffusionSolver_RegularQuad<1,3,9>;
template class DiffusionSolver_RegularQuad<1,4,4>;
template class DiffusionSolver_RegularQuad<1,4,8>;
template class DiffusionSolver_RegularQuad<1,4,12>;
template class DiffusionSolver_RegularQuad<1,5,5>;
template class DiffusionSolver_RegularQuad<1,5,10>;
template class DiffusionSolver_RegularQuad<1,5,15>;




//template class DiffusionSolver_RegularQuad<2>;
//template class DiffusionSolver_RegularQuad<3>;


vector<double> dyadic_product(const vector<double>& left, const vector<double>& right)
{
  vector<double> result(left.size() * right.size());
  for (unsigned int i = 0; i < right.size(); ++i)
    for (unsigned int j = 0; j < left.size(); ++j)
      result[i * left.size() + j] = left[j] * right[i];
  return result;
}


template<unsigned int left_size, unsigned int right_size>
array<double, left_size * right_size> dyadic_product(const array<double, left_size> left, const array<double, right_size> right)
{
  array<double, left_size * right_size> result;
  for (unsigned int i = 0; i < right_size; ++i)
    for (unsigned int j = 0; j < left_size; ++j)
      result[i * left_size + j] = left[j] * right[i];
  return result;
}


vector< vector<double> > double_dyadic_product(const vector< vector<double> >& left, const vector< vector<double> >& right)
{
  vector< vector<double> > result(left.size() * right.size());
  for (unsigned int i = 0; i < right.size(); ++i)
    for (unsigned int j = 0; j < left.size(); ++j)
      result[i * left.size() + j] = dyadic_product(left[j], right[i]);
  return result;
}


template<unsigned int left_outer, unsigned int left_inner, unsigned int right_outer, unsigned int right_inner>
array< array<double, left_inner * right_inner> , left_outer * right_outer >
double_dyadic_product(const array< array<double, left_inner> , left_outer > left, const array< array<double, right_inner> , right_outer > right)
{
  array< array<double, left_inner * right_inner> , left_outer * right_outer > result;
  for(unsigned int i = 0; i < right_outer; ++i)
    for(unsigned int j = 0; j < left_outer; ++j)
      result[i * left_outer + j] = dyadic_product(left[j], right[i]);
  return result;
}


vector<double> get_relevant_coeffs_indicator(const unsigned int hyperedge_dim, const unsigned int max_poly_degree, const unsigned int dimension, const unsigned int ansatz)
{
  vector<double> unity_vec(max_poly_degree+1, 1.);
  vector<double> ansatz_vec(max_poly_degree+1, 0.);
  assert( ansatz < ansatz_vec.size() );
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

/*
template<unsigned int hyperedge_dim, unsigned int max_poly_degree>
vector<double> get_relevant_coeffs_indicator(const unsigned int dimension, const unsigned int ansatz)
{
  array<double, max_poly_degree + 1> unity_vec, ansatz_vec;
  unity_vec.fill(1.);
  ansatz_vec.fill(0.);
  assert( ansatz < ansatz_vec.size() );
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
*/

template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
DiffusionSolver_RegularQuad(const double tau)
: tau_(tau)
{ 
//  constexpr unsigned int num_1D_quad = quadrature_points_amount(max_quad_degree);
  
  array<double, quadrature_points_amount(max_quad_degree)>
    quad_weights1D = quadrature_weights<max_quad_degree>();
  array<double, 1>
    quad_weights_bdr1D = { 1. };
  array< array<double, quadrature_points_amount(max_quad_degree)> , max_poly_degree + 1 >
    trials_at_quad1D = trial_functions_at_quadrature_points<max_poly_degree, max_quad_degree>();
  array< array<double, quadrature_points_amount(max_quad_degree)> , max_poly_degree + 1 >
    derivs_at_quad1D = derivs_of_trial_at_quadrature_points<max_poly_degree, max_quad_degree>();
  array< array<double, 2> , max_poly_degree + 1 >
    trials_at_bdr1D = trial_functions_at_boundaries<max_poly_degree>();
//  array< array<double, 2> , max_poly_degree + 1 >
//    derivs_at_bdr1D = derivs_of_trial_at_boundaries<max_poly_degree>();
  cout << num_ansatz_fct_<< " " << trials_bound_.size() << " " << trials_bound_[0].size() <<  endl;
  if constexpr (hyperedge_dim == 1)
  {
    quad_weights_ = quad_weights1D;
    quad_bdr_ = quad_weights_bdr1D;
    trials_quad_ = trials_at_quad1D;
    bound_trials_quad_[0] = { 1. };
    trials_bound_1D_ = trials_at_bdr1D;
    derivs_quad_[0] = derivs_at_quad1D;
    for (unsigned int side = 0; side < 2; ++side)
      for (unsigned int i = 0; i < max_poly_degree + 1; ++i)
        trials_bound_[side][i][0] = trials_bound_1D_[i][side];
  }
  else
    assert( 0 == 1 );
    
/*  
  vector < vector<double> > db_helper = derivs_of_trial_at_boundaries(max_poly_degree);
  unsigned int num_of_one_dimensional_quad = num_of_quad_;
  

  
  trials_bound_.resize(2);
  for (unsigned int side = 0; side < 2; ++side)
  {
    trials_bound_[side].resize(max_poly_degree + 1);
    for (unsigned int i = 0; i < max_poly_degree + 1; ++i)
      trials_bound_[side][i] = vector<double>(1, trials_bound_1D_[i][side]);
  }
  
  // TODO: QUADRATURE POINTS AND WEIGHTS IN SEVERAL DIMENSIONS
  if (hyperedge_dim > 1)
  {
    vector<double> quad_weights_helper = quad_weights_; // quad weights remain a vector
     quad_bdr_ = quad_weights_;
    for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
    {
      if (dim == hyperedge_dim - 2)  quad_bdr_ = quad_weights_;
      quad_weights_ = dyadic_product(quad_weights_, quad_weights_helper);
    }
    num_of_quad_ = quad_weights_.size();
    num_quad_bdr_ = quad_bdr_.size();
    
    vector< vector<double> > trials_helper = trials_quad_; // trials at quad remain a vector (func) of vectors (quad)
    for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
      trials_quad_ = double_dyadic_product(trials_quad_, trials_helper);
    
    bound_trials_quad_ = trials_helper; // trials of lower dimensional functions are vector (func) of vectors (quad)
    for (unsigned int dim = 1; dim < hyperedge_dim - 1; ++dim)
      bound_trials_quad_ = double_dyadic_product(bound_trials_quad_, trials_helper);
    
    vector< vector<double> >derivs_helper = derivs_quad_[0]; // derivs at quad becomes of vector (deriv) of vectors (func) of vectors (quad)
    derivs_quad_.resize(hyperedge_dim);
    for (unsigned int dim_deriv = 0; dim_deriv < hyperedge_dim; ++dim_deriv)
    {
      derivs_quad_[dim_deriv] = trials_helper;
      for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
        if(dim_deriv == 0 && dim == 1)  derivs_quad_[dim_deriv] = double_dyadic_product(derivs_helper, trials_helper);
        else if (dim_deriv == dim)      derivs_quad_[dim_deriv] = double_dyadic_product(derivs_quad_[dim_deriv], derivs_helper);
        else                            derivs_quad_[dim_deriv] = double_dyadic_product(derivs_quad_[dim_deriv], trials_helper);
    }
    
    // trials at boundary becomes a vector (boundary) of vectors (function) of vectors (quad)
    trials_bound_.clear();
    trials_bound_.resize(2 * hyperedge_dim);
    for (unsigned int side = 0; side < 2; ++side)
    {
      vector< vector<double> > helperling(trials_bound_1D_.size());
      for (unsigned int i = 0; i < helperling.size(); ++i)
        helperling[i] = vector<double>(1, trials_bound_1D_[i][side]);
      
      for (unsigned int dim_face = 0; dim_face < hyperedge_dim; ++dim_face)
      {
        trials_bound_[2*dim_face+side] = trials_helper;
        for (unsigned int dim = 1; dim < hyperedge_dim; ++dim)
          if (dim_face == 0 && dim == 1)  trials_bound_[2*dim_face+side] = double_dyadic_product(helperling, trials_helper);
          else if (dim_face == dim)       trials_bound_[2*dim_face+side] = double_dyadic_product(trials_bound_[2*dim_face+side], helperling);
          else                            trials_bound_[2*dim_face+side] = double_dyadic_product(trials_bound_[2*dim_face+side], trials_helper);
      }
    }
  }
  */
  assert( quad_weights_.size() == num_of_quad_ );
//  assert( num_of_quad_ == pow(num_of_one_dimensional_quad,hyperedge_dim) );
  assert( quad_bdr_.size() == num_quad_bdr_ );
//  assert( num_quad_bdr_ == pow( num_of_one_dimensional_quad, hyperedge_dim - 1 ) );
  
  assert( trials_quad_.size() == num_ansatz_fct_ );
  for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
    assert( trials_quad_[ansatz].size() == num_of_quad_ );
  
  assert( bound_trials_quad_.size() == num_ansatz_bdr_ );
  for (unsigned int ansatz = 0; ansatz < num_ansatz_bdr_; ++ansatz)
    assert( bound_trials_quad_[ansatz].size() == num_quad_bdr_ );
  
  assert( derivs_quad_.size() == hyperedge_dim );
  for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
  {
    assert( derivs_quad_[dim].size() == num_ansatz_fct_ );
    for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
      assert( derivs_quad_[dim][ansatz].size() == num_of_quad_ );
  }
  
  assert( trials_bound_.size() == 2 * hyperedge_dim );
  for (unsigned int bdr = 0; bdr < 2 * hyperedge_dim; ++bdr)
  {
    cout << trials_bound_[bdr].size() << "  " << num_ansatz_fct_ << endl;
    assert( trials_bound_[bdr].size() == num_ansatz_fct_ );
    for (unsigned int ansatz = 0; ansatz < num_ansatz_bdr_; ++ansatz)
      assert( trials_bound_[bdr][ansatz].size() == num_quad_bdr_ );
  }
}


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
unsigned int DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
loc_matrix_index(const unsigned int row, const unsigned int column) const
{
  assert( 0 <= row && row < (hyperedge_dim + 1) * pow((max_poly_degree + 1), hyperedge_dim) );
  assert( 0 <= column && column < (hyperedge_dim + 1) * pow((max_poly_degree + 1), hyperedge_dim) );
  return column * (hyperedge_dim + 1) * pow((max_poly_degree + 1), hyperedge_dim) + row;  // Transposed for LAPACK
}


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_>
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
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


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array<double, (hyperedge_dim+1) * num_ansatz_fct_>
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
assemble_rhs(const array< array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const
{
  array<double, (hyperedge_dim+1) * num_ansatz_fct_> right_hand_side;
  right_hand_side.fill(0.);
  
  assert(lambda_values.size() == 2 * hyperedge_dim);
  for (unsigned int i = 0; i < 2 * hyperedge_dim; ++i)
    assert(lambda_values[i].size() == num_ansatz_bdr_);
  
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


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array<double, (hyperedge_dim+1) * num_ansatz_fct_>
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
solve_local_system_of_eq(array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_>& loc_matrix,
                         array<double, (hyperedge_dim+1) * num_ansatz_fct_>& loc_rhs) const
{
  assert( loc_matrix.size() == loc_rhs.size() * loc_rhs.size() );
  const int system_size = loc_rhs.size();
  double *mat_a=loc_matrix.data(), *rhs_b = loc_rhs.data();
  int info = -1;  
  lapack_solve(system_size, mat_a, rhs_b, &info);
  assert( info == 0 );
  return loc_rhs;
}


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array<double, (hyperedge_dim+1) * num_ansatz_fct_>
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
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


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
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


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim > 
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
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


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
auto // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
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


template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
array< array<double, local_dof_amount_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim > // array< array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
DiffusionSolver_RegularQuad<hyperedge_dim, max_poly_degree, max_quad_degree>::
numerical_flux_from_lambda(const array< array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const
{
  return numerical_flux_at_boundary(lambda_values, solve_local_problem(lambda_values));
}
