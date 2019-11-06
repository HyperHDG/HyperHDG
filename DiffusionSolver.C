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
#include "FuncAndQuad.h"
#include "LapackWrapper.h"
#include <cassert>
#include <cmath>
// #include <Eigen/Dense>


#include <iostream>


using namespace std;


template class DiffusionSolver_RegularQuad<1>;
template class DiffusionSolver_RegularQuad<2>;
template class DiffusionSolver_RegularQuad<3>;


vector<double> dyadic_product(const vector<double>& left, const vector<double>& right)
{
  vector<double> result(left.size() * right.size());
  for (unsigned int i = 0; i < right.size(); ++i)
    for (unsigned int j = 0; j < left.size(); ++j)
      result[i * left.size() + j] = left[j] * right[i];
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


vector<double> get_relevant_coeffs_indicator(const unsigned int connector_dim, const unsigned int max_poly_degree, const unsigned int dimension, const unsigned int ansatz)
{
  vector<double> unity_vec(max_poly_degree+1, 1.);
  vector<double> ansatz_vec(max_poly_degree+1, 0.);
  assert( ansatz < ansatz_vec.size() );
  ansatz_vec[ansatz] = 1;
  vector<double> result;
  if( connector_dim == 1 )  result = unity_vec;
  else                      result = ansatz_vec;
  
  for (unsigned int dim = 1; dim < connector_dim; ++dim)
    if (dim == 1 && dimension == 0)  result = dyadic_product(unity_vec, result);
    else if (dim == dimension)       result = dyadic_product(result, unity_vec);
    else                             result = dyadic_product(result, ansatz_vec);
  
  return result;
}


template<unsigned int connector_dim>
DiffusionSolver_RegularQuad<connector_dim>::DiffusionSolver_RegularQuad
(const unsigned int max_poly_degree, const unsigned int num_of_quad, const double tau)
: max_poly_degree_(max_poly_degree), num_of_quad_(num_of_quad), num_quad_bdr_(1),
  tau_(tau), quad_weights_(quadrature_weights(num_of_quad)), quad_bdr_(1,1.),
  trials_quad_(trial_functions_at_quadrature_points(max_poly_degree, num_of_quad)),
  bound_trials_quad_(1,vector<double>(1,1.)),
  trials_bound_1D_(trial_functions_at_boundaries(max_poly_degree)),
  derivs_quad_(1,derivs_of_trial_at_quadrature_points(max_poly_degree, num_of_quad))
{ 
  vector < vector<double> > db_helper = derivs_of_trial_at_boundaries(max_poly_degree);
  unsigned int num_of_one_dimensional_quad = num_of_quad_;
/*  
  trials_bound_.resize(2);
  derivs_bound_.resize(2);
  for (unsigned int side = 0; side < 2; ++side)
  {
    trials_bound_[side].resize(1);
    derivs_bound_[side].resize(1);
    trials_bound_[side][0].resize(1);
    derivs_bound_[side][0].resize(1);
    trials_bound_[side][0][0] = trials_bound_1D_[0][side];
    derivs_bound_[side][0][0] = db_helper[0][side];
  }
  
  cout << trials_bound_[0][0][0] << "  " << trials_bound_[1][0][0] << endl;
*/  
  
  const unsigned int num_ansatz_fct = pow ( max_poly_degree_ + 1, connector_dim );
  const unsigned int num_ansatz_bdr = pow ( max_poly_degree_ + 1, connector_dim - 1 );
  const unsigned int num_1D_quad = num_of_quad_;
  
  trials_bound_.resize(2);
  for (unsigned int side = 0; side < 2; ++side)
  {
    trials_bound_[side].resize(max_poly_degree_ + 1);
    for (unsigned int i = 0; i < max_poly_degree_ + 1; ++i)
      trials_bound_[side][i] = vector<double>(1, trials_bound_1D_[i][side]);
  }
  
  // TODO: QUADRATURE POINTS AND WEIGHTS IN SEVERAL DIMENSIONS
  if (connector_dim > 1)
  {
    vector<double> quad_weights_helper = quad_weights_; // quad weights remain a vector
     quad_bdr_ = quad_weights_;
    for (unsigned int dim = 1; dim < connector_dim; ++dim)
    {
      if (dim == connector_dim - 2)  quad_bdr_ = quad_weights_;
      quad_weights_ = dyadic_product(quad_weights_, quad_weights_helper);
    }
    num_of_quad_ = quad_weights_.size();
    num_quad_bdr_ = quad_bdr_.size();
    
    vector< vector<double> > trials_helper = trials_quad_; // trials at quad remain a vector (func) of vectors (quad)
    for (unsigned int dim = 1; dim < connector_dim; ++dim)
      trials_quad_ = double_dyadic_product(trials_quad_, trials_helper);
    
    bound_trials_quad_ = trials_helper; // trials of lower dimensional functions are vector (func) of vectors (quad)
    for (unsigned int dim = 1; dim < connector_dim - 1; ++dim)
      bound_trials_quad_ = double_dyadic_product(bound_trials_quad_, trials_helper);
    
    vector< vector<double> >derivs_helper = derivs_quad_[0]; // derivs at quad becomes of vector (deriv) of vectors (func) of vectors (quad)
    derivs_quad_.resize(connector_dim);
    for (unsigned int dim_deriv = 0; dim_deriv < connector_dim; ++dim_deriv)
    {
      derivs_quad_[dim_deriv] = trials_helper;
      for (unsigned int dim = 1; dim < connector_dim; ++dim)
        if(dim_deriv == 0 && dim == 1)  derivs_quad_[dim_deriv] = double_dyadic_product(derivs_helper, trials_helper);
        else if (dim_deriv == dim)      derivs_quad_[dim_deriv] = double_dyadic_product(derivs_quad_[dim_deriv], derivs_helper);
        else                            derivs_quad_[dim_deriv] = double_dyadic_product(derivs_quad_[dim_deriv], trials_helper);
    }
    
    // trials at boundary becomes a vector (boundary) of vectors (function) of vectors (quad)
    trials_bound_.clear();
    trials_bound_.resize(2 * connector_dim);
    for (unsigned int side = 0; side < 2; ++side)
    {
      vector< vector<double> > helperling(trials_bound_1D_.size());
      for (unsigned int i = 0; i < helperling.size(); ++i)
        helperling[i] = vector<double>(1, trials_bound_1D_[i][side]);
      
      for (unsigned int dim_face = 0; dim_face < connector_dim; ++dim_face)
      {
        trials_bound_[2*dim_face+side] = trials_helper;
        for (unsigned int dim = 1; dim < connector_dim; ++dim)
          if (dim_face == 0 && dim == 1)  trials_bound_[2*dim_face+side] = double_dyadic_product(helperling, trials_helper);
          else if (dim_face == dim)       trials_bound_[2*dim_face+side] = double_dyadic_product(trials_bound_[2*dim_face+side], helperling);
          else                            trials_bound_[2*dim_face+side] = double_dyadic_product(trials_bound_[2*dim_face+side], trials_helper);
      }
    }
  }
  
  assert( quad_weights_.size() == num_of_quad_ );
  assert( num_of_quad_ == pow(num_of_one_dimensional_quad,connector_dim) );
  assert( quad_bdr_.size() == num_quad_bdr_ );
  assert( num_quad_bdr_ == pow( num_of_one_dimensional_quad, connector_dim - 1 ) );
  
  assert( trials_quad_.size() == num_ansatz_fct );
  for (unsigned int ansatz = 0; ansatz < num_ansatz_fct; ++ansatz)
    assert( trials_quad_[ansatz].size() == num_of_quad_ );
  
  assert( bound_trials_quad_.size() == num_ansatz_bdr );
  for (unsigned int ansatz = 0; ansatz < num_ansatz_bdr; ++ansatz)
    assert( bound_trials_quad_[ansatz].size() == num_quad_bdr_ );
  
  assert( derivs_quad_.size() == connector_dim );
  for (unsigned int dim = 0; dim < connector_dim; ++dim)
  {
    assert( derivs_quad_[dim].size() == num_ansatz_fct );
    for (unsigned int ansatz = 0; ansatz < num_ansatz_fct; ++ansatz)
      assert( derivs_quad_[dim][ansatz].size() == num_of_quad_ );
  }
  
  assert( trials_bound_.size() == 2 * connector_dim );
  for (unsigned int bdr = 0; bdr < 2 * connector_dim; ++bdr)
  {
    assert( trials_bound_[bdr].size() == num_ansatz_fct );
    for (unsigned int ansatz = 0; ansatz < num_ansatz_bdr; ++ansatz)
      assert( trials_bound_[bdr][ansatz].size() == num_quad_bdr_ );
  }
}


template<unsigned int connector_dim>
unsigned int DiffusionSolver_RegularQuad<connector_dim>::loc_matrix_index(const unsigned int row, const unsigned int column) const
{
  assert( 0 <= row && row < (connector_dim + 1) * pow((max_poly_degree_ + 1), connector_dim) );
  assert( 0 <= column && column < (connector_dim + 1) * pow((max_poly_degree_ + 1), connector_dim) );
  return column * (connector_dim + 1) * pow((max_poly_degree_ + 1), connector_dim) + row;  // Transposed for LAPACK
}


template<unsigned int connector_dim>
vector<double> DiffusionSolver_RegularQuad<connector_dim>::assemble_loc_mat() const
{
  double connector_area = 1.;
  unsigned int num_ansatz_fct = pow ( max_poly_degree_ + 1, connector_dim );
  unsigned int num_ansatz_bdr = pow ( max_poly_degree_ + 1, connector_dim - 1 );
  vector<double> local_mat( (connector_dim+1) * num_ansatz_fct * (connector_dim+1) * num_ansatz_fct , 0.);
  
  for (unsigned int dim = 0; dim < connector_dim; ++dim)
    for (unsigned int i = 0; i < num_ansatz_fct; ++i)
    {
      for (unsigned int j = 0; j < num_ansatz_fct; ++j)
        for (unsigned int q = 0; q < num_of_quad_; ++q)
          local_mat[loc_matrix_index( dim * num_ansatz_fct + i , dim * num_ansatz_fct + j )] += 
            quad_weights_[q] * connector_area * trials_quad_[i][q] * trials_quad_[j][q];
      for (unsigned int j = 0; j < num_ansatz_fct; ++j)
        for (unsigned int q = 0; q < num_of_quad_; ++q)
          local_mat[loc_matrix_index(  dim * num_ansatz_fct + i , connector_dim * num_ansatz_fct + j )] -=
            quad_weights_[q] * derivs_quad_[dim][i][q] * trials_quad_[j][q];
    }
  
  for (unsigned int i = 0; i < num_ansatz_fct; ++i)
  {
    for (unsigned int dim = 0; dim < connector_dim; ++dim)
    {
      for (unsigned int j = 0; j < num_ansatz_fct; ++j)
      {
        for (unsigned int q = 0; q < num_of_quad_; ++q)
          local_mat[loc_matrix_index( connector_dim * num_ansatz_fct + i , dim * num_ansatz_fct + j )] -=
            quad_weights_[q] * derivs_quad_[dim][i][q] * trials_quad_[j][q];
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          local_mat[loc_matrix_index( connector_dim * num_ansatz_fct + i , dim * num_ansatz_fct + j )] +=
            + quad_bdr_[q] * trials_bound_[2*dim+1][i][q] * trials_bound_[2*dim+1][j][q]
            - quad_bdr_[q] * trials_bound_[2*dim+0][i][q] * trials_bound_[2*dim+0][j][q];
      }
      for (unsigned int j = 0; j < num_ansatz_fct; ++j)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          local_mat[loc_matrix_index( connector_dim * num_ansatz_fct + i , connector_dim * num_ansatz_fct + j )] +=
            tau_ * quad_bdr_[q] * (trials_bound_[2*dim+0][i][q] * trials_bound_[2*dim+0][j][q]
                                   + trials_bound_[2*dim+1][i][q] * trials_bound_[2*dim+1][j][q]);
    }
  }
/*  
  cout << endl << endl;
  for (unsigned int i = 0; i < (connector_dim + 1) * pow((max_poly_degree_ + 1), connector_dim); ++i)
  {
    for(unsigned int j = 0; j < (connector_dim + 1) * pow((max_poly_degree_ + 1), connector_dim); ++j)
      cout << local_mat[loc_matrix_index(i,j)] << "  ";
    cout << endl;
  }
*/  
  return local_mat;
}


template<unsigned int connector_dim>
vector<double> DiffusionSolver_RegularQuad<connector_dim>::assemble_rhs(const vector< vector<double> >& lambda_values) const
{
  unsigned int num_ansatz_fct = pow ( max_poly_degree_ + 1, connector_dim );
  unsigned int num_ansatz_bdr = pow ( max_poly_degree_ + 1, connector_dim - 1 );
  vector<double> right_hand_side((connector_dim+1) * num_ansatz_fct, 0.);
  
  assert(lambda_values.size() == 2 * connector_dim);
  for (unsigned int i = 0; i < 2 * connector_dim; ++i)
    assert(lambda_values[i].size() == num_ansatz_bdr);
  
  for (unsigned int dim = 0; dim < connector_dim; ++dim)
    for (unsigned int ansatz = 0; ansatz < num_ansatz_fct; ++ansatz)
      for (unsigned int bdr_ans = 0; bdr_ans < num_ansatz_bdr; ++bdr_ans)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          right_hand_side[dim * num_ansatz_fct + ansatz] += 
            + quad_bdr_[q] * lambda_values[2*dim+0][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+0][ansatz][q]
            - quad_bdr_[q] * lambda_values[2*dim+1][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+1][ansatz][q];
  
  for (unsigned int ansatz = 0; ansatz < num_ansatz_fct; ++ansatz)
    for (unsigned int dim = 0; dim < connector_dim; ++dim)
      for (unsigned int bdr_ans = 0; bdr_ans < num_ansatz_bdr; ++bdr_ans)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          right_hand_side[connector_dim * num_ansatz_fct + ansatz] += tau_ * quad_bdr_[q] *
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


template<unsigned int connector_dim>
vector<double> DiffusionSolver_RegularQuad<connector_dim>::solve_local_system_of_eq(vector<double>& loc_matrix, vector<double>& loc_rhs) const
{
  assert( loc_matrix.size() == loc_rhs.size() * loc_rhs.size() );
  const int system_size = loc_rhs.size();
  double *mat_a=loc_matrix.data(), *rhs_b = loc_rhs.data();
  int info = -1;  
  lapack_solve(system_size, mat_a, rhs_b, &info);
  assert( info == 0 );
  return loc_rhs;
}


template<unsigned int connector_dim>
vector<double> DiffusionSolver_RegularQuad<connector_dim>::solve_local_problem(const vector< vector<double> >& lambda_values) const
{
  vector<double> right_hand_side, local_matrix;
  right_hand_side = assemble_rhs(lambda_values);
  local_matrix = assemble_loc_mat();
  
  vector<double> solution = solve_local_system_of_eq(local_matrix, right_hand_side);
/*  
  cout << endl;
  for (unsigned int i = 0; i < solution.size(); ++i) cout << "  " << solution[i];
  cout << endl;
*/  
  return solution;
}


template<unsigned int connector_dim>
vector< vector<double> > DiffusionSolver_RegularQuad<connector_dim>::dual_at_boundary(const vector<double>& coeffs) const
{
  unsigned int num_ansatz_fct = pow ( max_poly_degree_ + 1, connector_dim );
  unsigned int num_ansatz_bdr = pow ( max_poly_degree_ + 1, connector_dim - 1 );
  
  vector< vector<double> > bdr_values;
  bdr_values.resize(2 * connector_dim);
  
  for (unsigned int dim = 0; dim < connector_dim; ++dim)
  {
    bdr_values[2*dim].resize(num_ansatz_bdr);
    bdr_values[2*dim+1].resize(num_ansatz_bdr);
    for(unsigned int bdr_ansatz = 0; bdr_ansatz < num_ansatz_bdr; ++bdr_ansatz)
      for (unsigned int ansatz = 0; ansatz < num_ansatz_fct; ++ansatz)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
        {
          bdr_values[2*dim+0][bdr_ansatz] -= quad_bdr_[q] * coeffs[dim * num_ansatz_fct + ansatz] * trials_bound_[2*dim+0][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
          bdr_values[2*dim+1][bdr_ansatz] += quad_bdr_[q] * coeffs[dim * num_ansatz_fct + ansatz] * trials_bound_[2*dim+1][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
        }
  }
  
  return bdr_values;
}


template<unsigned int connector_dim>
vector< vector<double> > DiffusionSolver_RegularQuad<connector_dim>::primal_at_boundary(const vector<double>& coeffs) const
{
  unsigned int num_ansatz_fct = pow ( max_poly_degree_ + 1, connector_dim );
  unsigned int num_ansatz_bdr = pow ( max_poly_degree_ + 1, connector_dim - 1 );
  
  vector<double> relevant_dofs;
  vector< vector<double> > bdr_values;
  bdr_values.resize(2 * connector_dim);
  
  for (unsigned int dim = 0; dim < connector_dim; ++dim)
  {
    bdr_values[2*dim].resize(num_ansatz_bdr);
    bdr_values[2*dim+1].resize(num_ansatz_bdr);
    for (unsigned int bdr_ansatz = 0; bdr_ansatz < num_ansatz_bdr; ++bdr_ansatz)
      for (unsigned int ansatz = 0; ansatz < num_ansatz_fct; ++ansatz)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
        {
          bdr_values[2*dim+0][bdr_ansatz] += quad_bdr_[q] * coeffs[connector_dim * num_ansatz_fct + ansatz] * trials_bound_[2*dim+0][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
          bdr_values[2*dim+1][bdr_ansatz] += quad_bdr_[q] * coeffs[connector_dim * num_ansatz_fct + ansatz] * trials_bound_[2*dim+1][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
        }
  }
  
  return bdr_values;
}


template<unsigned int connector_dim>
vector< vector<double> > DiffusionSolver_RegularQuad<connector_dim>::numerical_flux_at_boundary(const vector< vector<double> >& lambda_values, const vector<double>& coeffs) const
{
  vector< vector<double> > bdr_values;
  vector< vector<double> > primals = primal_at_boundary(coeffs);
  vector< vector<double> > duals   = dual_at_boundary(coeffs);
  
  bdr_values.resize(lambda_values.size());
  for (unsigned int i = 0; i < lambda_values.size(); ++i)
  {
	bdr_values[i].resize(lambda_values[i].size());
    for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
      bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];
  }
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


template<unsigned int connector_dim>
vector< vector<double> > DiffusionSolver_RegularQuad<connector_dim>::numerical_flux_from_lambda(const vector< vector<double> >& lambda_values) const
{
  return numerical_flux_at_boundary(lambda_values, solve_local_problem(lambda_values));
}

/*
template<unsigned int dim, unsigned int unknown_dim>
void DiffusionSolver<dim,unknown_dim>::update_auxiliaries(Connector<dim,unknown_dim>& connector, const bool hom_rhs) const
{
  vector<double> fluxes = numerical_flux_from_lambda(connector, hom_rhs);
  
  if (connector.get_left().get_joint_type() != 1)  connector.get_left().aux_dof(0) += fluxes[0];
  if (connector.get_right().get_joint_type() != 1)  connector.get_right().aux_dof(0) += fluxes[1];
}
*/
