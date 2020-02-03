/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "LSol_Diffusion.hxx"
#include "HyAssert.hxx"
#include "LapackWrapper.hxx"
#include <cmath>
#include <vector>


#include <iostream>


using namespace std;
using namespace FuncQuad;
#include "LSol_Diffusion.inst"


template<unsigned int hyEdge_dim, unsigned int poly_deg>
inline unsigned int loc_matrix_index(const unsigned int row, const unsigned int column)
{
  hy_assert( 0 <= row ,
             "Row index should be larger than or equal to zero." );
  hy_assert( row < (hyEdge_dim + 1) * pow((poly_deg + 1), hyEdge_dim) ,
             "Row index should be smaller than total amount of rows." );
  hy_assert( 0 <= column ,
             "Column index should be larger than or equal to zero." );
  hy_assert( column < (hyEdge_dim + 1) * pow((poly_deg + 1), hyEdge_dim) ,
             "Column index should smaller than total amount of columns." );
  return column * (hyEdge_dim + 1) * pow((poly_deg + 1), hyEdge_dim) + row;  // Transposed for LAPACK
}

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


inline vector<double> get_relevant_coeffs_indicator(const unsigned int hyEdge_dim, const unsigned int poly_deg, const unsigned int dimension, const unsigned int ansatz)
{
  vector<double> unity_vec(poly_deg+1, 1.);
  vector<double> ansatz_vec(poly_deg+1, 0.);
  hy_assert( ansatz < ansatz_vec.size() , "Ansatz function index in one dimension must be smaller than maximal degree + 1." );
  ansatz_vec[ansatz] = 1;
  vector<double> result;
  if( hyEdge_dim == 1 )  result = unity_vec;
  else                      result = ansatz_vec;
  
  for (unsigned int dim = 1; dim < hyEdge_dim; ++dim)
    if (dim == 1 && dimension == 0)  result = dyadic_product(unity_vec, result);
    else if (dim == dimension)       result = dyadic_product(result, unity_vec);
    else                             result = dyadic_product(result, ansatz_vec);
  
  return result;
}

/*
template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
DiffusionSolverNaive_RegularQuad(const constructor_value_type& tau)
: tau_(tau)
{ 
  array<double, compute_n_quad_points(quad_deg)>
    quad_weights1D = quad_weights<quad_deg>();
  array< array<double, compute_n_quad_points(quad_deg)> , poly_deg + 1 >
    trials_at_quad1D = shape_fcts_at_quad_points<poly_deg, quad_deg>();
  array< array<double, compute_n_quad_points(quad_deg)> , poly_deg + 1 >
    derivs_at_quad1D = shape_ders_at_quad_points<poly_deg, quad_deg>();
  array< array<double, 2> , poly_deg + 1 >
    trials_at_bdr1D = shape_fcts_at_bdrs<poly_deg>();
//  array< array<double, 2> , poly_deg + 1 >
//    derivs_at_bdr1D = shape_ders_at_bdrs<poly_deg>();
    
  // In the one-dimensional case, we are done now.
  if constexpr (hyEdge_dim == 1)
  {
    quad_weights_ = quad_weights1D;
    quad_bdr_ = { 1. };
    trials_quad_ = trials_at_quad1D;
    bound_trials_quad_[0] = { 1. };
    derivs_quad_[0] = derivs_at_quad1D;
    trials_in_corners_ = trials_at_bdr1D;
    for (unsigned int side = 0; side < 2; ++side)
      for (unsigned int i = 0; i < poly_deg + 1; ++i)
        trials_bound_[side][i][0] = trials_at_bdr1D[i][side];
  }
  else // There is more to be done and we need to use growing containers.
  {
    // Deal with quadrature points that remain a vector.
    vector<double> quad_weights1D_vec(quad_weights1D.begin(), quad_weights1D.end());
    vector<double> quad_weights_vec(quad_weights1D.begin(), quad_weights1D.end());
    if (quad_bdr_.size() == quad_weights1D_vec.size())
      for (unsigned int i = 0; i < quad_bdr_.size(); ++i)  quad_bdr_[i] = quad_weights1D_vec[i];
    for (unsigned int dim = 1; dim < hyEdge_dim; ++dim)
    {
      if (dim == hyEdge_dim - 2)
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
    for (unsigned int dim = 1; dim < hyEdge_dim; ++dim)
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
    for (unsigned int dim = 1; dim < hyEdge_dim; ++dim)
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
    for (unsigned int dim = 1; dim < hyEdge_dim - 1; ++dim)
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
    vector< vector< vector<double> > > derivs_at_quad_vec(hyEdge_dim);
    for (unsigned int dim_deriv = 0; dim_deriv < hyEdge_dim; ++dim_deriv)
    {
      derivs_at_quad_vec[dim_deriv] = trials_at_quad1D_vec;
      for (unsigned int dim = 1; dim < hyEdge_dim; ++dim)
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
    vector< vector< vector<double> > > trials_bound_vec(2 * hyEdge_dim);
    for (unsigned int side = 0; side < 2; ++side)
    {
      vector< vector<double> > helperling(trials_at_bdr1D.size());
      for (unsigned int i = 0; i < helperling.size(); ++i)
        helperling[i] = vector<double>(1, trials_at_bdr1D[i][side]);
      for (unsigned int dim_face = 0; dim_face < hyEdge_dim; ++dim_face)
      {
        trials_bound_vec[2*dim_face+side] = trials_at_quad1D_vec;
        for (unsigned int dim = 1; dim < hyEdge_dim; ++dim)
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



template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
inline auto // array<double, (hyEdge_dim+1) * num_ansatz_fct_ * (hyEdge_dim+1) * num_ansatz_fct_>
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
assemble_loc_mat() const
{
  double hyEdge_area = 1.;
  array<double, (hyEdge_dim+1) * num_ansatz_fct_ * (hyEdge_dim+1) * num_ansatz_fct_> local_mat;
  local_mat.fill(0.);
  
  for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
    for (unsigned int i = 0; i < num_ansatz_fct_; ++i)
    {
      for (unsigned int j = 0; j < num_ansatz_fct_; ++j)
        for (unsigned int q = 0; q < n_quads_; ++q)
          local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( dim * num_ansatz_fct_ + i , dim * num_ansatz_fct_ + j )] += 
            quad_weights_[q] * hyEdge_area * trials_quad_[i][q] * trials_quad_[j][q];
      for (unsigned int j = 0; j < num_ansatz_fct_; ++j)
        for (unsigned int q = 0; q < n_quads_; ++q)
          local_mat[loc_matrix_index<hyEdge_dim,poly_deg>(  dim * num_ansatz_fct_ + i , hyEdge_dim * num_ansatz_fct_ + j )] -=
            quad_weights_[q] * derivs_quad_[dim][i][q] * trials_quad_[j][q];
    }
  
  for (unsigned int i = 0; i < num_ansatz_fct_; ++i)
  {
    for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
    {
      for (unsigned int j = 0; j < num_ansatz_fct_; ++j)
      {
        for (unsigned int q = 0; q < n_quads_; ++q)
          local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * num_ansatz_fct_ + i , dim * num_ansatz_fct_ + j )] -=
            quad_weights_[q] * derivs_quad_[dim][i][q] * trials_quad_[j][q];
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * num_ansatz_fct_ + i , dim * num_ansatz_fct_ + j )] +=
            + quad_bdr_[q] * trials_bound_[2*dim+1][i][q] * trials_bound_[2*dim+1][j][q]
            - quad_bdr_[q] * trials_bound_[2*dim+0][i][q] * trials_bound_[2*dim+0][j][q];
      }
      for (unsigned int j = 0; j < num_ansatz_fct_; ++j)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * num_ansatz_fct_ + i , hyEdge_dim * num_ansatz_fct_ + j )] +=
            tau_ * quad_bdr_[q] * (trials_bound_[2*dim+0][i][q] * trials_bound_[2*dim+0][j][q]
                                   + trials_bound_[2*dim+1][i][q] * trials_bound_[2*dim+1][j][q]);
    }
  }

  return local_mat;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
inline auto // array<double, (hyEdge_dim+1) * num_ansatz_fct_>
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
assemble_rhs(const array< array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * num_ansatz_fct_> right_hand_side;
  right_hand_side.fill(0.);
  
  hy_assert( lambda_values.size() == 2 * hyEdge_dim ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
    hy_assert( lambda_values[i].size() == num_ansatz_bdr_ ,
               "The size of the lambda values should be the amount of ansatz functions ar boundary." );
  
  for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
    for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
      for (unsigned int bdr_ans = 0; bdr_ans < num_ansatz_bdr_; ++bdr_ans)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          right_hand_side[dim * num_ansatz_fct_ + ansatz] += 
            + quad_bdr_[q] * lambda_values[2*dim+0][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+0][ansatz][q]
            - quad_bdr_[q] * lambda_values[2*dim+1][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+1][ansatz][q];
  
  for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
    for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      for (unsigned int bdr_ans = 0; bdr_ans < num_ansatz_bdr_; ++bdr_ans)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
          right_hand_side[hyEdge_dim * num_ansatz_fct_ + ansatz] += tau_ * quad_bdr_[q] *
            ( lambda_values[2*dim+0][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+0][ansatz][q]
            + lambda_values[2*dim+1][bdr_ans] * bound_trials_quad_[bdr_ans][q] * trials_bound_[2*dim+1][ansatz][q] );

  return right_hand_side;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
inline auto // array<double, (hyEdge_dim+1) * num_ansatz_fct_>
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
solve_local_problem(const array< array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * num_ansatz_fct_ * (hyEdge_dim+1) * num_ansatz_fct_> local_matrix;
  array<double, (hyEdge_dim+1) * num_ansatz_fct_> right_hand_side;
  right_hand_side = assemble_rhs(lambda_values);
  local_matrix = assemble_loc_mat();
  
  array<double, (hyEdge_dim+1) * num_ansatz_fct_> solution = lapack_solve<(hyEdge_dim+1) * num_ansatz_fct_>(local_matrix, right_hand_side);
  return solution;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
inline auto // array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
dual_at_boundary(const array<double, (hyEdge_dim+1) * num_ansatz_fct_>& coeffs) const
{
  array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim > bdr_values;
  
  for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
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


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
inline auto // array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim > 
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
primal_at_boundary(const array<double, (hyEdge_dim+1) * num_ansatz_fct_>& coeffs) const
{
  array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim > bdr_values;
  
  for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
  {
    bdr_values[2*dim].fill(0.);
    bdr_values[2*dim+1].fill(0.);
    for (unsigned int bdr_ansatz = 0; bdr_ansatz < num_ansatz_bdr_; ++bdr_ansatz)
      for (unsigned int ansatz = 0; ansatz < num_ansatz_fct_; ++ansatz)
        for (unsigned int q = 0; q < num_quad_bdr_; ++q)
        {
          bdr_values[2*dim+0][bdr_ansatz] += quad_bdr_[q] * coeffs[hyEdge_dim * num_ansatz_fct_ + ansatz] * trials_bound_[2*dim+0][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
          bdr_values[2*dim+1][bdr_ansatz] += quad_bdr_[q] * coeffs[hyEdge_dim * num_ansatz_fct_ + ansatz] * trials_bound_[2*dim+1][ansatz][q] * bound_trials_quad_[bdr_ansatz][q];
        }
  }
  
  return bdr_values;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
auto // array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
numerical_flux_at_boundary(const array< array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values, const array<double, (hyEdge_dim+1) * num_ansatz_fct_>& coeffs) const
{
  array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim > bdr_values;
  array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim > primals = primal_at_boundary(coeffs);
  array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim > duals   = dual_at_boundary(coeffs);
  
  for (unsigned int i = 0; i < lambda_values.size(); ++i)
    for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
      bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];

  return bdr_values;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
array<double, Hypercube<hyEdge_dim>::n_vertices()>
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
primal_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * num_ansatz_fct_> coefficients = solve_local_problem(lambda_values);
  array<double, Hypercube<hyEdge_dim>::n_vertices()> primal_in_corners;
  primal_in_corners.fill(0.);
  for (unsigned int corner = 0; corner < Hypercube<hyEdge_dim>::n_vertices(); ++corner)
    for (unsigned int ansatz_fct = 0; ansatz_fct < num_ansatz_fct_; ++ansatz_fct)
      primal_in_corners[corner] += coefficients[hyEdge_dim * num_ansatz_fct_ + ansatz_fct] * trials_in_corners_[ansatz_fct][corner];
  return primal_in_corners;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
array< array<double, hyEdge_dim> , Hypercube<hyEdge_dim>::n_vertices() >
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
dual_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * num_ansatz_fct_> coefficients = solve_local_problem(lambda_values);
  array< array<double, hyEdge_dim> , Hypercube<hyEdge_dim>::n_vertices() > dual_in_corners;
  for (unsigned int corner = 0; corner < Hypercube<hyEdge_dim>::n_vertices(); ++corner)
  {
    dual_in_corners[corner].fill(0.);
    for (unsigned int ansatz_fct = 0; ansatz_fct < num_ansatz_fct_; ++ansatz_fct)
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
        dual_in_corners[corner][dim] += coefficients[dim * num_ansatz_fct_ + ansatz_fct] * trials_in_corners_[ansatz_fct][corner];
  }
  return dual_in_corners;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
auto // array< array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
DiffusionSolverNaive_RegularQuad<hyEdge_dim, poly_deg, quad_deg>::
numerical_flux_from_lambda(const array< array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  return numerical_flux_at_boundary(lambda_values, solve_local_problem(lambda_values));
}
*/



















template <unsigned int hyEdge_dim>
inline void index_decompose(unsigned int index, unsigned int range, array<unsigned int, max(hyEdge_dim,1U)>& decomposition)
{
  if ( decomposition.size() == 1 )  decomposition[0] = index;
  else
  {
    for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
    {
      decomposition[dim] = index % range;
      index /= range;
    }
  }
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg, unsigned int n_dofs_per_node>
array<double, (hyEdge_dim+1) * n_dofs_per_node * (poly_deg + 1) * (hyEdge_dim+1) * n_dofs_per_node * (poly_deg + 1)>
assemble_loc_matrix(const double tau)
{ 
  const unsigned int n_quads = compute_n_quad_points(quad_deg);
  const unsigned int n_shape_fct = n_dofs_per_node * (poly_deg + 1);
  const std::array<double, n_quads> q_weights = quad_weights<quad_deg>();
  const std::array< std::array<double, n_quads > , poly_deg + 1 > trial = shape_fcts_at_quad_points<poly_deg, quad_deg>();
  const std::array< std::array<double, n_quads > , poly_deg + 1 > deriv = shape_ders_at_quad_points<poly_deg, quad_deg>();
  const std::array< std::array<double, 2> , poly_deg + 1 > trial_bdr = shape_fcts_at_bdrs<poly_deg>();
  
  array<unsigned int, hyEdge_dim> dec_i, dec_j;
  double integral, integral1D;
  
  array<double, (hyEdge_dim+1) * n_shape_fct * (hyEdge_dim+1) * n_shape_fct> local_mat;
  local_mat.fill(0.);
  
  for (unsigned int i = 0; i < n_shape_fct; ++i)
  {
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int j = 0; j < n_shape_fct; ++j)
    {
      index_decompose<hyEdge_dim>(j, poly_deg+1, dec_j);
      
      // Integral_element phi_i phi_j dx in diagonal blocks
      integral = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
      {
        integral1D = 0.;
        for (unsigned int q = 0; q < n_quads; ++q)
          integral1D += q_weights[q] * trial[dec_j[dim_fct]][q] * trial[dec_i[dim_fct]][q];
        integral *= integral1D;
      }
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( dim * n_shape_fct + i , dim * n_shape_fct + j )] += integral;
      
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      { 
        // Integral_element - nabla phi_i \vec phi_j dx = Integral_element - div \vec phi_i phi_j dx
        // in right upper and left lower blocks
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          integral1D = 0.;
          for (unsigned int q = 0; q < n_quads; ++q)
            integral1D += q_weights[q] * trial[dec_j[dim_fct]][q] * 
                          ( ( dim == dim_fct ) ? deriv[dec_i[dim_fct]][q] : trial[dec_i[dim_fct]][q] );
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct + i , dim * n_shape_fct + j )] -= integral;
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>(  dim * n_shape_fct + i , hyEdge_dim * n_shape_fct + j )] -= integral;
        
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr[dec_i[dim_fct]][1] * trial_bdr[dec_j[dim_fct]][1];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads; ++q)
              integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct + i , dim * n_shape_fct + j )] += integral;
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr[dec_i[dim_fct]][0] * trial_bdr[dec_j[dim_fct]][0];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads; ++q)
              integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct + i , dim * n_shape_fct + j )] -= integral;
        
        // Penalty in lower right diagonal block
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr[dec_i[dim_fct]][0] * trial_bdr[dec_j[dim_fct]][0];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads; ++q)
              integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct + i , hyEdge_dim * n_shape_fct + j )] += tau * integral;
        // Penalty in lower right diagonal block
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr[dec_i[dim_fct]][1] * trial_bdr[dec_j[dim_fct]][1];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads; ++q)
              integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct + i , hyEdge_dim * n_shape_fct + j )] += tau * integral;
      }
    }
  }
  
  return local_mat;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::
DiffusionSolverTensorStruc(const constructor_value_type& tau)
: tau_(tau), q_weights_(quad_weights<quad_deg>()),
  trial_(shape_fcts_at_quad_points<poly_deg, quad_deg>()),
  trial_bdr_(shape_fcts_at_bdrs<poly_deg>()),
  loc_mat_(assemble_loc_matrix<hyEdge_dim,poly_deg,quad_deg, n_glob_dofs_per_node()>(tau))
{ } 


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
array<double, (hyEdge_dim+1) * DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::n_shape_fct_>
DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::
assemble_rhs(const array< array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<unsigned int, hyEdge_dim> dec_i;
  array<unsigned int, max(hyEdge_dim-1,1U)> dec_j;
  double integral, integral1D;
  
  array<double, (hyEdge_dim+1) * n_shape_fct_> right_hand_side;
  right_hand_side.fill(0.);
  
  hy_assert( lambda_values.size() == 2 * hyEdge_dim ,
             "The size of the lambda values should be twice the dimension of a hyperedge." );
  for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
    hy_assert( lambda_values[i].size() == n_shape_bdr_ ,
               "The size of the lambda values should be the amount of ansatz functions ar boundary." );
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      {
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_; ++q)
              integral1D += q_weights_[q] * trial_[dec_i[dim_fct]][q] * trial_[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        right_hand_side[dim * n_shape_fct_ + i] += lambda_values[2*dim+0][j] * integral;
        right_hand_side[hyEdge_dim * n_shape_fct_ + i] += tau_ * lambda_values[2*dim+0][j] * integral;
        
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_; ++q)
              integral1D += q_weights_[q] * trial_[dec_i[dim_fct]][q] * trial_[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        right_hand_side[dim * n_shape_fct_ + i] -= lambda_values[2*dim+1][j] * integral;
        right_hand_side[hyEdge_dim * n_shape_fct_ + i] += tau_ * lambda_values[2*dim+1][j] * integral;
      }
    }
  }
  
  return right_hand_side;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
inline array<double, (hyEdge_dim+1) * DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::n_shape_fct_>
DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::
solve_local_problem(const array< array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * n_shape_fct_ * (hyEdge_dim+1) * n_shape_fct_> local_matrix = loc_mat_;
  array<double, (hyEdge_dim+1) * n_shape_fct_> right_hand_side = assemble_rhs(lambda_values);
  
  return lapack_solve<(hyEdge_dim+1) * n_shape_fct_>(local_matrix, right_hand_side);
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
inline array< array<double, DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::n_shape_bdr_> , 2 * hyEdge_dim >
DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::
dual_at_boundary(const array<double, (hyEdge_dim+1) * n_shape_fct_>& coeffs) const
{
  array<unsigned int, hyEdge_dim> dec_i;
  array<unsigned int, max(hyEdge_dim-1,1U)> dec_j;
  double integral, integral1D;
  
  array< array<double, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
  
  for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
  {
    bdr_values[2*dim+0].fill(0.);
    bdr_values[2*dim+1].fill(0.);
  }
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      {
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_; ++q)
              integral1D += q_weights_[q] * trial_[dec_i[dim_fct]][q] * trial_[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        bdr_values[2*dim+0][j] -= coeffs[dim * n_shape_fct_ + i] * integral;
        
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_; ++q)
              integral1D += q_weights_[q] * trial_[dec_i[dim_fct]][q] * trial_[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        bdr_values[2*dim+1][j] += coeffs[dim * n_shape_fct_ + i] * integral;
      }
    }
  }
  
  return bdr_values;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
inline array< array<double, DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::n_shape_bdr_> , 2 * hyEdge_dim > 
DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::
primal_at_boundary(const array<double, (hyEdge_dim+1) * n_shape_fct_>& coeffs) const
{
  array<unsigned int, hyEdge_dim> dec_i;
  array<unsigned int, max(hyEdge_dim-1,1U)> dec_j;
  double integral, integral1D;
  
  array< array<double, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
  
  for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
  {
    bdr_values[2*dim+0].fill(0.);
    bdr_values[2*dim+1].fill(0.);
  }
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
    {
      index_decompose<hyEdge_dim - 1>(j, poly_deg+1, dec_j);
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      {
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][0];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_; ++q)
              integral1D += q_weights_[q] * trial_[dec_i[dim_fct]][q] * trial_[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        bdr_values[2*dim+0][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
        
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_; ++q)
              integral1D += q_weights_[q] * trial_[dec_i[dim_fct]][q] * trial_[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        bdr_values[2*dim+1][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
      }
    }
  }
  
  return bdr_values;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
array< array<double, DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::n_shape_bdr_> , 2 * hyEdge_dim >
DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::
numerical_flux_from_lambda(const array< array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * n_shape_fct_> coeffs = solve_local_problem(lambda_values);
  
  array< array<double, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
  array< array<double, n_shape_bdr_> , 2 * hyEdge_dim > primals = primal_at_boundary(coeffs);
  array< array<double, n_shape_bdr_> , 2 * hyEdge_dim > duals   = dual_at_boundary(coeffs);
  
  for (unsigned int i = 0; i < lambda_values.size(); ++i)
    for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
      bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];
       
  return bdr_values;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
template<unsigned int sizeT>
array<double, Hypercube<hyEdge_dim>::pow(sizeT)>
DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::primal_at_dyadic
(const array<double, sizeT>& abscissas, const array< array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * n_shape_fct_> coefficients = solve_local_problem(lambda_values);
  array<double, Hypercube<hyEdge_dim>::pow(sizeT)> values;
  array< array<double, abscissas.size()>, poly_deg+1 > values1D;
  double fct_value;
  
  values.fill(0.);
  array<unsigned int, hyEdge_dim> dec_i, dec_q;
  
  for (unsigned int deg = 0; deg < poly_deg+1; ++deg)
  {
    for (unsigned int q = 0; q < abscissas.size(); ++q)
      values1D[deg][q] = FuncQuad::shape_fct_eval(deg, abscissas[q]);
  }
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dim>::pow(sizeT); ++q)
    {
      index_decompose<hyEdge_dim>(q, abscissas.size(), dec_q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      values[q] += coefficients[hyEdge_dim * n_shape_fct_ + i] * fct_value;
    }
  }

  return values;
}


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
template<unsigned int sizeT>
array< array<double,hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) >
DiffusionSolverTensorStruc<hyEdge_dim, poly_deg, quad_deg>::dual_at_dyadic
(const array<double, sizeT>& abscissas, const array< array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * n_shape_fct_> coefficients = solve_local_problem(lambda_values);
  array< array<double, hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > values;
  array< array<double, abscissas.size()> , poly_deg+1 > values1D;
  double fct_value;
  
  for (unsigned int i = 0; i < values.size(); ++i)  values[i].fill(0.);
  array<unsigned int, hyEdge_dim> dec_i, dec_q;
  
  for (unsigned int deg = 0; deg < poly_deg+1; ++deg)
  { 
    for (unsigned int q = 0; q < abscissas.size(); ++q)
      values1D[deg][q] = FuncQuad::shape_fct_eval(deg, abscissas[q]);
  }
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  { 
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int q = 0; q < Hypercube<hyEdge_dim>::pow(sizeT); ++q)
    {
      index_decompose<hyEdge_dim>(q, abscissas.size(), dec_q);
      fct_value = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        fct_value *= values1D[dec_i[dim_fct]][dec_q[dim_fct]];
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
        values[q][dim] += coefficients[dim * n_shape_fct_ + i] * fct_value;
    }
  }
  
  return values;
}
