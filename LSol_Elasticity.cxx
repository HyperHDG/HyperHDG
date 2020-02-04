/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include <LSol_Elasticity.hxx>
#include <HyAssert.hxx>
#include <LapackWrapper.hxx>
#include <cmath>

#include <Geom_File.hxx>

using namespace std;
using namespace FuncQuad;
using namespace Geometry;
#include "LSol_Elasticity.inst"


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


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
ElasticRods<poly_deg, quad_deg, GeomT>::
ElasticRods(const constructor_value_type& tau)
: tau_(tau), q_weights_(quad_weights<quad_deg>()),
  trial_(shape_fcts_at_quad_points<poly_deg, quad_deg>()),
  trial_bdr_(shape_fcts_at_bdrs<poly_deg>()),
  loc_mat_(assemble_loc_matrix<hyEdge_dim,poly_deg,quad_deg, n_glob_dofs_per_node()>(tau))
{ } 


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
array<double, (ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim+1) * ElasticRods<poly_deg, quad_deg, GeomT>::n_shape_fct_>
ElasticRods<poly_deg, quad_deg, GeomT>::
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


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
inline array<double, (ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim+1) * ElasticRods<poly_deg, quad_deg, GeomT>::n_shape_fct_>
ElasticRods<poly_deg, quad_deg, GeomT>::
solve_local_problem(const array< array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
{
  array<double, (hyEdge_dim+1) * n_shape_fct_ * (hyEdge_dim+1) * n_shape_fct_> local_matrix = loc_mat_;
  array<double, (hyEdge_dim+1) * n_shape_fct_> right_hand_side = assemble_rhs(lambda_values);
  
  return lapack_solve<(hyEdge_dim+1) * n_shape_fct_>(local_matrix, right_hand_side);
}


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
inline array< array<double, ElasticRods<poly_deg, quad_deg, GeomT>::n_shape_bdr_> , 2 * ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim >
ElasticRods<poly_deg, quad_deg, GeomT>::
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


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
inline array< array<double, ElasticRods<poly_deg, quad_deg, GeomT>::n_shape_bdr_> , 2 * ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim > 
ElasticRods<poly_deg, quad_deg, GeomT>::
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


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
inline array< array<double, ElasticRods<poly_deg, quad_deg, GeomT>::n_shape_bdr_>, 2 * ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim >
ElasticRods<poly_deg, quad_deg, GeomT>::
node_dof_to_edge_dof(const array< array<double, GeomT::space_dim() * n_shape_bdr_>, 2 * ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim > lambda, const GeomT& geom) const
{
  array< array<double, n_shape_bdr_> , 2*hyEdge_dim > result;
  hy_assert( result.size() == 2 , "Only implemented in one dimension!" );
  for (unsigned int i = 0; i < result.size(); ++i)
  {
    hy_assert( result[i].size() == 1 , "Only implemented in one dimension!" );
    result[i].fill(0.);
  }
  
  for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
  {
    Point<space_dim> normal_vector = geom.normal(1);
    for (unsigned int dim = 0; dim < space_dim; ++dim)  result[i][0] += normal_vector[dim] * lambda[i][dim];
  }
  
  return result;
}


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
inline array< array<double, GeomT::space_dim() * ElasticRods<poly_deg, quad_deg, GeomT>::n_shape_bdr_>, 2 * ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim >
ElasticRods<poly_deg, quad_deg, GeomT>::
edge_dof_to_node_dof(const array< array<double, n_shape_bdr_>, 2 * hyEdge_dim > lambda, const GeomT& geom) const
{
  std::array< std::array<double, space_dim * n_shape_bdr_> , 2*hyEdge_dim > result;
  for (unsigned int i = 0; i < result.size(); ++i)  result[i].fill(0.);
  
  for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
  {
    Point<space_dim> normal_vector = geom.normal(1);
    for (unsigned int dim = 0; dim < space_dim; ++dim)  result[i][dim] += normal_vector[dim] * lambda[i][0];
  }
  
  return result;
}


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
array< array<double, GeomT::space_dim() * ElasticRods<poly_deg, quad_deg, GeomT>::n_shape_bdr_> , 2 * ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim >
ElasticRods<poly_deg, quad_deg, GeomT>::
numerical_flux_from_lambda(const array< array<double, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >& lambda_values, const GeomT& geom) const
{
  auto lambda = node_dof_to_edge_dof(lambda_values, geom);
  array<double, (hyEdge_dim+1) * n_shape_fct_> coeffs = solve_local_problem(lambda);
  
  array< array<double, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
  array< array<double, n_shape_bdr_> , 2 * hyEdge_dim > primals = primal_at_boundary(coeffs);
  array< array<double, n_shape_bdr_> , 2 * hyEdge_dim > duals   = dual_at_boundary(coeffs);
  
  for (unsigned int i = 0; i < lambda_values.size(); ++i)
    for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
      bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];
       
  return edge_dof_to_node_dof(bdr_values, geom);
}


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
template<unsigned int sizeT>
array<double, Hypercube<ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim>::pow(sizeT)>
ElasticRods<poly_deg, quad_deg, GeomT>::primal_at_dyadic
(const array<double, sizeT>& abscissas, const array< array<double, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >& lambda_values, const GeomT& geom) const
{
  auto lambda = node_dof_to_edge_dof(lambda_values, geom);
  array<double, (hyEdge_dim+1) * n_shape_fct_> coefficients = solve_local_problem(lambda);
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


template<unsigned int poly_deg, unsigned int quad_deg, class GeomT>
template<unsigned int sizeT>
array< array<double,ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim> , Hypercube<ElasticRods<poly_deg, quad_deg, GeomT>::hyEdge_dim>::pow(sizeT) >
ElasticRods<poly_deg, quad_deg, GeomT>::dual_at_dyadic
(const array<double, sizeT>& abscissas, const array< array<double, GeomT::space_dim() * n_shape_bdr_> , 2*hyEdge_dim >& lambda_values, const GeomT& geom) const
{
  auto lambda = node_dof_to_edge_dof(lambda_values, geom);
  array<double, (hyEdge_dim+1) * n_shape_fct_> coefficients = solve_local_problem(lambda);
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
