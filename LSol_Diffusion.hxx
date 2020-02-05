/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 * 
 * Definition local solver class: In this case for diffusion.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef LSOL_DIFFUSION_HXX
#define LSOL_DIFFUSION_HXX

#include <FuncAndQuad.hxx>
#include <Hypercube.hxx>
#include <LapackWrapper.hxx>

#include <array>
#include <cmath>
/*
template<unsigned int hyEdge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
class DiffusionSolverNaive_RegularQuad
{
  public:
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    static constexpr bool use_geometry() { return false; }
    static constexpr unsigned int n_glob_dofs_per_node()
    {
      unsigned int amount = 1;
      for (unsigned int iteration = 0; iteration < hyEdge_dim - 1; ++ iteration)  amount *= max_poly_degree + 1;
      return amount;
    }
  private:
    static constexpr unsigned int n_quads_1D_        = FuncQuad::compute_n_quad_points(max_quad_degree, hyEdge_dim),
                                  num_quad_bdr_   = FuncQuad::compute_n_quad_points(max_quad_degree, hyEdge_dim - 1),
                                  num_ansatz_fct_ = n_glob_dofs_per_node() * (max_poly_degree + 1),
                                  num_ansatz_bdr_ = n_glob_dofs_per_node();
    const double tau_;
    std::array<double, n_quads_1D_> quad_weights_;
    std::array<double, num_quad_bdr_> quad_bdr_;
    std::array< std::array<double, n_quads_1D_> , num_ansatz_fct_ > trials_quad_;
    std::array< std::array<double, num_quad_bdr_> , num_ansatz_bdr_ > bound_trials_quad_;
    std::array< std::array<double, (1 << hyEdge_dim)> , num_ansatz_fct_ > trials_in_corners_;
    std::array< std::array< std::array<double, n_quads_1D_> , num_ansatz_fct_ > , hyEdge_dim > derivs_quad_;
    std::array< std::array< std::array<double, num_quad_bdr_> , num_ansatz_fct_ > , 2 * hyEdge_dim > trials_bound_;
     

    inline auto assemble_loc_mat() const; // std::array<double, (hyEdge_dim+1) * num_ansatz_fct_ * (hyEdge_dim+1) * num_ansatz_fct_>
    inline auto assemble_rhs(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const; // std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>
//    auto solve_local_system_of_eq // std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>
//      (std::array<double, (hyEdge_dim+1) * num_ansatz_fct_ * (hyEdge_dim+1) * num_ansatz_fct_>& loc_matrix, std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>& loc_rhs) const;
    inline auto solve_local_problem(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const; // std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>
    inline auto dual_at_boundary(const std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>& coeffs) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
    inline auto primal_at_boundary(const std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>& coeffs) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
    auto numerical_flux_at_boundary // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
      (const std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >& lambda_values, const std::array<double, (hyEdge_dim + 1) * num_ansatz_fct_>& coeffs) const;
  public:
    typedef double constructor_value_type;
    DiffusionSolverNaive_RegularQuad(const constructor_value_type& tau);
    std::array<double, Hypercube<hyEdge_dim>::n_vertices()> primal_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    std::array< std::array<double, hyEdge_dim> , Hypercube<hyEdge_dim>::n_vertices() > dual_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    auto//std::array< std::array<double, n_glob_dofs_per_node()> , 2 * hyEdge_dim >
      numerical_flux_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim 
};
*/


/*!
 * \todo  Will diffusivity depend on spatial (local or global) coordinates and/or the index of the hyperedge? How should this be encoded?
 */

template<unsigned int hyEdge_dim, unsigned int poly_deg>
inline unsigned int loc_matrix_index(const unsigned int row, const unsigned int column)
{
  hy_assert( 0 <= row ,
             "Row index should be larger than or equal to zero." );
  hy_assert( row < (hyEdge_dim + 1) * std::pow((poly_deg + 1), hyEdge_dim) ,
             "Row index should be smaller than total amount of rows." );
  hy_assert( 0 <= column ,
             "Column index should be larger than or equal to zero." );
  hy_assert( column < (hyEdge_dim + 1) * std::pow((poly_deg + 1), hyEdge_dim) ,
             "Column index should smaller than total amount of columns." );
  return column * (hyEdge_dim + 1) * std::pow((poly_deg + 1), hyEdge_dim) + row;  // Transposed for LAPACK
}

template <unsigned int hyEdge_dim>
inline void index_decompose(unsigned int index, unsigned int range, std::array<unsigned int, std::max(hyEdge_dim,1U)>& decomposition)
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


template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
class DiffusionSolverTensorStruc
{
  public:
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    static constexpr bool use_geometry() { return false; }
    static constexpr unsigned int n_glob_dofs_per_node()
    { return Hypercube<hyEdge_dim-1>::pow(poly_deg + 1); }
  private:
    static constexpr unsigned int n_quads_1D_  = FuncQuad::compute_n_quad_points(quad_deg);
    static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1);
    static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node();
    static constexpr unsigned int n_loc_dofs_  = (hyEdge_dim+1) * n_shape_fct_;


static std::array<double, n_loc_dofs_ * n_loc_dofs_ > assemble_loc_matrix(const double tau)
{ 
  const std::array<double, n_quads_1D_> q_weights = FuncQuad::quad_weights<quad_deg>();
  const std::array< std::array<double, n_quads_1D_ > , poly_deg + 1 > trial = FuncQuad::shape_fcts_at_quad_points<poly_deg, quad_deg>();
  const std::array< std::array<double, n_quads_1D_ > , poly_deg + 1 > deriv = FuncQuad::shape_ders_at_quad_points<poly_deg, quad_deg>();
  const std::array< std::array<double, 2> , poly_deg + 1 > trial_bdr = FuncQuad::shape_fcts_at_bdrs<poly_deg>();
  
  std::array<unsigned int, hyEdge_dim> dec_i, dec_j;
  double integral, integral1D;
  
  std::array<double, n_loc_dofs_ * n_loc_dofs_> local_mat;
  local_mat.fill(0.);
  
  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    index_decompose<hyEdge_dim>(i, poly_deg+1, dec_i);
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      index_decompose<hyEdge_dim>(j, poly_deg+1, dec_j);
      
      // Integral_element phi_i phi_j dx in diagonal blocks
      integral = 1.;
      for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
      {
        integral1D = 0.;
        for (unsigned int q = 0; q < n_quads_1D_; ++q)
          integral1D += q_weights[q] * trial[dec_j[dim_fct]][q] * trial[dec_i[dim_fct]][q];
        integral *= integral1D;
      }
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( dim * n_shape_fct_ + i , dim * n_shape_fct_ + j )] += integral;
      
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      { 
        // Integral_element - nabla phi_i \vec phi_j dx = Integral_element - div \vec phi_i phi_j dx
        // in right upper and left lower blocks
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          integral1D = 0.;
          for (unsigned int q = 0; q < n_quads_1D_; ++q)
            integral1D += q_weights[q] * trial[dec_j[dim_fct]][q] * 
                          ( ( dim == dim_fct ) ? deriv[dec_i[dim_fct]][q] : trial[dec_i[dim_fct]][q] );
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct_ + i , dim * n_shape_fct_ + j )] -= integral;
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>(  dim * n_shape_fct_ + i , hyEdge_dim * n_shape_fct_ + j )] -= integral;
        
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr[dec_i[dim_fct]][1] * trial_bdr[dec_j[dim_fct]][1];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_1D_; ++q)
              integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct_ + i , dim * n_shape_fct_ + j )] += integral;
        // Corresponding boundary integrals from integration by parts in left lower blocks
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr[dec_i[dim_fct]][0] * trial_bdr[dec_j[dim_fct]][0];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_1D_; ++q)
              integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct_ + i , dim * n_shape_fct_ + j )] -= integral;
        
        // Penalty in lower right diagonal block
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr[dec_i[dim_fct]][0] * trial_bdr[dec_j[dim_fct]][0];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_1D_; ++q)
              integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct_ + i , hyEdge_dim * n_shape_fct_ + j )] += tau * integral;
        // Penalty in lower right diagonal block
        integral = 1.;
        for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
        {
          if (dim == dim_fct)  integral1D = trial_bdr[dec_i[dim_fct]][1] * trial_bdr[dec_j[dim_fct]][1];
          else
          {
            integral1D = 0.;
            for (unsigned int q = 0; q < n_quads_1D_; ++q)
              integral1D += q_weights[q] * trial[dec_i[dim_fct]][q] * trial[dec_j[dim_fct]][q];
          }
          integral *= integral1D;
        }
        local_mat[loc_matrix_index<hyEdge_dim,poly_deg>( hyEdge_dim * n_shape_fct_ + i , hyEdge_dim * n_shape_fct_ + j )] += tau * integral;
      }
    }
  }
  
  return local_mat;
}

    const double tau_;
    const std::array<double, n_quads_1D_> q_weights_;
    const std::array< std::array<double, n_quads_1D_ > , poly_deg + 1 > trial_;
    const std::array< std::array<double, 2 > , poly_deg + 1 > trial_bdr_;
    const std::array<double, n_loc_dofs_ * n_loc_dofs_ > loc_mat_;
    
    inline std::array<double, n_loc_dofs_ > assemble_rhs
    ( const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values ) const
    {
      std::array<unsigned int, hyEdge_dim> dec_i;
      std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
      double integral, integral1D;
  
      std::array<double, (hyEdge_dim+1) * n_shape_fct_> right_hand_side;
      right_hand_side.fill(0.);
  
      hy_assert( lambda_values.size() == 2 * hyEdge_dim ,
                  "The size of the lambda values should be twice the dimension of a hyperedge." );
      for (unsigned int i = 0; i < 2 * hyEdge_dim; ++i)
        hy_assert( lambda_values[i].size() == n_shape_bdr_ ,
                   "The size of lambda should be the amount of ansatz functions at boundary." );
  
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
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            right_hand_side[dim * n_shape_fct_ + i] += lambda_values[2*dim+0][j] * integral;
            right_hand_side[hyEdge_dim*n_shape_fct_ + i] += tau_*lambda_values[2*dim+0][j]*integral;
        
            integral = 1.;
            for (unsigned int dim_fct = 0; dim_fct < hyEdge_dim; ++dim_fct)
            {
              if (dim == dim_fct)  integral1D = trial_bdr_[dec_i[dim_fct]][1];
              else
              {
                integral1D = 0.;
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q] * trial_[dec_i[dim_fct]][q] * trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            right_hand_side[dim*n_shape_fct_ + i] -= lambda_values[2*dim+1][j] * integral;
            right_hand_side[hyEdge_dim*n_shape_fct_ + i] += tau_*lambda_values[2*dim+1][j]*integral;
          }
        }
      }
  
      return right_hand_side;
    }


    inline std::array<double, n_loc_dofs_ > solve_local_problem
    ( const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values ) const
    {
      std::array<double, n_loc_dofs_ * n_loc_dofs_> local_matrix = loc_mat_;
      std::array<double, n_loc_dofs_> right_hand_side = assemble_rhs(lambda_values);
  
      return lapack_solve<n_loc_dofs_>(local_matrix, right_hand_side);
    }

    inline std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > primal_at_boundary
    ( const std::array<double, n_loc_dofs_ >& coeffs ) const
    {
      std::array<unsigned int, hyEdge_dim> dec_i;
      std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
      std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
      double integral, integral1D;
    
      for (unsigned int dim_n = 0; dim_n < 2 * hyEdge_dim; ++dim_n)  bdr_values[dim_n].fill(0.);
  
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
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
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
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            bdr_values[2*dim+1][j] += coeffs[hyEdge_dim * n_shape_fct_ + i] * integral;
          }
        }
      }
      
      return bdr_values;
    }

    inline std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > dual_at_boundary
    ( const std::array<double, (hyEdge_dim+1) * n_shape_fct_>& coeffs ) const
    {
      std::array<unsigned int, hyEdge_dim> dec_i;
      std::array<unsigned int, std::max(hyEdge_dim-1,1U)> dec_j;
      std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > bdr_values;
      double integral, integral1D;
    
      for (unsigned int dim_n = 0; dim_n < 2*hyEdge_dim; ++dim_n)  bdr_values[dim_n].fill(0.);

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
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
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
                for (unsigned int q = 0; q < n_quads_1D_; ++q)
                  integral1D += q_weights_[q]*trial_[dec_i[dim_fct]][q]*trial_[dec_j[dim_fct]][q];
              }
              integral *= integral1D;
            }
            bdr_values[2*dim+1][j] += coeffs[dim * n_shape_fct_ + i] * integral;
          }
        }
      }
  
      return bdr_values;
    }

  public:
    typedef double constructor_value_type;

    DiffusionSolverTensorStruc(const constructor_value_type& tau)
    : tau_(tau), q_weights_(FuncQuad::quad_weights<quad_deg>()),
      trial_(FuncQuad::shape_fcts_at_quad_points<poly_deg, quad_deg>()),
      trial_bdr_(FuncQuad::shape_fcts_at_bdrs<poly_deg>()),
      loc_mat_(assemble_loc_matrix(tau))
    { } 
    
    // Function for matrix--vector multiply    
    std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > numerical_flux_from_lambda
    (const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const
    {
      std::array<double, n_loc_dofs_ > coeffs = solve_local_problem(lambda_values);
      
      std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > 
        bdr_values, primals(primal_at_boundary(coeffs)), duals(dual_at_boundary(coeffs));
  
      for (unsigned int i = 0; i < lambda_values.size(); ++i)
        for (unsigned int j = 0; j < lambda_values[i].size(); ++j)
          bdr_values[i][j] = duals[i][j] + tau_ * primals[i][j] - tau_ * lambda_values[i][j];
       
      return bdr_values;
    }
    
    // Plotting functions
    template<unsigned int sizeT> std::array<double, Hypercube<hyEdge_dim>::pow(sizeT)>
    primal_at_dyadic
    ( const std::array<double, sizeT>& abscissas,
      const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values ) const
    {
      std::array<double, n_loc_dofs_ > coefficients = solve_local_problem(lambda_values);
      std::array<double, Hypercube<hyEdge_dim>::pow(sizeT)> values;
      std::array< std::array<double, abscissas.size()>, poly_deg+1 > values1D;
      std::array<unsigned int, hyEdge_dim> dec_i, dec_q;
      double fct_value;
  
      for (unsigned int deg = 0; deg < poly_deg+1; ++deg)
        for (unsigned int q = 0; q < abscissas.size(); ++q)
          values1D[deg][q] = shape_fct_eval(deg, abscissas[q]);
      
      values.fill(0.);
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
    /*!*********************************************************************************************
     * \brief   Evaluate discrete function at given points.
     *
     * Function to evaluate dual variable / part of the solution at dyadic product of abscissas. The
     * main purpose of the function in closely related to plotting.
     *
     * \tparam  sizeT         Size of the passed \c std::array containing the abscissas.
     * \param   abscissas     Coordinates at whose tensor products the function is evaluated.
     * \param   lambda_values Coefficients of the associated skeletal function.
     * \retval  fct_val       Evaluation of dual variable at prescribed points.
     **********************************************************************************************/
    template<unsigned int sizeT>
    std::array< std::array<double,hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > dual_at_dyadic
    ( const std::array<double, sizeT>& abscissas,
      const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values ) const
    {
      std::array<double, n_loc_dofs_> coefficients = solve_local_problem(lambda_values);
      std::array< std::array<double, hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > values;
      std::array< std::array<double, abscissas.size()> , poly_deg+1 > values1D;
      std::array<unsigned int, hyEdge_dim> dec_i, dec_q;
      double fct_value;
      
      for (unsigned int deg = 0; deg < poly_deg+1; ++deg)
        for (unsigned int q = 0; q < abscissas.size(); ++q)
          values1D[deg][q] = FuncQuad::shape_fct_eval(deg, abscissas[q]);

      for (unsigned int i = 0; i < values.size(); ++i)  values[i].fill(0.);
  
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
}; // end of class DiffusionSolverTensorStruc

#endif