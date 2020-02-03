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

#include <array>

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
    static constexpr unsigned int n_quads_        = FuncQuad::compute_n_quad_points(max_quad_degree, hyEdge_dim),
                                  num_quad_bdr_   = FuncQuad::compute_n_quad_points(max_quad_degree, hyEdge_dim - 1),
                                  num_ansatz_fct_ = n_glob_dofs_per_node() * (max_poly_degree + 1),
                                  num_ansatz_bdr_ = n_glob_dofs_per_node();
    const double tau_;
    std::array<double, n_quads_> quad_weights_;
    std::array<double, num_quad_bdr_> quad_bdr_;
    std::array< std::array<double, n_quads_> , num_ansatz_fct_ > trials_quad_;
    std::array< std::array<double, num_quad_bdr_> , num_ansatz_bdr_ > bound_trials_quad_;
    std::array< std::array<double, (1 << hyEdge_dim)> , num_ansatz_fct_ > trials_in_corners_;
    std::array< std::array< std::array<double, n_quads_> , num_ansatz_fct_ > , hyEdge_dim > derivs_quad_;
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

template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
class DiffusionSolverTensorStruc
{
  public:
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    static constexpr bool use_geometry() { return false; }
    static constexpr unsigned int n_glob_dofs_per_node()
    {
      unsigned int amount = 1;
      for (unsigned int iteration = 0; iteration < hyEdge_dim - 1; ++ iteration)  amount *= poly_deg + 1;
      return amount;
    }
  private:
    static constexpr unsigned int n_quads_     = FuncQuad::compute_n_quad_points(quad_deg),
                                  n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1),
                                  n_shape_bdr_ = n_glob_dofs_per_node();
    const double tau_;
    const std::array<double, n_quads_> q_weights_;
    const std::array< std::array<double, n_quads_ > , poly_deg + 1 > trial_;
    const std::array< std::array<double, 2 > , poly_deg + 1 > trial_bdr_;
    const std::array<double, (hyEdge_dim+1) * n_shape_fct_ * (hyEdge_dim+1) * n_shape_fct_> loc_mat_;
    
    inline std::array<double, (hyEdge_dim+1) * n_shape_fct_> assemble_rhs(const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    inline std::array<double, (hyEdge_dim+1) * n_shape_fct_> solve_local_problem(const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    
    inline std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > dual_at_boundary(const std::array<double, (hyEdge_dim+1) * n_shape_fct_>& coeffs) const;
    inline std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim > primal_at_boundary(const std::array<double, (hyEdge_dim+1) * n_shape_fct_>& coeffs) const;
  public:
    typedef double constructor_value_type;
    DiffusionSolverTensorStruc(const constructor_value_type& tau);
    
    // Function for matrix--vector multiply    
    std::array< std::array<double, n_shape_bdr_> , 2 * hyEdge_dim >
      numerical_flux_from_lambda(const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    
    // Plotting functions
    template<unsigned int sizeT> std::array<double, Hypercube<hyEdge_dim>::pow(sizeT)> primal_at_dyadic
      (const std::array<double, sizeT>& abscissas, const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    template<unsigned int sizeT> std::array< std::array<double,hyEdge_dim> , Hypercube<hyEdge_dim>::pow(sizeT) > dual_at_dyadic
      (const std::array<double, sizeT>& abscissas, const std::array< std::array<double, n_shape_bdr_> , 2*hyEdge_dim >& lambda_values) const;
};

#endif
