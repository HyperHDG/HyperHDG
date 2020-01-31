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


#ifndef COMPUTE_CORNERS
#define COMPUTE_CORNERS
// Naive implementation finding the amount of corners for a hypersquare
constexpr const unsigned int compute_n_corners_of_cube(const unsigned int hyEdge_dim)
{
  unsigned int amount = 1;
  for (unsigned int dim = 0; dim < hyEdge_dim; ++dim) amount *= 2;
  return amount;
}
#endif


#ifndef DIFFUSIONSOLVERNAIVE_HXX
#define DIFFUSIONSOLVERNAIVE_HXX

#include "FuncAndQuad.hxx"
#include "HyperNodeFactory.hxx"
#include <array>
#include <vector>


template<unsigned int hyEdge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
class DiffusionSolverNaive_RegularQuad
{
  private:
    static constexpr unsigned int n_quads_    = FuncQuad::compute_n_quad_points(max_quad_degree, hyEdge_dim),
                                  num_quad_bdr_   = FuncQuad::compute_n_quad_points(max_quad_degree, hyEdge_dim - 1),
                                  num_ansatz_fct_ = compute_n_dofs_per_node(hyEdge_dim, max_poly_degree) * (max_poly_degree + 1),
                                  num_ansatz_bdr_ = compute_n_dofs_per_node(hyEdge_dim, max_poly_degree);
    const double tau_;
    std::array<double, n_quads_> quad_weights_;
    std::array<double, num_quad_bdr_> quad_bdr_;
    std::array< std::array<double, n_quads_> , num_ansatz_fct_ > trials_quad_;
    std::array< std::array<double, num_quad_bdr_> , num_ansatz_bdr_ > bound_trials_quad_;
    std::array< std::array<double, compute_n_corners_of_cube(hyEdge_dim)> , num_ansatz_fct_ > trials_in_corners_;
    std::array< std::array< std::array<double, n_quads_> , num_ansatz_fct_ > , hyEdge_dim > derivs_quad_;
    std::array< std::array< std::array<double, num_quad_bdr_> , num_ansatz_fct_ > , 2 * hyEdge_dim > trials_bound_;
     
    inline unsigned int loc_matrix_index(const unsigned int row, const unsigned int column) const;
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
    std::vector<double> primal_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    std::vector< std::array<double, hyEdge_dim> > dual_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    std::array< std::array<double, compute_n_dofs_per_node(hyEdge_dim, max_poly_degree)> , 2 * hyEdge_dim >
      numerical_flux_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    static constexpr unsigned int polynomial_degree() { return max_poly_degree; }
    static constexpr unsigned int solution_dimension_hyEdge() { return 1; }
    static constexpr unsigned int solution_dimension_hyNode() { return 1; }
    static constexpr bool need_geometry_processing() { return false; }
};



template<unsigned int hyEdge_dim, unsigned int poly_deg, unsigned int quad_deg>
class DiffusionSolverTensorStruc
{
  private:
    static constexpr unsigned int n_quads_        = FuncQuad::compute_n_quad_points(quad_deg),
                                  num_ansatz_fct_ = compute_n_dofs_per_node(hyEdge_dim, poly_deg) * (poly_deg + 1),
                                  num_ansatz_bdr_ = compute_n_dofs_per_node(hyEdge_dim, poly_deg);
    const double tau_;
    const std::array<double, n_quads_> q_weights_;
    const std::array< std::array<double, n_quads_ > , poly_deg + 1 > trial_, deriv_;
    const std::array< std::array<double, 2> , poly_deg + 1 > trial_bdr_;
    
    inline unsigned int loc_matrix_index(const unsigned int row, const unsigned int column) const;
    inline auto assemble_loc_mat() const; // std::array<double, (hyEdge_dim+1) * num_ansatz_fct_ * (hyEdge_dim+1) * num_ansatz_fct_>
    inline auto assemble_rhs(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const; // std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>
    inline auto solve_local_problem(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const; // std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>
    inline auto dual_at_boundary(const std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>& coeffs) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
    inline auto primal_at_boundary(const std::array<double, (hyEdge_dim+1) * num_ansatz_fct_>& coeffs) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
    auto numerical_flux_at_boundary // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
      (const std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >& lambda_values, const std::array<double, (hyEdge_dim + 1) * num_ansatz_fct_>& coeffs) const;
  public:
    typedef double constructor_value_type;
    DiffusionSolverTensorStruc(const constructor_value_type& tau);
    
    std::vector<double> primal_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    std::vector< std::array<double, hyEdge_dim> > dual_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    std::array< std::array<double, compute_n_dofs_per_node(hyEdge_dim, poly_deg)> , 2 * hyEdge_dim >
      numerical_flux_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyEdge_dim >
    
    std::vector<double> primal_at_dyadic(const std::vector<double>& abscissas, const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    std::vector< std::array<double,hyEdge_dim> > dual_at_dyadic(const std::vector<double>& abscissas, const std::array< std::array<double, num_ansatz_bdr_> , 2*hyEdge_dim >& lambda_values) const;
    
    static constexpr unsigned int hyEdge_dimension() { return hyEdge_dim; }
    static constexpr unsigned int polynomial_degree() { return poly_deg; }
    static constexpr unsigned int solution_dimension_hyEdge() { return 1; }
    static constexpr unsigned int solution_dimension_hyNode() { return 1; }
    static constexpr bool need_geometry_processing() { return false; }
};

#endif
