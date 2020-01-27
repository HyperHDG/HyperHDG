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
constexpr const unsigned int compute_n_corners_of_cube(const unsigned int hyperedge_dim)
{
  unsigned int amount = 1;
  for (unsigned int dim = 0; dim < hyperedge_dim; ++dim) amount *= 2;
  return amount;
}
#endif

#ifndef ELASTICITYSOLVER_H
#define ELASTICITYSOLVER_H

#include "FuncAndQuad.h"
#include "HyperEdge_Geometry.h"
#include "HyperNodeFactory.h"
#include <array>

template<unsigned int hyperedge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
class ElasticitySolver_RegularQuad
{
  private:
    static constexpr unsigned int num_of_quad_    = FuncQuad::compute_n_quad_points(max_quad_degree, hyperedge_dim),
                                  num_quad_bdr_   = FuncQuad::compute_n_quad_points(max_quad_degree, hyperedge_dim - 1),
                                  num_ansatz_fct_ = compute_n_dofs_per_node(hyperedge_dim, max_poly_degree) * (max_poly_degree + 1),
                                  num_ansatz_bdr_ = compute_n_dofs_per_node(hyperedge_dim, max_poly_degree);
    const double tau_;
    std::array<double, num_of_quad_> quad_weights_;
    std::array<double, num_quad_bdr_> quad_bdr_;
    std::array< std::array<double, num_of_quad_> , num_ansatz_fct_ > trials_quad_;
    std::array< std::array<double, num_quad_bdr_> , num_ansatz_bdr_ > bound_trials_quad_;
    std::array< std::array<double, compute_n_corners_of_cube(hyperedge_dim)> , max_poly_degree + 1 > trials_in_corners_;
    std::array< std::array< std::array<double, num_of_quad_> , num_ansatz_fct_ > , hyperedge_dim > derivs_quad_;
    std::array< std::array< std::array<double, num_quad_bdr_> , num_ansatz_fct_ > , 2 * hyperedge_dim > trials_bound_;
     
    inline unsigned int loc_matrix_index(const unsigned int row, const unsigned int column) const;
    inline auto assemble_loc_mat() const; // std::array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_>
    inline auto assemble_rhs(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const; // std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>
    auto solve_local_system_of_eq // std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>
      (std::array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_>& loc_matrix, std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>& loc_rhs) const;
    inline auto solve_local_problem(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const; // std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>
    inline auto dual_at_boundary(const std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>& coeffs) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
    inline auto primal_at_boundary(const std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>& coeffs) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
    auto numerical_flux_at_boundary // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
      (const std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >& lambda_values, const std::array<double, (hyperedge_dim + 1) * num_ansatz_fct_>& coeffs) const;
  public:
    typedef double constructor_value_type;
    ElasticitySolver_RegularQuad(const constructor_value_type& tau);
    std::array<double, compute_n_corners_of_cube(hyperedge_dim)> primal_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const;
    std::array< std::array<double, hyperedge_dim> , compute_n_corners_of_cube(hyperedge_dim) > dual_in_corners_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const;
    std::array< std::array<double, compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim >
      numerical_flux_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
    
    std::array< std::array<double, compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim >
      preprocess_data(std::array< std::array<double, space_dim * compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2*hyperedge_dim >& hyperedge_dofs, Geometry::HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim>& geometry ) const;
    std::array< std::array<double, space_dim * compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim >
      postprocess_data(std::array< std::array<double, compute_n_dofs_per_node(hyperedge_dim, max_poly_degree)> , 2*hyperedge_dim >& hyperedge_dofs, Geometry::HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim>& geometry ) const;
    
    static constexpr unsigned int hyperedge_dimension() { return hyperedge_dim; };
    static constexpr unsigned int polynomial_degree() { return max_poly_degree; };
    static constexpr unsigned int solution_dimension_hyperedge() { return hyperedge_dim; };
    static constexpr unsigned int solution_dimension_hypernode() { return space_dim; };
    static constexpr bool need_geometry_processing() { return true; };
};

#endif
