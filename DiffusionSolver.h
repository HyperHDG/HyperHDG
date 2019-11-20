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


#ifndef DIFFUSIONSOLVER_H
#define DIFFUSIONSOLVER_H

#include "FuncAndQuad.h"
#include "HyperNodeFactory.h"
#include <array>

template<unsigned int hyperedge_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
class DiffusionSolver_RegularQuad
{
  private:
    static constexpr unsigned int num_of_quad_    = quadrature_points_amount(max_quad_degree, hyperedge_dim),
                                  num_quad_bdr_   = quadrature_points_amount(max_quad_degree, hyperedge_dim - 1),
                                  num_ansatz_fct_ = local_dof_amount_node(hyperedge_dim, max_poly_degree) * (max_poly_degree + 1),
                                  num_ansatz_bdr_ = local_dof_amount_node(hyperedge_dim, max_poly_degree);
    const double tau_;
    std::array<double, num_of_quad_> quad_weights_;
    std::array<double, num_quad_bdr_> quad_bdr_;
    std::array< std::array<double, num_of_quad_> , num_ansatz_fct_ > trials_quad_;
    std::array< std::array<double, num_quad_bdr_> , num_ansatz_bdr_ > bound_trials_quad_;
    std::array< std::array<double, 2> , max_poly_degree + 1 > trials_bound_1D_;
    std::array< std::array< std::array<double, num_of_quad_> , num_ansatz_fct_ > , hyperedge_dim > derivs_quad_;
    std::array< std::array< std::array<double, num_quad_bdr_> , num_ansatz_fct_ > , 2 * hyperedge_dim > trials_bound_;
    
    unsigned int loc_matrix_index(const unsigned int row, const unsigned int column) const;
    auto assemble_loc_mat() const; // std::array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_>
    auto assemble_rhs(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const; // std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>
    auto solve_local_system_of_eq // std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>
      (std::array<double, (hyperedge_dim+1) * num_ansatz_fct_ * (hyperedge_dim+1) * num_ansatz_fct_>& loc_matrix, std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>& loc_rhs) const;
    auto solve_local_problem(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const; // std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>
    auto dual_at_boundary(const std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>& coeffs) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
    auto primal_at_boundary(const std::array<double, (hyperedge_dim+1) * num_ansatz_fct_>& coeffs) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
    auto numerical_flux_at_boundary // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
      (const std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >& lambda_values, const std::array<double, (hyperedge_dim + 1) * num_ansatz_fct_>& coeffs) const;
  public:
    DiffusionSolver_RegularQuad(const double tau);
    std::array< std::array<double, local_dof_amount_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim > 
      primal_at_boundary_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
    std::array< std::array<double, local_dof_amount_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim > 
      flux_at_boundary_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
    std::array< std::array<double, local_dof_amount_node(hyperedge_dim, max_poly_degree)> , 2 * hyperedge_dim >
      numerical_flux_from_lambda(const std::array< std::array<double, num_ansatz_bdr_> , 2*hyperedge_dim >& lambda_values) const; // std::array< std::array<double, num_ansatz_bdr_> , 2 * hyperedge_dim >
};

#endif
