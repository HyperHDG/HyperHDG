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

#include "HyperEdge.h"
#include <vector>

template<unsigned int hyperedge_dim>
class DiffusionSolver_RegularQuad
{
  private:
    unsigned int max_poly_degree_, num_of_quad_, num_quad_bdr_;
    double tau_;
    std::vector<double> quad_weights_, quad_bdr_;
    std::vector< std::vector<double> > trials_quad_, bound_trials_quad_, trials_bound_1D_;
    std::vector< std::vector< std::vector<double> > > derivs_quad_, trials_bound_;
    
    unsigned int loc_matrix_index(const unsigned int row, const unsigned int column) const;
    std::vector<double> assemble_loc_mat() const;
    std::vector<double> assemble_rhs(const std::vector< std::vector<double> >& lambda_values) const;
    std::vector<double> solve_local_system_of_eq(std::vector<double>& loc_matrix, std::vector<double>& loc_rhs) const;
    std::vector<double> solve_local_problem(const std::vector< std::vector<double> >& lambda_values) const;
    std::vector< std::vector<double> > dual_at_boundary(const std::vector<double>& coeffs) const;
    std::vector< std::vector<double> > primal_at_boundary(const std::vector<double>& coeffs) const;
    std::vector< std::vector<double> > numerical_flux_at_boundary(const std::vector< std::vector<double> >& lambda_values, const std::vector<double>& coeffs) const;
  public:
    DiffusionSolver_RegularQuad(const unsigned int max_poly_degree, const unsigned int num_of_quad, const double tau);
    std::vector< std::vector<double> >numerical_flux_from_lambda(const std::vector< std::vector<double> >& lambda_values) const;
};

#endif
