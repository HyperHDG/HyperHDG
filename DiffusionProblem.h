/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef DIFFUSIONPROBLEM_H
#define DIFFUSIONPROBLEM_H

#include "HyperGraph.h"
#include "DiffusionSolver.h"

/**
 * This is an example problem.
 *
 * \todo The loop in matrix_vector_multiply() only combines properties of HyperGraph with local solvers, right? Dirichlet boundary conditions? Post filtering!
 */
template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
class DiffusionProblemRegular
{
  private:
    HyperGraph < local_dof_amount_node(hyperedge_dim, polynomial_degree),
                 HyperGraph_Cubic< hyperedge_dim, space_dim >,
                 HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >
               hyper_graph_;
    std::vector<unsigned int> dirichlet_indices_;
    DiffusionSolver_RegularQuad<hyperedge_dim, polynomial_degree, 2 * polynomial_degree> local_solver_;
  public:
    DiffusionProblemRegular(std::vector<int> num_elements);
    void read_dirichlet_indices(std::vector<int> indices);
    std::vector<double> return_zero_vector();
    std::vector<double> matrix_vector_multiply(std::vector<double> x_vec);
    int size_of_system();
    void plot_solution(std::vector<double> lambda);
};

#endif
