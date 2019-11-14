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

template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
class DiffusionProblemRegular
{
  private:
    HyperGraph < local_dof_amount(hyperedge_dim, polynomial_degree),
                 HyperGraph_Cubic< hyperedge_dim, space_dim > >
               hyper_graph_topology;
    std::vector<int> dirichlet_indices;
    DiffusionSolver_RegularQuad<hyperedge_dim> local_solver;
    std::vector<unsigned int> Dirichlet_indices;
  public:
    DiffusionProblemRegular(std::vector<int> num_elements);
    std::vector<double> return_zero_vector();
    std::vector<double> matrix_vector_multiply(std::vector<double> x_vec);
    int size_of_system();
};

#endif
