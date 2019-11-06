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

#include "HyperGraphTopology.h"
#include "DiffusionSolver.h"

template <unsigned int connector_dim, unsigned int space_dim>
class DiffusionProblemRegular
{
  private:
    HyperGraphTopology < JointGetter_RegularQuad<connector_dim,space_dim>,
                         ConnectorGetter_RegularQuad<connector_dim,space_dim>,
                         Joint_RegularQuad,
                         Connector_RegularQuad<connector_dim,space_dim> >
                       hyper_graph_topology;
    std::vector<int> dirichlet_indices;
    DiffusionSolver_RegularQuad<connector_dim> local_solver;
    std::vector<unsigned int> Dirichlet_indices;
  public:
    DiffusionProblemRegular(std::vector<int> num_elements, int polynomials);
    std::vector<double> return_zero_vector();
    std::vector<double> matrix_vector_multiply(std::vector<double> x_vec);
    int size_of_system();
};

#endif
