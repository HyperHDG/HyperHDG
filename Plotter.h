/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 * 
 * Definition of a class for point objects providing getter functions only.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef PLOTTER_H
#define PLOTTER_H

#include "HyperGraph.h"
#include "DiffusionSolver.h"

template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
class Plotter
{
  private:
    HyperGraph < local_dof_amount_node(hyperedge_dim, polynomial_degree),
                 HyperGraph_Cubic< hyperedge_dim, space_dim >,
                 HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >&
               hyper_graph_;
    DiffusionSolver_RegularQuad<hyperedge_dim, polynomial_degree, 2 * polynomial_degree>& local_solver_;
  public:
    Plotter(HyperGraph< local_dof_amount_node(hyperedge_dim, polynomial_degree), HyperGraph_Cubic< hyperedge_dim, space_dim >, HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >& hyper_graph,
            DiffusionSolver_RegularQuad<hyperedge_dim, polynomial_degree, 2 * polynomial_degree>& local_solver);
    void plot(std::vector<double> lambda);
};

#endif
