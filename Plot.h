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


#ifndef PLOT_H
#define PLOT_H

#include "HDGHyperGraph.h"
#include "DiffusionSolver.h"

/**
 * \todo I do not see the point of having this as an object. A
 * function template should do it. Or could we output several datasets
 * on the same mesh?
 *
 * \todo It must be possible to give it a file handle or file name
 *
 * \todo The name should indicate VTU
 */

template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
class PlotOptions
{
  private:
    const HDGHyperGraph < compute_n_dofs_per_node(hyperedge_dim, polynomial_degree),
                          Topology::HyperGraph_Cubic< hyperedge_dim, space_dim >,
                          Geometry::HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >&
                        hyper_graph_;
    const DiffusionSolver_RegularQuad<hyperedge_dim, polynomial_degree, 2 * polynomial_degree>& local_solver_;
  public:
    std::string outputDir, fileName, fileEnding;
    unsigned int fileNumber;
    bool printFileNumber, incrementFileNumber;
    
    PlotOptions(HDGHyperGraph< compute_n_dofs_per_node(hyperedge_dim, polynomial_degree), Topology::HyperGraph_Cubic< hyperedge_dim, space_dim >, Geometry::HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >& hyper_graph,
            DiffusionSolver_RegularQuad<hyperedge_dim, polynomial_degree, 2 * polynomial_degree>& local_solver);
    
    const HDGHyperGraph < compute_n_dofs_per_node(hyperedge_dim, polynomial_degree),
                          Topology::HyperGraph_Cubic< hyperedge_dim, space_dim >,
                          Geometry::HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >&
                        hyper_graph();
    const DiffusionSolver_RegularQuad<hyperedge_dim, polynomial_degree, 2 * polynomial_degree>& local_solver();
};


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
void plot(std::vector<double> lambda, PlotOptions<hyperedge_dim,space_dim,polynomial_degree>& plotOpt);

#endif
