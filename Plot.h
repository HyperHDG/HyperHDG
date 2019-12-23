#ifndef PLOT_H
#define PLOT_H

#include "HDGHyperGraph.h"
#include "DiffusionSolver.h"

/*!*************************************************************************************************
 * @brief   A class to encode the chosen options for plotting.
 *
 * \todo  I do not see the point of having this as an object. A function template should do it. Or
 *        could we output several datasets on the same mesh?
 *        -> This has been changed and a PlotOptions class has been introduced, where the output
 *        type could be configured. This can still be extended.
 *
 * \todo  It must be possible to give it a file handle or file name
 *        -> This is now possible. One can give a directory and a file name.
 *
 * \todo  The name should indicate VTU
 *        -> This has not been implemented directly, but the PlotOptions class contains a flag
 *        discriminating which plot_vtu, ... function should be chosen in the .C file. However, at
 *        the moment only vtu output is possible.
 * 
 * One of the advantages of this software package is the strict discrimination between the topology
 * and the geometry of the domain @f$\Omega@f$. Thus, one can exemplarily define a single topology
 * (the one of a cube) to approximate PDEs that live on the cube's boundary and PDEs that live on a
 * sphere, since their topology is the same. However, different geometries have to be defined, since
 * these obviously are not equal. Thus, all parts of the code that involve communication and/or
 * solving systems of equations are reusable in a much more general (than the standard) sense.
 * Beyond that, absurd (on first sight) domains can be defined easily. This also covers variously
 * periodic domains, for example.
 *
 * @tparam  hyedge_dim   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * @tparam  space_dim       The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyperedge_dim.
 * @tparam  poly_degree     Maximal polynomial degree of utilized test and trial functions.
 *
 * @authors   Guido Kanschat, University of Heidelberg, 2019.
 * @authors   Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************/
template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
class PlotOptions
{
  private:
    const HDGHyperGraph < compute_n_dofs_per_node(hyperedge_dim, poly_degree),
                          Topology::HyperGraph_Cubic< hyperedge_dim, space_dim >,
                          Geometry::HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >&
                        hyper_graph_;
    const DiffusionSolver_RegularQuad<hyperedge_dim, poly_degree, 2 * poly_degree>&
                        local_solver_;
  public:
    std::string outputDir, fileName, fileEnding;
    unsigned int fileNumber;
    bool printFileNumber, incrementFileNumber;
    
    PlotOptions(HDGHyperGraph < compute_n_dofs_per_node(hyperedge_dim, poly_degree),
                                Topology::HyperGraph_Cubic< hyperedge_dim, space_dim >,
                                Geometry::HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >&
                              hyper_graph,
                DiffusionSolver_RegularQuad<hyperedge_dim, poly_degree, 2*poly_degree>&
                  local_solver);
    
    const HDGHyperGraph < compute_n_dofs_per_node(hyperedge_dim, poly_degree),
                          Topology::HyperGraph_Cubic< hyperedge_dim, space_dim >,
                          Geometry::HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >&
                        hyper_graph();
    const DiffusionSolver_RegularQuad<hyperedge_dim, poly_degree, 2 * poly_degree>&
            local_solver();
}; // end class PlotOptions


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
void plot
(std::vector<double> lambda, PlotOptions<hyperedge_dim,space_dim,poly_degree>& plotOpt);

#endif
