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
 * @todo  Turn PlotOptions class into fully templated class where the HDGHyperGraph and the
 *        LocalSolver are template classes!
 * 
 * This class contains all information needed for creating a plot from a std::vector, i.e., it
 * contains an \c HDGHyperGraph and a \c LocalSolver and several fields to encode a file name, the
 * name of the computed unknowns, etc. It needs to be passed to the main \c plot function to
 * configure its behavior.
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
template <class HyperGraphT, class LocalSolverT>
class PlotOptions
{
  private:
    const HyperGraphT& hyper_graph_;
    const LocalSolverT& local_solver_;
  public:
    std::string outputDir, fileName, fileEnding;
    unsigned int fileNumber;
    bool printFileNumber, incrementFileNumber;
    
    PlotOptions(HyperGraphT& hyper_graph, LocalSolverT& local_solver);
    
    const HyperGraphT& hyper_graph();
    const LocalSolverT& local_solver();
}; // end class PlotOptions


template <class HyperGraphT, class LocalSolverT>
void plot
(std::vector<double> lambda, PlotOptions<HyperGraphT, LocalSolverT>& plotOpt);

#endif
