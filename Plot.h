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
 *        LocalSolver are template classes! -> Done, but are you ok with this?
 * 
 * This class contains all information needed for creating a plot from a @c std::vector, i.e., it
 * contains an @c HDGHyperGraph and a @c LocalSolver and several fields to encode a file name, the
 * name of the computed unknowns, etc. It needs to be passed to the main @c plot function to
 * configure its behavior.
 *
 * @tparam  HyperGraphT     Template parameter describing the precise class of the @c HDGHyperGraph,
 *                          i.e., it contains an @c HDGHyperGraph with chosen template parameters
 *                          describing its topology, geometry, etc.
 * @tparam  LocalSolverT    Template parameter describing the precise class of the local solver,
 *                          i.e., it contains an local solver for a specific equation living on the
 *                          hypergraph.
 *
 * @authors   Guido Kanschat, University of Heidelberg, 2020.
 * @authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template <class HyperGraphT, class LocalSolverT>
class PlotOptions
{
  private:
    /*!*********************************************************************************************
     * @brief   Reference to the hypergraph to be plotted.
     *
     * A @c HyperGraphT& containing a @c const reference to a @c HDGHyperGraph.
     **********************************************************************************************/
    const HyperGraphT& hyper_graph_;
    /*!*********************************************************************************************
     * @brief   Reference to the local olver used to solve the equation to be plotted.
     *
     * A @c LocalSolverT& containing a @c const reference to a local solver class that represents
     * the equation to be plotted and whose domain is the hypergraph.
     **********************************************************************************************/
    const LocalSolverT& local_solver_;
  public:
    /*!*********************************************************************************************
     * @brief   Name of the directory to put the output into.
     *
     * This @c std::string describes the directory the output is created in. Default is "output".
     **********************************************************************************************/
    std::string outputDir;
    /*!*********************************************************************************************
     * @brief   Name of the file emcoding the plot.
     *
     * This @c std::string describes the name of the file for the output plot. Default is "example".
     **********************************************************************************************/
    std::string fileName;
    /*!*********************************************************************************************
     * @brief   File ending and also way of plotting.
     *
     * This @c std::string describes the name of the file ending for the output plot. Thus, it also
     * characterizes which applications can read the output files properly. Default and currently
     * only option is "vtu" for Paraview.
     **********************************************************************************************/
    std::string fileEnding;
    /*!*********************************************************************************************
     * @brief   Number of the plot file.
     *
     * This @c unsigned @c int describes the number of the plot file which is created. If a problem
     * is solved repeatedly (e.g. parabolic problem, local refinement, ...) this number indicates
     * the number of the file (e.g. time step, refinement step, ...). Default is 0.
     **********************************************************************************************/
    unsigned int fileNumber;
    /*!*********************************************************************************************
     * @brief   Decicde whether @c fileNumber is part of the file name.
     *
     * This @c boolean discriminates whether the @c fileNumber should appear within the name of the
     * file (true) or not (false). Default is true.
     **********************************************************************************************/
    bool printFileNumber;
    /*!*********************************************************************************************
     * @brief   Decicde whether @c fileNumber is incremented after plotting.
     *
     * This @c boolean discriminates whether the @c fileNumber should be incremented after a file
     * has been written (true) or not (false). Default is true.
     **********************************************************************************************/
    bool incrementFileNumber;
    /*!*********************************************************************************************
     * @brief   Construct a @c PlotOptions class object containing default values.
     *
     * Constructs a @c PlotOptions object containing default options and information about the
     * hypergraph and the local solver.
     * 
     * @param   hyper_graph     Reference to @c HyperGraphT& representing the hypergraph.
     * @param   local_solver    Reference to @c LocalSolverT& representing the local solver.
     **********************************************************************************************/
    PlotOptions(HyperGraphT& hyper_graph, LocalSolverT& local_solver);
    /*!*********************************************************************************************
     * @brief   Return reference to hypergraph.
     *
     * Return a reference to the @c HDGHyperGraph where the solution should be plotted on.
     *
     * @retval  hyper_graph     Reference to @c HyperGraphT& representing the hypergraph.
     **********************************************************************************************/
    const HyperGraphT& hyper_graph();
    /*!*********************************************************************************************
     * @brief   Return reference to local solver.
     *
     * Return a reference to the local solver representing the equation solved on the hypergraph.
     *
     * @retval  local_solver    Reference to @c LocalSolverT& representing the local solver.
     **********************************************************************************************/
    const LocalSolverT& local_solver();
}; // end of class PlotOptions

/*!*************************************************************************************************
 * @brief   Function plotting the solution of an equation on a hypergraph.
 *
 * Creates a file according to set plotting options in @c plotOpt. This file contains the solution
 * of the PDE defined in @c plotOpt having the representation @c lambda in terms of its skeletal
 * degrees of freedom (related to skeletal variable lambda).
 *
 * @tparam  HyperGraphT     Template parameter describing the precise class of the @c HDGHyperGraph,
 *                          i.e., it contains an @c HDGHyperGraph with chosen template parameters
 *                          describing its topology, geometry, etc.
 * @tparam  LocalSolverT    Template parameter describing the precise class of the local solver,
 *                          i.e., it contains an local solver for a specific equation living on the
 *                          hypergraph.
 * @param   lambda          @c std::vector containing the skeletal degrees of freedom encoding the
 *                          representation of the unique solution.
 * @param   plotOpt         @c PlotOptions object containing plotting options.
 *
 * @authors   Guido Kanschat, University of Heidelberg, 2020.
 * @authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template <class HyperGraphT, class LocalSolverT>
void plot
(std::vector<double> lambda, PlotOptions<HyperGraphT, LocalSolverT>& plotOpt);

#endif // end of ifndef PLOT_H
