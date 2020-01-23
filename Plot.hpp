/*!*************************************************************************************************
 * @file    Plot.hpp
 * @brief   This file provides the function plot and the class PlotOptions.
 *
 * This file provides a file @c plot which takes a @c PlotOptions instatiation to plot given data
 * encoded in a @c std::vector. Additionally, it defines the @c PlotOptions class which is closely
 * rekated to the function @c plot.
 * 
 * This file is an .hpp file, since class and function are compiled "dynamically" depending on the
 * considered problem in Python or C++ code. Dynamically means, that either, when the C++ problem
 * or Python's ClassWrapper are compiled, the relevant template parameters for the respective class
 * and functions of this file are deduced and the needed versions are compiled.
 *
 * @authors   Guido Kanschat, University of Heidelberg, 2020.
 * @authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#ifndef PLOT_HPP
#define PLOT_HPP

// Includes needed for external communication.
// These also would ben included when splitted in .C and .h files.
#include "HDGHyperGraph.h"
#include "LocalSolvers.h"

// Includes solely needed for implementation of the different functions.
// These would not be included when splitted in .C and .h files.
#include "HyAssert.h"
#include "TypeDefs.h"
#include <fstream>
#include <iomanip>
#include <cmath>

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
    PlotOptions(HyperGraphT& hyper_graph, LocalSolverT& local_solver)
    : hyper_graph_(hyper_graph), local_solver_(local_solver),
      outputDir("output"), fileName("example"), fileEnding("vtu"), fileNumber(0),
      printFileNumber(true), incrementFileNumber(true)  { };
    /*!*********************************************************************************************
     * @brief   Return reference to hypergraph.
     *
     * Return a reference to the @c HDGHyperGraph where the solution should be plotted on.
     *
     * @retval  hyper_graph     Reference to @c HyperGraphT& representing the hypergraph.
     **********************************************************************************************/
    const HyperGraphT& hyper_graph()
    {
      return hyper_graph_;
    };
    /*!*********************************************************************************************
     * @brief   Return reference to local solver.
     *
     * Return a reference to the local solver representing the equation solved on the hypergraph.
     *
     * @retval  local_solver    Reference to @c LocalSolverT& representing the local solver.
     **********************************************************************************************/
    const LocalSolverT& local_solver()
    {
      return local_solver_;
    };
}; // end of class PlotOptions

/*!*************************************************************************************************
 * @brief   Function plotting the solution of an equation on a hypergraph in vtu format.
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
void plot_vtu(std::vector<double> lambda, PlotOptions<HyperGraphT,LocalSolverT>& plotOpt)
{
  constexpr unsigned int hyperedge_dim = HyperGraphT::hyperedge_dimension();
  constexpr unsigned int space_dim = HyperGraphT::space_dimension();
  constexpr unsigned int points_per_hyperedge = pow(2, hyperedge_dim);
  
  hyperedge_index_type num_of_hyperedges = plotOpt.hyper_graph().num_of_hyperedges();
  point_index_type num_of_points = points_per_hyperedge * num_of_hyperedges;
  
  static_assert (hyperedge_dim <= 3);
  unsigned int element_id;
  if constexpr (hyperedge_dim == 1)       element_id = 3;
  else if constexpr (hyperedge_dim == 2)  element_id = 8;
  else if constexpr (hyperedge_dim == 3)  element_id = 11;
  
  Point<space_dim> point;
  std::ofstream myfile;
  
  std::string filename = plotOpt.outputDir;
  filename.append("/"); filename.append(plotOpt.fileName);
  if(plotOpt.printFileNumber)
  {
    filename.append("."); filename.append(std::to_string(plotOpt.fileNumber));
  }
  filename.append(".vtu");
  if(plotOpt.incrementFileNumber)  ++(plotOpt.fileNumber);
  
	myfile.open(filename);
	myfile << "<?xml version=\"1.0\"?>"  << std::endl;
	myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
         << "compressor=\"vtkZLibDataCompressor\">" << std::endl;
	myfile << "  <UnstructuredGrid>" << std::endl;
	myfile << "    <Piece NumberOfPoints=\"" << num_of_points << "\" NumberOfCells= \""
         << num_of_hyperedges << "\">" << std::endl;
	myfile << "      <Points>" << std::endl;
	myfile << "        <Datastd::array type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"
         << std::endl;
  
  for (hyperedge_index_type he_number = 0; he_number < num_of_hyperedges; ++he_number)
  {
    auto hyperedge_geometry = plotOpt.hyper_graph().hyperedge_geometry(he_number);
    for (unsigned int pt_number = 0; pt_number < points_per_hyperedge; ++pt_number)
    {
      myfile << "        ";
      point = hyperedge_geometry.point(pt_number);
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        myfile << "  " << std::fixed << std::scientific << std::setprecision(3) << point[dim];
      for (unsigned int dim = space_dim; dim < 3; ++dim)
        myfile << "  0.0";
      myfile << std::endl;
    }
  }
  
  myfile << "        </Datastd::array>" << std::endl;
  myfile << "      </Points>" << std::endl;
	myfile << "      <Cells>" << std::endl;
	myfile << "        <Datastd::array type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
         << std::endl;
  myfile << "        ";
  
	for (point_index_type pt_number = 0; pt_number < num_of_points; ++pt_number)
    myfile << "  " << pt_number;
  myfile << std::endl;
  
  myfile << "        </Datastd::array>" << std::endl;
  myfile << "        <Datastd::array type=\"Int32\" Name=\"offsets\" format=\"ascii\">"
         << std::endl;
  myfile << "        ";
  
  for (point_index_type pt_number = points_per_hyperedge; pt_number <= num_of_points;
       pt_number += points_per_hyperedge)
    myfile << "  " << pt_number;
  
  myfile << std::endl;
  
  myfile << "        </Datastd::array>" << std::endl;
	myfile << "        <Datastd::array type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
  myfile << "        ";
  
  for (hyperedge_index_type he_number = 0; he_number < num_of_hyperedges; ++he_number)
    myfile << "  " << element_id;
  myfile << std::endl;
  
  myfile << "        </Datastd::array>" << std::endl;
	myfile << "      </Cells>" << std::endl;
  
  
  myfile << "      <PointData Scalars=\"" << "example_scalar" << "\" std::vectors=\""
         << "example_std::vector" << "\">" << std::endl;
  myfile << "        <Datastd::array type=\"Float32\" Name=\"" << "dual" 
         << "\" NumberOfComponents=\"" << hyperedge_dim << "\" format=\"ascii\">" << std::endl;
    
  std::array< std::array<double, HyperGraphT::n_dof_per_node() > , 2*hyperedge_dim > hyperedge_dofs;
  std::array<unsigned int, 2*hyperedge_dim> hyperedge_hypernodes;
  std::array<double, compute_n_corners_of_cube(hyperedge_dim)> local_primal;
  std::array< std::array<double, hyperedge_dim> , compute_n_corners_of_cube(hyperedge_dim) >
    local_dual;
  
  for (hyperedge_index_type he_number = 0; he_number < num_of_hyperedges; ++he_number)
  {
    hyperedge_hypernodes = 
      plotOpt.hyper_graph().hyperedge_topology(he_number).get_hypernode_indices();
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyperedge_dofs[hypernode] = plotOpt.hyper_graph().hypernode_factory().
        get_dof_values(hyperedge_hypernodes[hypernode], lambda);
    local_dual = plotOpt.local_solver().dual_in_corners_from_lambda(hyperedge_dofs);
    myfile << "      ";
    for (unsigned int corner = 0; corner < compute_n_corners_of_cube(hyperedge_dim); ++corner)
    {
      myfile << "  ";
      for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
        myfile << "  " << local_dual[corner][dim];
    }
    myfile << std::endl;
  }

  myfile << "        </Datastd::array>" << std::endl;
  myfile << "        <Datastd::array type=\"Float32\" Name=\"" << "primal" 
         << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  
  
  for (hyperedge_index_type he_number = 0; he_number < num_of_hyperedges; ++he_number)
  {
    hyperedge_hypernodes = plotOpt.hyper_graph().hyperedge_topology(he_number).
      get_hypernode_indices();
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyperedge_dofs[hypernode] = plotOpt.hyper_graph().hypernode_factory().
        get_dof_values(hyperedge_hypernodes[hypernode], lambda);
    local_primal = plotOpt.local_solver().primal_in_corners_from_lambda(hyperedge_dofs);
    myfile << "        ";
    for (unsigned int corner = 0; corner < compute_n_corners_of_cube(hyperedge_dim); ++corner)
      myfile << "  " << local_primal[corner];
    myfile << std::endl;
  }

  myfile << "        </Datastd::array>" << std::endl;
  myfile << "      </PointData>" << std::endl;
  

	myfile << "    </Piece>" << std::endl;
	myfile << "  </UnstructuredGrid>" << std::endl;
	myfile << "</VTKFile>" << std::endl;
  myfile.close();
}; // end of void plot_vtu

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
(std::vector<double> lambda, PlotOptions<HyperGraphT, LocalSolverT>& plotOpt)
{
  hy_assert( plotOpt.fileEnding == "vtu" , 
             "Only file ending vtu is supported at the moment. Your choice has been "
             << plotOpt.fileEnding << ", which is invalid.");
  hy_assert( !plotOpt.fileName.empty() , "File name must not be empty!" );
  hy_assert( !plotOpt.outputDir.empty() , "Ouput directory must not be empty!" );
  plot_vtu<HyperGraphT,LocalSolverT>(lambda, plotOpt);
}; // end of void plot

#endif // end of ifndef PLOT_HPP
