/*!*************************************************************************************************
 * \file    Plot.hpp
 * \brief   This file provides the function plot and the class PlotOptions.
 *
 * \todo    Plot has been deactivated to implement elasticity!
 * 
 * This file provides a file \c plot which takes a \c PlotOptions instatiation to plot given data
 * encoded in a \c std::vector. Additionally, it defines the \c PlotOptions class which is closely
 * rekated to the function \c plot.
 * 
 * This file is an .hpp file, since class and function are compiled "dynamically" depending on the
 * considered problem in Python or C++ code. Dynamically means, that either, when the C++ problem
 * or Python's ClassWrapper are compiled, the relevant template parameters for the respective class
 * and functions of this file are deduced and the needed versions are compiled.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#ifndef PLOT_HPP
#define PLOT_HPP

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
 * \brief   A class storing options for plotting.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
class PlotOptions
{
  public:
    /*!*********************************************************************************************
     * \brief   Name of the directory to put the output into.
     *
     * This \c std::string describes the directory the output is created in. Default is "output".
     **********************************************************************************************/
    std::string outputDir;
    /*!*********************************************************************************************
     * \brief   Name of the file plotted.
     *
     * This \c std::string describes the name of the file for the output plot. Default is "example".
     **********************************************************************************************/
    std::string fileName;
    /*!*********************************************************************************************
     * \brief   File ending and also way of plotting.
     *
     * This \c std::string describes the name of the file ending for the output plot. Thus, it also
     * characterizes which applications can read the output files properly. Default and currently
     * only option is "vtu" for Paraview.
     * \todo G recommends to make this an enum. Text matching is always apita
     **********************************************************************************************/
    std::string fileEnding;
    /*!*********************************************************************************************
     * \brief   Number of the plot file.
     *
     * This \c unsigned \c int describes the number of the plot file which is created. If a problem
     * is solved repeatedly (e.g. parabolic problem, local refinement, ...) this number indicates
     * the number of the file (e.g. time step, refinement step, ...). Default is 0.
     **********************************************************************************************/
    unsigned int fileNumber;
    /*!*********************************************************************************************
     * \brief   Decide whether \c fileNumber is part of the file name.
     *
     * This \c boolean discriminates whether the \c fileNumber should appear within the name of the
     * file (true) or not (false). Default is true.
     **********************************************************************************************/
    bool printFileNumber;
    /*!*********************************************************************************************
     * \brief   Decide whether \c fileNumber is incremented after plotting.
     *
     * This \c boolean discriminates whether the \c fileNumber should be incremented after a file
     * has been written (true) or not (false). Default is true.
     *
     * \todo This could not be implemented in your version with call
     * by value, nor can it be implemented in mine with a const
     * reference.
     **********************************************************************************************/
    bool incrementFileNumber;
    /*!
     * \brief Include the nodes with their function values into the plot
     *
     * Defaults to `false`.
     */
    bool plot_nodes;
    /*!
     * \brief Include the hyperedges with their function values into the plot
     *
     * Defaults to `true`.
     */
    bool plot_edges;
    /*!
     * \brief Number of subintervals for plotting
     *
     * When plotting an interval, it is split into n_subintervals
     * intervals such that higher order polynomials and other
     * functions can be displayed appropriately. When plotting higher
     * dimensional objects, this subdivision is applied accordingly.
     *
     * Defaults to 1.
     */
    unsigned int n_subintervals;  
    /*!*********************************************************************************************
     * \brief   Construct a \c PlotOptions class object containing default values.
     *
     * Constructs a \c PlotOptions object containing default options and information about the
     * hypergraph and the local solver.
     * 
     * \param   hyper_graph     Reference to \c HyperGraphT& representing the hypergraph.
     * \param   local_solver    Reference to \c LocalSolverT& representing the local solver.
     **********************************************************************************************/
    PlotOptions();
}; // end of class PlotOptions

/*!*************************************************************************************************
 * \brief   Function plotting the solution of an equation on a hypergraph in vtu format.
 *
 * Creates a file according to set plotting options in \c plotOpt. This file contains the solution
 * of the PDE defined in \c plotOpt having the representation \c lambda in terms of its skeletal
 * degrees of freedom (related to skeletal variable lambda).
 *
 * \tparam  HyperGraphT     Template parameter describing the precise class of the \c HDGHyperGraph,
 *                          i.e., it contains an \c HDGHyperGraph with chosen template parameters
 *                          describing its topology, geometry, etc.
 * \tparam  LocalSolverT    Template parameter describing the precise class of the local solver,
 *                          i.e., it contains an local solver for a specific equation living on the
 *                          hypergraph.
 * \param   lambda          \c std::vector containing the skeletal degrees of freedom encoding the
 *                          representation of the unique solution.
 * \param   plotOpt         \c PlotOptions object containing plotting options.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template <class HyperGraphT, class LocalSolverT>
void plot(const HyperGraphT& hyper_graph,
	  const LocalSolverT& local_solver,
	  const std::vector<double>& lambda,
	  const PlotOptions& plotOpt);


PlotOptions::PlotOptions()
  : outputDir("output"), fileName("example"), fileEnding("vtu"), fileNumber(0),
    printFileNumber(true), incrementFileNumber(true), plot_nodes(false), plot_edges(true),
    n_subintervals(1)
{}


template <class HyperGraphT, class LocalSolverT>
void plot_vtu(const HyperGraphT& hyper_graph,
	      const LocalSolverT& local_solver,
	      const std::vector<double>& lambda,
	      const PlotOptions& plotOpt)
{
  constexpr unsigned int hyperedge_dim = HyperGraphT::hyperedge_dimension();
  constexpr unsigned int space_dim = HyperGraphT::space_dimension();
  
  const hyperedge_index_type n_hyperedges = hyper_graph.n_hyperedges();
  const unsigned int points_per_hyperedge = 1 << hyperedge_dim;
  
  point_index_type n_points = points_per_hyperedge * n_hyperedges;
  
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
  
  myfile.open(filename);
  myfile << "<?xml version=\"1.0\"?>"  << std::endl;
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
  myfile << "  <UnstructuredGrid>" << std::endl;
  myfile << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells= \"" << n_hyperedges << "\">" << std::endl;
  myfile << "      <Points>" << std::endl;
  myfile << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
  
  for (hyperedge_index_type he_number = 0; he_number < n_hyperedges; ++he_number)
  {
    auto hyperedge_geometry = hyper_graph.hyperedge_geometry(he_number);
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
  
  myfile << "        </DataArray>" << std::endl;
  myfile << "      </Points>" << std::endl;
	myfile << "      <Cells>" << std::endl;
	myfile << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
  myfile << "        ";
  
	for (point_index_type pt_number = 0; pt_number < n_points; ++pt_number)
    myfile << "  " << pt_number;
  myfile << std::endl;
  
  myfile << "        </DataArray>" << std::endl;
  myfile << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
  myfile << "        ";
  
  for (point_index_type pt_number = points_per_hyperedge; pt_number <= n_points;
       pt_number += points_per_hyperedge)
    myfile << "  " << pt_number;
  myfile << std::endl;
  
  myfile << "        </DataArray>" << std::endl;
	myfile << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
  myfile << "        ";
  
  for (hyperedge_index_type he_number = 0; he_number < n_hyperedges; ++he_number)
    myfile << "  " << element_id;
  myfile << std::endl;
  
  myfile << "        </DataArray>" << std::endl;
	myfile << "      </Cells>" << std::endl;
  
  
  myfile << "      <PointData Scalars=\"" << "example_scalar" << "\" Vectors=\"" << "example_vector" << "\">" << std::endl;
  myfile << "        <DataArray type=\"Float32\" Name=\"" << "dual" << "\" NumberOfComponents=\"" << hyperedge_dim << "\" format=\"ascii\">" << std::endl;
    
  std::array< std::array<double, HyperGraphT::n_dof_per_node() > , 2*hyperedge_dim > hyperedge_dofs;
  std::array<unsigned int, 2*hyperedge_dim> hyperedge_hypernodes;
  std::array<double, compute_n_corners_of_cube(hyperedge_dim)> local_primal;
  std::array< std::array<double, hyperedge_dim> , compute_n_corners_of_cube(hyperedge_dim) >
    local_dual;
  
  for (hyperedge_index_type he_number = 0; he_number < n_hyperedges; ++he_number)
  {
    hyperedge_hypernodes = hyper_graph.hyperedge_topology(he_number).get_hypernode_indices();
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyperedge_dofs[hypernode] = hyper_graph.hypernode_factory().get_dof_values(hyperedge_hypernodes[hypernode], lambda);
    local_dual = local_solver.dual_in_corners_from_lambda(hyperedge_dofs);
    myfile << "      ";
    for (unsigned int corner = 0; corner < compute_n_corners_of_cube(hyperedge_dim); ++corner)
    {
      myfile << "  ";
      for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
        myfile << "  " << local_dual[corner][dim];
    }
    myfile << std::endl;
  }

  myfile << "        </DataArray>" << std::endl;
  myfile << "        <DataArray type=\"Float32\" Name=\"" << "primal" << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  
  
  for (hyperedge_index_type he_number = 0; he_number < n_hyperedges; ++he_number)
  {
    hyperedge_hypernodes = hyper_graph.hyperedge_topology(he_number).get_hypernode_indices();
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyperedge_dofs[hypernode] = hyper_graph.hypernode_factory().get_dof_values(hyperedge_hypernodes[hypernode], lambda);
    local_primal = local_solver.primal_in_corners_from_lambda(hyperedge_dofs);
    myfile << "        ";
    for (unsigned int corner = 0; corner < compute_n_corners_of_cube(hyperedge_dim); ++corner)
      myfile << "  " << local_primal[corner];
    myfile << std::endl;
  }

  myfile << "        </DataArray>" << std::endl;
  myfile << "      </PointData>" << std::endl;
  

	myfile << "    </Piece>" << std::endl;
	myfile << "  </UnstructuredGrid>" << std::endl;
	myfile << "</VTKFile>" << std::endl;
  myfile.close();
}; // end of void plot_vtu

/*!*************************************************************************************************
 * \brief   Function plotting the solution of an equation on a hypergraph.
 * 
 * \todo    Plot function has been deactivated for elasticity!
 *
 * Creates a file according to set plotting options in \c plotOpt. This file contains the solution
 * of the PDE defined in \c plotOpt having the representation \c lambda in terms of its skeletal
 * degrees of freedom (related to skeletal variable lambda).
 *
 * \tparam  HyperGraphT     Template parameter describing the precise class of the \c HDGHyperGraph,
 *                          i.e., it contains an \c HDGHyperGraph with chosen template parameters
 *                          describing its topology, geometry, etc.
 * \tparam  LocalSolverT    Template parameter describing the precise class of the local solver,
 *                          i.e., it contains an local solver for a specific equation living on the
 *                          hypergraph.
 * \param   lambda          \c std::vector containing the skeletal degrees of freedom encoding the
 *                          representation of the unique solution.
 * \param   plotOpt         \c PlotOptions object containing plotting options.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template <class HyperGraphT, class LocalSolverT>
void plot(const HyperGraphT& hyper_graph,
	  const LocalSolverT& local_solver,
	  const std::vector<double>& lambda,
	  const PlotOptions& plotOpt)
{
  hy_assert( plotOpt.fileEnding == "vtu" , 
             "Only file ending vtu is supported at the moment. Your choice has been "
             << plotOpt.fileEnding << ", which is invalid.");
  hy_assert( !plotOpt.fileName.empty() , "File name must not be empty!" );
  hy_assert( !plotOpt.outputDir.empty() , "Ouput directory must not be empty!" );
  plot_vtu(hyper_graph, local_solver, lambda, plotOpt);
}; // end of void plot


template<class HyperGraphT, unsigned int hdim, unsigned int sdim, unsigned int pd, unsigned int qd>
void plot(const HyperGraphT& hyper_graph,
	  const ElasticitySolver_RegularQuad<hdim,sdim,pd,qd>& local_solver,
	  const std::vector<double>& lambda,
	  const PlotOptions& plotOpt)
{
}


#endif // end of ifndef PLOT_HPP
