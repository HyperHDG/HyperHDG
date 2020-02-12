/*!*************************************************************************************************
 * \file    Plot.hxx
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

#ifndef PLOT_HXX
#define PLOT_HXX

#include <HDGHyperGraph.hxx>
#include <Hypercube.hxx>

// Includes solely needed for implementation of the different functions.
// These would not be included when splitted in .C and .h files.
#include <HyAssert.hxx>
#include <TypeDefs.hxx>
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
    /*!*********************************************************************************************
     * \brief Include the nodes with their function values into the plot
     *
     * Defaults to `false`.
     **********************************************************************************************/
    bool plot_nodes;
    /*!*********************************************************************************************
     * \brief Include the hyperedges with their function values into the plot
     *
     * Defaults to `true`.
     **********************************************************************************************/
    bool plot_edges;
    /*!*********************************************************************************************
     * \brief Number of subintervals for plotting
     *
     * When plotting an interval, it is split into n_subintervals intervals such that higher order
     * polynomials and other functions can be displayed appropriately. When plotting higher
     * dimensional objects, this subdivision is applied accordingly.
     *
     * Defaults to 1.
     **********************************************************************************************/
    unsigned int n_subintervals;  
    /*!*********************************************************************************************
     * \brief   Construct a \c PlotOptions class object containing default values.
     *
     * Constructs a \c PlotOptions object containing default options and information about the
     * hypergraph and the local solver.
     **********************************************************************************************/
    PlotOptions();
}; // end of class PlotOptions

/*!*************************************************************************************************
 * \brief   Function plotting the solution of an equation on a hypergraph in vtu format.
 *
 * Creates a file according to set plotting options in \c plot_options. This file contains the solution
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
 * \param   plot_options    PlotOptions object containing plotting options.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template <class HyperGraphT, class LocalSolverT>
void plot(const HyperGraphT& hyper_graph,
	  const LocalSolverT& local_solver,
	  const std::vector<double>& lambda,
	  const PlotOptions& plot_options);


PlotOptions::PlotOptions()
  : outputDir("output"), fileName("example"), fileEnding("vtu"), fileNumber(0),
    printFileNumber(true), incrementFileNumber(true), plot_nodes(false), plot_edges(true),
    n_subintervals(1)
{}

/*!**********************************************************************
 * \brief Auxiliary functions for writing graphics files
 ***********************************************************************/
namespace PlotFunctions
{
  /*!**********************************************************************
   * \brief Auxiliary function for writing geometry section VTU files
   *
   * This function plots the geometry part of an unstructured mesh in
   * a VTU file. The typical file structure is
   *
   * ```
   * <Preamble/>
   * <UnstructuredGrid>
   * <Geometry/>
   * <Data/>
   * </UnstructuredGrid>
   *
   * This function writes the `<Geometry>` part of the structure.
   ************************************************************************/
  template
  <class HyperGraphT, typename hyEdge_index_t = unsigned int, typename pt_index_t = unsigned int >
  void plot_vtu_unstructured_geometry(std::ostream& output,
				      const HyperGraphT& hyper_graph,
				      const PlotOptions& plot_options)
  {
    constexpr unsigned int hyEdge_dim = HyperGraphT::hyEdge_dimension();
    constexpr unsigned int space_dim = HyperGraphT::space_dimension();
    
    const hyEdge_index_t n_hyEdges = hyper_graph.n_hyEdges();
    const unsigned int points_per_hyEdge = 1 << hyEdge_dim;
    
    pt_index_t n_points = points_per_hyEdge * n_hyEdges;
    
    static_assert (hyEdge_dim <= 3);
    unsigned int element_id;
    if constexpr (hyEdge_dim == 1)       element_id = 3;
    else if constexpr (hyEdge_dim == 2)  element_id = 8;
    else if constexpr (hyEdge_dim == 3)  element_id = 11;
    
    output << "    <Piece NumberOfPoints=\"" << n_points << "\" NumberOfCells= \"" << n_hyEdges << "\">" << std::endl;
    output << "      <Points>" << std::endl;
    output << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;
    
    for (hyEdge_index_t he_number = 0; he_number < n_hyEdges; ++he_number)
      {
	auto hyEdge_geometry = hyper_graph.hyEdge_geometry(he_number);
	for (unsigned int pt_number = 0; pt_number < points_per_hyEdge; ++pt_number)
	  {
	    output << "        ";
	    const Point<space_dim> point = hyEdge_geometry.point(pt_number);
	    for (unsigned int dim = 0; dim < space_dim; ++dim)
	      output << "  " << std::fixed << std::scientific << std::setprecision(3) << point[dim];
	    for (unsigned int dim = space_dim; dim < 3; ++dim)
	      output << "  0.0";
	    output << std::endl;
	  }
      }
    
    output << "        </DataArray>" << std::endl;
    output << "      </Points>" << std::endl;
    output << "      <Cells>" << std::endl;
    output << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << std::endl;
    output << "        ";
    
    for (pt_index_t pt_number = 0; pt_number < n_points; ++pt_number)
      output << "  " << pt_number;
    output << std::endl;
    
    output << "        </DataArray>" << std::endl;
    output << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    output << "        ";
    
    for (pt_index_t pt_number = points_per_hyEdge; pt_number <= n_points;
	 pt_number += points_per_hyEdge)
      output << "  " << pt_number;
    output << std::endl;
    
    output << "        </DataArray>" << std::endl;
    output << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    output << "        ";
    
    for (hyEdge_index_t he_number = 0; he_number < n_hyEdges; ++he_number)
      output << "  " << element_id;
    output << std::endl;
    
    output << "        </DataArray>" << std::endl;
    output << "      </Cells>" << std::endl;
  }
}

template <class HyperGraphT, class LocalSolverT, typename hyEdge_index_t = unsigned int>
void plot_vtu(const HyperGraphT& hyper_graph,
	      const LocalSolverT& local_solver,
	      const std::vector<double>& lambda,
	      const PlotOptions& plot_options)
{
  constexpr unsigned int hyEdge_dim = HyperGraphT::hyEdge_dimension();
  
  const hyEdge_index_t n_hyEdges = hyper_graph.n_hyEdges();
  //  const unsigned int points_per_hyEdge = 1 << hyEdge_dim;
  const std::array<double, 2> abscissas = {0., 1.};
  
  static_assert (hyEdge_dim <= 3);
  
  std::ofstream myfile;
  
  std::string filename = plot_options.outputDir;
  filename.append("/"); filename.append(plot_options.fileName);
  if(plot_options.printFileNumber)
  {
    filename.append("."); filename.append(std::to_string(plot_options.fileNumber));
  }
  filename.append(".vtu");
  
  myfile.open(filename);
  myfile << "<?xml version=\"1.0\"?>"  << std::endl;
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl;
  myfile << "  <UnstructuredGrid>" << std::endl;

  PlotFunctions::plot_vtu_unstructured_geometry(myfile, hyper_graph, plot_options);
  
  myfile << "      <PointData Scalars=\"" << "example_scalar" << "\" Vectors=\"" << "example_vector" << "\">" << std::endl;
  myfile << "        <DataArray type=\"Float32\" Name=\"" << "dual" << "\" NumberOfComponents=\"" << hyEdge_dim << "\" format=\"ascii\">" << std::endl;
    
  std::array< std::array<double, HyperGraphT::n_dofs_per_node() > , 2*hyEdge_dim > hyEdge_dofs;
  std::array<unsigned int, 2*hyEdge_dim> hyEdge_hyNodes;
  std::array<double, Hypercube<hyEdge_dim>::n_vertices() > local_primal; // CHANGED FOR ABSCISSAS!
  std::array< std::array<double, hyEdge_dim> , Hypercube<hyEdge_dim>::n_vertices() > local_dual; // CHANGED FOR ABSCISSAS!
  
  for (hyEdge_index_t he_number = 0; he_number < n_hyEdges; ++he_number)
  {
    hyEdge_hyNodes = hyper_graph.hyEdge_topology(he_number).get_hyNode_indices();
    for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
      hyEdge_dofs[hyNode] = hyper_graph.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], lambda);
    if constexpr ( LocalSolverT::use_geometry() )
      local_dual = local_solver.dual_at_dyadic(abscissas, hyper_graph.hyEdge_geometry(he_number), hyEdge_dofs);
    else
      local_dual = local_solver.dual_at_dyadic(abscissas, hyEdge_dofs);
    myfile << "      ";
    for (unsigned int corner = 0; corner < Hypercube<hyEdge_dim>::n_vertices(); ++corner)
    {
      myfile << "  ";
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
        myfile << "  " << local_dual[corner][dim];
    }
    myfile << std::endl;
  }

  myfile << "        </DataArray>" << std::endl;
  myfile << "        <DataArray type=\"Float32\" Name=\"" << "primal" << "\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
  
  for (hyEdge_index_t he_number = 0; he_number < n_hyEdges; ++he_number)
  {
    hyEdge_hyNodes = hyper_graph.hyEdge_topology(he_number).get_hyNode_indices();
    for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
      hyEdge_dofs[hyNode] = hyper_graph.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], lambda);
    if constexpr ( LocalSolverT::use_geometry() )
      local_primal = local_solver.primal_at_dyadic(abscissas, hyper_graph.hyEdge_geometry(he_number), hyEdge_dofs);
    else
      local_primal = local_solver.primal_at_dyadic(abscissas, hyEdge_dofs);
    myfile << "        ";
    for (unsigned int corner = 0; corner < Hypercube<hyEdge_dim>::n_vertices(); ++corner)
      myfile << "  " << local_primal[corner];
    myfile << std::endl;
  }

  myfile << "        </DataArray>" << std::endl;
  myfile << "      </PointData>" << std::endl;  
  
  myfile << "    </Piece>" << std::endl;
  myfile << "  </UnstructuredGrid>" << std::endl;
  myfile << "</VTKFile>" << std::endl;
  myfile.close();
} // end of void plot_vtu

/*!*************************************************************************************************
 * \brief   Function plotting the solution of an equation on a hypergraph.
 * 
 * \todo    Plot function has been deactivated for elasticity!
 *
 * Creates a file according to set plotting options in \c plot_options. This file contains the solution
 * of the PDE defined in \c plot_options having the representation \c lambda in terms of its skeletal
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
 * \param   plot_options         \c PlotOptions object containing plotting options.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template <class HyperGraphT, class LocalSolverT>
void plot(const HyperGraphT& hyper_graph,
	  const LocalSolverT& local_solver,
	  const std::vector<double>& lambda,
	  const PlotOptions& plot_options)
{
  hy_assert( plot_options.fileEnding == "vtu" , 
             "Only file ending vtu is supported at the moment. Your choice has been "
             << plot_options.fileEnding << ", which is invalid.");
  hy_assert( !plot_options.fileName.empty() , "File name must not be empty!" );
  hy_assert( !plot_options.outputDir.empty() , "Ouput directory must not be empty!" );
  plot_vtu(hyper_graph, local_solver, lambda, plot_options);
} // end of void plot

template<unsigned int hyEdge_dim, unsigned int space_dim, unsigned int max_poly_degree, unsigned int max_quad_degree>
class ElasticitySolver_RegularQuad;

template<class HyperGraphT, unsigned int hdim, unsigned int sdim, unsigned int pd, unsigned int qd>
void plot(const HyperGraphT& hyper_graph,
	  const ElasticitySolver_RegularQuad<hdim,sdim,pd,qd>& local_solver,
	  const std::vector<double>& lambda,
	  const PlotOptions& plot_options)
{
}


#endif // end of ifndef PLOT_HXX
