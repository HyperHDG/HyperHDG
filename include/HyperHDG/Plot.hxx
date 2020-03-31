#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/HDGHyperGraph.hxx>
#include <HyperHDG/Hypercube.hxx>
#include <HyperHDG/SmallVec.hxx>

#include <fstream>
#include <iomanip>
// #include <cmath>

/*!*************************************************************************************************
 * \brief   A class storing options for plotting.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
struct PlotOptions
{
  /*!***********************************************************************************************
   * \brief   Name of the directory to put the output into.
   *
   * This \c std::string describes the directory the output is created in. Default is "output".
   ************************************************************************************************/
  std::string outputDir;
  /*!***********************************************************************************************
   * \brief   Name of the file plotted.
   *
   * This \c std::string describes the name of the file for the output plot. Default is "example".
   ************************************************************************************************/
  std::string fileName;
  /*!***********************************************************************************************
   * \brief   File ending and also way of plotting.
   *
   * This \c std::string describes the name of the file ending for the output plot. Thus, it also
   * characterizes which applications can read the output files properly. Default and currently
   * only option is "vtu" for Paraview.
   * \todo    G recommends to make this an enum. Text matching is always apita
   *          -> A agrees! Should he do this or is G in this part of the file?
   ************************************************************************************************/
  std::string fileEnding;
  /*!***********************************************************************************************
   * \brief   Number of the plot file.
   *
   * This \c unsigned \c int describes the number of the plot file which is created. If a problem is
   * solved repeatedly (e.g. parabolic problem, local refinement, ...) this number indicates the
   * number of the file (e.g. time step, refinement step, ...). Default is 0.
   ************************************************************************************************/
  unsigned int fileNumber;
  /*!***********************************************************************************************
   * \brief   Decide whether \c fileNumber is part of the file name.
   *
   * This \c boolean discriminates whether the \c fileNumber should appear within the name of the
   * file (true) or not (false). Default is true.
   ************************************************************************************************/
  bool printFileNumber;
  /*!***********************************************************************************************
   * \brief   Decide whether \c fileNumber is incremented after plotting.
   *
   * This \c boolean discriminates whether the \c fileNumber should be incremented after a file has
   * been written (true) or not (false). Default is true.
   *
   * \todo    This could not be implemented in your version with call by value, nor can it be
              implemented in mine with a const reference.
   ************************************************************************************************/
  bool incrementFileNumber;
  /*!***********************************************************************************************
   * \brief   Include the nodes with their function values into the plot.
   *
   * Defaults to `false`.
   ************************************************************************************************/
  bool plot_nodes;
  /*!***********************************************************************************************
   * \brief   Include the hyperedges with their function values into the plot.
   *
   * Defaults to `true`.
   ************************************************************************************************/
  bool plot_edges;
  /*!***********************************************************************************************
   * \brief   Number of subintervals for plotting.
   *
   * When plotting an interval, it is split into n_subintervals intervals such that higher order
   * polynomials and other functions can be displayed appropriately. When plotting higher
   * dimensional objects, this subdivision is applied accordingly.
   *
   * Defaults to 1.
   ************************************************************************************************/
  unsigned int n_subintervals;  
  /*!***********************************************************************************************
   * \brief   A factor for scaling each object of the plot locally.
   *
   * This factor defaults to 1 in order to produce a plot of a contiguous domain. If it is chosen
   * less than 1, each edge or node is scaled by this factor around its center.
   ************************************************************************************************/
  float scale;
  /*!***********************************************************************************************
   * \brief   Output a cell data array with the number of each edge and/or each node.
   ************************************************************************************************/
  bool numbers;
  /*!***********************************************************************************************
   * \brief   Construct a \c PlotOptions class object containing default values.
   *
   * Constructs a \c PlotOptions object containing default options and information about the
   * hypergraph and the local solver.
   ************************************************************************************************/
  PlotOptions();
}; // end of class PlotOptions

/*!*************************************************************************************************
 * \brief   Function plotting the solution of an equation on a hypergraph in vtu format.
 *
 * Creates a file according to set plotting options in \c plot_options. This file contains the 
 * solution of the PDE defined in \c plotOpt having the representation \c lambda in terms of its
 * skeleta degrees of freedom (related to skeletal variable lambda).
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
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <class HyperGraphT, class LocalSolverT, typename dof_value_t = double>
void plot
( const HyperGraphT& hyper_graph, const LocalSolverT& local_solver,
  const std::vector<dof_value_t>& lambda, const PlotOptions& plot_options );


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF MEMBER FUNCTIONS OF PlotOptions
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// PlotOptions
// -------------------------------------------------------------------------------------------------

PlotOptions::PlotOptions()
  : outputDir("output"), fileName("example"), fileEnding("vtu"), fileNumber(0),
    printFileNumber(true), incrementFileNumber(true), plot_nodes(false), plot_edges(true),
    n_subintervals(1), scale(1.), numbers(false)
{}


// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF AUXILIARY FUNCTIONS FOR plot() and plot()
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------


// -------------------------------------------------------------------------------------------------
// Auxiliary functions for plot()
// -------------------------------------------------------------------------------------------------


/*!*************************************************************************************************
 * \brief   Auxiliary functions for writing graphics files.
 **************************************************************************************************/
namespace PlotFunctions
{
  /*!***********************************************************************************************
   * \brief   Output of the cubes of the subdivision of an edge in lexicographic order.
   *
   * \tparam  pt_index_t  The index type for global point numbers.
   *
   * \param   dim         The dimension of the cube.
   * \param   n           The number of subdivision points in each direction.
   ************************************************************************************************/
  template <unsigned int dim, typename pt_index_t>
  void  vtu_sub_cube_connectivity(std::ostream& output, unsigned int n, pt_index_t offset)
  {
    if constexpr (dim==1)
      for (unsigned int i=0;i<n-1;++i)  output << offset+i << ' ' << offset+i+1 << "\n";
    else if constexpr (dim==2)
      for (unsigned int i=0;i<n-1;++i)
        for (unsigned int j=0;j<n-1;++j)
          output << offset+i*n+j << ' ' << offset+i*n+j+1 << ' ' << offset+i*n+j+n << ' '
		             << offset+i*n+j+n+1 << "\n";
    else if constexpr (dim==3)
    {
      const unsigned int nn = n*n;
	    for (unsigned int i=0;i<n-1;++i)
        for (unsigned int j=0;j<n-1;++j)
          for (unsigned int k=0;j<n-1;++j)
            output << offset+(i*n+j)*n+k        << ' ' << offset+(i*n+j)*n+k+1      << ' '
                   << offset+(i*n+j)*n+k+n      << ' ' << offset+(i*n+j)*n+k+n+1    << ' '
		               << offset+(i*n+j)*n+k+nn     << ' ' << offset+(i*n+j)*n+k+nn+1   << ' '
		               << offset+(i*n+j)*n+k+nn+n   << ' ' << offset+(i*n+j)*n+k+nn+n+1 << "\n";
    }
  } // end of vtu_sub_cube_connectivity
  /*!***********************************************************************************************
   * \brief   Auxiliary function for writing geometry section VTU files.
   *
   * This function plots the geometry part of an unstructured mesh in* a VTU file. 
   * The typical file structure is
   *
   * ```
   * <Preamble/>
   * <UnstructuredGrid>
   * <Geometry/>
   * <Data/>
   * </UnstructuredGrid>
   *
   * This function writes the `<Geometry>` part of the structure.
   ************************************************************************************************/
  template
  < class HyperGraphT, std::size_t n_subpoints, typename hyEdge_index_t = unsigned int,
    typename pt_index_t = unsigned int >
  void plot_vtu_unstructured_geometry(std::ostream& output,
				      const HyperGraphT& hyper_graph,
				      const std::array<float, n_subpoints>& sub_points,
				      const PlotOptions& plot_options)
  {
    constexpr unsigned int edge_dim = HyperGraphT::hyEdge_dim();
    constexpr unsigned int space_dim = HyperGraphT::space_dim();
    
    const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
    const unsigned int points_per_hyEdge = Hypercube<edge_dim>::pow(n_subpoints);
    
    const pt_index_t n_plot_points = points_per_hyEdge * n_edges;
    const pt_index_t n_plot_edges = n_edges * Hypercube<edge_dim>::pow(n_subpoints-1);
    
    static_assert (edge_dim <= 3);
    unsigned int element_id;
    if constexpr (edge_dim == 1)       element_id = 3;
    else if constexpr (edge_dim == 2)  element_id = 8;
    else if constexpr (edge_dim == 3)  element_id = 11;
    
    output << "    <Piece NumberOfPoints=\"" << n_plot_points << "\" NumberOfCells= \"" 
           << n_plot_edges << "\">" << std::endl;
    output << "      <Points>" << std::endl;
    output << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"
           << std::endl;

    for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
    {
      auto edge = hyper_graph.hyEdge_geometry(he_number);
      auto mapping = edge.mapping_tensor(sub_points);

      for (unsigned int pt_number = 0; pt_number < points_per_hyEdge; ++pt_number)
      {
        output << "        ";
        const Point<space_dim>& point = mapping.lexicographic(pt_number);
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
    output << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
           << std::endl;

    pt_index_t offset = 0;
    for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
    {
      vtu_sub_cube_connectivity<edge_dim>(output, n_subpoints, offset);
      offset += points_per_hyEdge;
    }
    output << "        </DataArray>" << std::endl;
    output << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
    output << "        ";
    
    for (hyEdge_index_t he_number = 1; he_number <= n_plot_edges; ++he_number)
      output << "  " <<  Hypercube<edge_dim>::n_vertices() * he_number;
    output << std::endl;
    
    output << "        </DataArray>" << std::endl;
    output << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
    output << "        ";
    
    for (hyEdge_index_t he_number = 0; he_number < n_plot_edges; ++he_number)
      output << "  " << element_id;
    output << std::endl;
    
    output << "        </DataArray>" << std::endl;
    output << "      </Cells>" << std::endl;
  } // end of void plot_vtu_unstructured_geometry
} // end of namespace PlotFunctions

/*!*************************************************************************************************
 * \brief   Auxiliary function to plot solution values to vtu file.
 **************************************************************************************************/
template
< class HyperGraphT, class LocalSolverT, unsigned int n_subdivisions = 1,
 typename dof_value_t = double, typename hyEdge_index_t = unsigned int >
void plot_vtu
( const HyperGraphT& hyper_graph, const LocalSolverT& local_solver,
	const std::vector<dof_value_t>& lambda, const PlotOptions& plot_options )
{
  constexpr unsigned int edge_dim = HyperGraphT::hyEdge_dim();
  
  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  //  const unsigned int points_per_hyEdge = 1 << edge_dim;
  
  std::array<float, n_subdivisions+1> abscissas;
  for (unsigned int i=0;i<=n_subdivisions;++i)
    abscissas[i] = plot_options.scale*(1.*i/n_subdivisions-0.5)+0.5;
  
  static_assert (edge_dim <= 3 , "Plotting hyperedges with dimensions larger than 3 is hard.");
  
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
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
         << "compressor=\"vtkZLibDataCompressor\">" << std::endl;
  myfile << "  <UnstructuredGrid>" << std::endl;

  PlotFunctions::plot_vtu_unstructured_geometry(myfile, hyper_graph, abscissas, plot_options);

  myfile << "      <CellData>" << std::endl;
  if (plot_options.numbers)
  {
    myfile << "        <DataArray type=\"Int32\" Name=\"index\" NumberOfComponents=\"1\" "
           << "format=\"ascii\">\n";
    if (plot_options.plot_edges)
    {
      for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
        for (unsigned int i=0;i<Hypercube<edge_dim>::pow(n_subdivisions);++i)
          myfile << ' ' << he_number;
    }
    if (plot_options.plot_nodes)  hy_assert(false, "Not yet implemented");
      
    myfile << "        </DataArray>";
  }
  myfile << "      </CellData>" << std::endl;
  
  myfile << "      <PointData>" << std::endl;
  if (LocalSolverT::system_dimension() != 0)
  {
    myfile << "        <DataArray type=\"Float32\" Name=\"values"
	         << "\" NumberOfComponents=\"" << LocalSolverT::system_dimension()
	         << "\" format=\"ascii\">" << std::endl;
      
    std::array< std::array<dof_value_t, HyperGraphT::n_dofs_per_node() > , 2*edge_dim > hyEdge_dofs;
    std::array<unsigned int, 2*edge_dim> hyEdge_hyNodes;
      
    // std::array<std::array<dof_value_t,
    // 			    LocalSolverT::system_dimension()>,
    // 			    Hypercube<edge_dim>::n_vertices()> local_values;
      
    for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
    {
      hyEdge_hyNodes = hyper_graph.hyEdge_topology(he_number).get_hyNode_indices();
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyEdge_dofs[hyNode] = hyper_graph.hyNode_factory().get_dof_values
                                (hyEdge_hyNodes[hyNode], lambda);
      // if constexpr ( LocalSolverT::use_geometry() )
	    // 		 local_values = local_solver.bulk_values
      //                      (abscissas, hyper_graph.hyEdge_geometry(he_number), hyEdge_dofs);
	    // else
      auto local_values = local_solver.bulk_values(abscissas, hyEdge_dofs);
      myfile << "      ";
      for (unsigned int corner = 0; corner < Hypercube<edge_dim>::n_vertices(); ++corner)
	    {
	      myfile << "  ";
	      for (unsigned int d = 0; d < LocalSolverT::system_dimension(); ++d)
          myfile << "  " << local_values[corner][d];
	    }
      myfile << std::endl;
    }
      
    myfile << "        </DataArray>" << std::endl;
  }
  myfile << "      </PointData>" << std::endl;  
  myfile << "    </Piece>" << std::endl;
  myfile << "  </UnstructuredGrid>" << std::endl;
  myfile << "</VTKFile>" << std::endl;
  myfile.close();
} // end of void plot_vtu


// -------------------------------------------------------------------------------------------------
// plot()
// -------------------------------------------------------------------------------------------------

template <class HyperGraphT, class LocalSolverT, typename dof_value_t = double>
void plot
( const HyperGraphT& hyper_graph, const LocalSolverT& local_solver,
	const std::vector<dof_value_t>& lambda, const PlotOptions& plot_options )
{
  hy_assert( plot_options.fileEnding == "vtu" , 
             "Only file ending vtu is supported at the moment. Your choice has been "
             << plot_options.fileEnding << ", which is invalid.");
  hy_assert( !plot_options.fileName.empty() , "File name must not be empty!" );
  hy_assert( !plot_options.outputDir.empty() , "Ouput directory must not be empty!" );
  plot_vtu(hyper_graph, local_solver, lambda, plot_options);
} // end of void plot
