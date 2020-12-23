#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hypercube.hxx>

#ifndef NOFILEOUT
#include <filesystem>
#endif

#include <fstream>
#include <iomanip>
#include <iostream>
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
  std::string outputDir = "output";
  /*!***********************************************************************************************
   * \brief   Name of the file plotted.
   *
   * This \c std::string describes the name of the file for the output plot. Default is "example".
   ************************************************************************************************/
  std::string fileName = "example";
  /*!***********************************************************************************************
   * \brief   Enum for possible file endings.
   ************************************************************************************************/
  enum fileType
  {
    vtu
  };
  /*!***********************************************************************************************
   * \brief   File ending and also way of plotting.
   *
   * This \c std::string describes the name of the file ending for the output plot. Thus, it also
   * characterizes which applications can read the output files properly. Default and currently
   * only option is "vtu" for Paraview.
   ************************************************************************************************/
  fileType fileEnding = vtu;
  /*!***********************************************************************************************
   * \brief   Number of the plot file.
   *
   * This \c unsigned \c int describes the number of the plot file which is created. If a problem is
   * solved repeatedly (e.g. parabolic problem, local refinement, ...) this number indicates the
   * number of the file (e.g. time step, refinement step, ...). Default is 0.
   ************************************************************************************************/
  unsigned int fileNumber = 0;
  /*!***********************************************************************************************
   * \brief   Decide whether \c fileNumber is part of the file name.
   *
   * This \c boolean discriminates whether the \c fileNumber should appear within the name of the
   * file (true) or not (false). Default is true.
   ************************************************************************************************/
  bool printFileNumber = true;
  /*!***********************************************************************************************
   * \brief   Decide whether \c fileNumber is incremented after plotting.
   *
   * This \c boolean discriminates whether the \c fileNumber should be incremented after a file has
   * been written (true) or not (false). Default is true.
   ************************************************************************************************/
  bool incrementFileNumber = true;
  /*!***********************************************************************************************
   * \brief   Include the edge boundaries with their function values into the plot.
   *
   * Defaults to false.
   ************************************************************************************************/
  bool plot_edge_boundaries = false;
  /*!***********************************************************************************************
   * \brief   Include the hyperedges with their function values into the plot.
   *
   * Defaults to true.
   ************************************************************************************************/
  bool plot_edges = true;
  /*!***********************************************************************************************
   * \brief   Number of subintervals for plotting.
   *
   * When plotting an interval, it is split into #n_subintervals intervals such that higher order
   * polynomials and other functions can be displayed appropriately. When plotting higher
   * dimensional objects, this subdivision is applied accordingly in each direction.
   *
   * This functionality is implemented such that higher order polynomials can be output as piecewise
   * linear functions giving them sufficient meaning. It will increase the number of cells seen by
   * the visualization tool, such that the cell boundaries there are not the actual cell boundaries
   * anymore. You can still use the parameter #scale below to see the separate edges.
   *
   * Defaults to 1.
   ************************************************************************************************/
  unsigned int n_subintervals = 1;
  /*!***********************************************************************************************
   * \brief   A factor for scaling each object of the plot locally.
   *
   * This factor defaults to 1 in order to produce a plot of a contiguous domain. If it is chosen
   * less than 1, each edge or node is scaled by this factor around its center.
   ************************************************************************************************/
  float scale = 1.;
  /*!***********************************************************************************************
   * \brief   A factor by which to separate two boundaries belonging to different edges.
   *
   * This factor defaults to 1 which results in boundaries of adjacent edges being drawn on top of
   * each other, if scale is less than one a suitable choice would be (1+scale)/2 in order separate
   * different boundaries. A smaller factor results in a wider gap.
   ************************************************************************************************/
  float boundary_scale = 1.;
  /*!***********************************************************************************************
   * \brief   Output a cell data array with the number of each edge and/or each node.
   ************************************************************************************************/
  bool numbers = false;
};  // end of class PlotOptions

/*!*************************************************************************************************
 * \brief Set a plot option and return the new value of this option as std::string.
 **************************************************************************************************/
std::string set_plot_option(PlotOptions& plot_options,
                            const std::string& option,
                            std::string value = "")
{
  if (value == "")
    ;
  else if (option == "outputDir")
    plot_options.outputDir = value;
  else if (option == "fileName")
    plot_options.fileName = value;
  else if (option == "fileEnding")
  {
    if (value == "vtu")
      plot_options.fileEnding = PlotOptions::vtu;
    else
      hy_assert(false, "You have chosen an invalid file type!");
  }
  else if (option == "fileNumber")
    plot_options.fileNumber = stoi(value);
  else if (option == "printFileNumber")
    plot_options.printFileNumber = (value == "true" || value == "1");
  else if (option == "incrementFileNumber")
    plot_options.incrementFileNumber = (value == "true" || value == "1");
  else if (option == "plotEdges")
    plot_options.plot_edges = (value == "true" || value == "1");
  else if (option == "plotEdgeBoundaries")
    plot_options.plot_edge_boundaries = (value == "true" || value == "1");
  else if (option == "boundaryScale")
    plot_options.boundary_scale = std::stof(value);
  else if (option == "scale")
    plot_options.scale = stof(value);
  else
    hy_assert(false, "This plot option has not been defined (yet).");

  std::string return_value;
  if (option == "outputDir")
    return_value = plot_options.outputDir;
  else if (option == "fileName")
    return_value = plot_options.fileName;
  else if (option == "fileEnding")
    return_value = plot_options.fileEnding;
  else if (option == "fileNumber")
    return_value = std::to_string(plot_options.fileNumber);
  else if (option == "printFileNumber")
    return_value = std::to_string(plot_options.printFileNumber);
  else if (option == "incrementFileNumber")
    return_value = std::to_string(plot_options.incrementFileNumber);
  else if (option == "plotEdges")
    return_value = std::to_string(plot_options.plot_edges);
  else if (option == "plotEdgeBoundaries")
    return_value = std::to_string(plot_options.plot_edge_boundaries);
  else if (option == "scale")
    return_value = std::to_string(plot_options.scale);
  else if (option == "boundaryScale")
    return_value = std::to_string(plot_options.boundary_scale);
  else
    hy_assert(false, "This plot option has not been defined (yet).");

  return return_value;
}
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
 * \tparam  LargeVecT       The typename of the large vector.
 * \tparam  floatT          The floating point type in wihch time is given of dof values.
 * \param   hyper_graph     The hypergraph.
 * \param   local_solver    The local solver.
 * \param   lambda          Large vector containing the skeletal degrees of freedom encoding the
 *                          representation of the unique solution.
 * \param   plot_options    PlotOptions object containing plotting options.
 * \param   time            The time stamp.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename HyperGraphT, typename LocalSolverT, typename LargeVecT, typename floatT = float>
void plot(HyperGraphT& hyper_graph,
          const LocalSolverT& local_solver,
          const LargeVecT& lambda,
          PlotOptions& plot_options,
          const floatT time = 0.);

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------
//
// IMPLEMENTATION OF AUXILIARY FUNCTIONS FOR plot()
//
// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   Auxiliary functions for writing graphics files.
 **************************************************************************************************/
namespace PlotFunctions
{
/*!*************************************************************************************************
 * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
 **************************************************************************************************/
HAS_MEMBER_FUNCTION(bulk_values, has_bulk_values);
/*!*************************************************************************************************
 * \brief   Output of the cubes of the subdivision of an edge in lexicographic order.
 *
 * \tparam  dim         The dimension of the cube.
 * \tparam  pt_index_t  The index type for global point numbers.
 * \param   output      The output.
 * \param   n           The number of subdivision points in each direction.
 * \param   offset      The index of the first vertex
 **************************************************************************************************/
template <unsigned int dim, typename pt_index_t>
void vtu_sub_cube_connectivity(std::ostream& output, unsigned int n, pt_index_t offset)
{
  if constexpr (dim == 0)
    for (unsigned int i = 0; i < n - 1; ++i)
      output << offset + i << "\n";
  else if constexpr (dim == 1)
    for (unsigned int i = 0; i < n - 1; ++i)
      output << offset + i << ' ' << offset + i + 1 << "\n";
  else if constexpr (dim == 2)
    for (unsigned int i = 0; i < n - 1; ++i)
      for (unsigned int j = 0; j < n - 1; ++j)
        output << offset + i * n + j << ' ' << offset + i * n + j + 1 << ' '
               << offset + i * n + j + n << ' ' << offset + i * n + j + n + 1 << "\n";
  else if constexpr (dim == 3)
  {
    const unsigned int nn = n * n;
    for (unsigned int i = 0; i < n - 1; ++i)
      for (unsigned int j = 0; j < n - 1; ++j)
        for (unsigned int k = 0; k < n - 1; ++k)
          output << offset + (i * n + j) * n + k << ' ' << offset + (i * n + j) * n + k + 1 << ' '
                 << offset + (i * n + j) * n + k + n << ' ' << offset + (i * n + j) * n + k + n + 1
                 << ' ' << offset + (i * n + j) * n + k + nn << ' '
                 << offset + (i * n + j) * n + k + nn + 1 << ' '
                 << offset + (i * n + j) * n + k + nn + n << ' '
                 << offset + (i * n + j) * n + k + nn + n + 1 << "\n";
  }
}  // end of vtu_sub_cube_connectivity
/*!*************************************************************************************************
 * \brief   Auxiliary function for writing geometry section VTU files.
 *
 * This function plots the geometry part of an unstructured mesh in* a VTU file.
 * The typical file structure is
 *
 * \< Preamble/\>
 * \< UnstructuredGrid\>
 * \< Geometry/\>
 * \< Data/\>
 * \< /UnstructuredGrid\>
 *
 * This function writes the \< Geometry\> part of the structure.
 **************************************************************************************************/
template <class HyperGraphT,
          unsigned int n_subpoints,
          typename hyEdge_index_t = unsigned int,
          typename pt_index_t = unsigned int>
void plot_vtu_unstructured_geometry(std::ostream& output,
                                    HyperGraphT& hyper_graph,
                                    const SmallVec<n_subpoints, float>& sub_points,
                                    const SmallVec<n_subpoints, float>& boundary_sub_points,
                                    const PlotOptions& plot_options)
{
  constexpr unsigned int edge_dim = HyperGraphT::hyEdge_dim();
  constexpr unsigned int space_dim = HyperGraphT::space_dim();

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  // The number of cells which are actually plotted. This
  // is the number of edges in the graph times the number of
  // cells in an edge due to n_subdivisions
  const pt_index_t n_plot_edges = n_edges * Hypercube<edge_dim>::pow(n_subpoints - 1);
  const unsigned int points_per_edge = Hypercube<edge_dim>::pow(n_subpoints);

  const hyEdge_index_t n_boundaries_per_edge = 2 * edge_dim;
  const hyEdge_index_t n_edge_boundaries = n_edges * n_boundaries_per_edge;
  const unsigned int points_per_boundary = Hypercube<edge_dim - 1>::pow(n_subpoints);
  const pt_index_t n_plot_boundaries =
    n_edge_boundaries * Hypercube<edge_dim - 1>::pow(n_subpoints - 1);

  // The total number of points and cells (using VTK nomenclature)
  // is the sum of the points for edges and the points for
  // edge boundaries. Equally, the total number of cells is the sum of the
  // cells of dimension edge_dim plus those of dimension node_dim.
  const pt_index_t n_plot_points =
    (plot_options.plot_edges ? (points_per_edge * n_edges) : 0) +
    (plot_options.plot_edge_boundaries ? (points_per_boundary * n_edge_boundaries) : 0);
  const pt_index_t n_plot_cells = (plot_options.plot_edges ? n_plot_edges : 0) +
                                  (plot_options.plot_edge_boundaries ? n_plot_boundaries : 0);

  // The element id can be found in the VTK file format documentation.
  static_assert(edge_dim <= 3);
  unsigned int element_id;
  if constexpr (edge_dim == 1)
    element_id = 3;
  else if constexpr (edge_dim == 2)
    element_id = 8;
  else if constexpr (edge_dim == 3)
    element_id = 11;
  unsigned int boundary_element_id = 0;
  if constexpr (edge_dim == 1)
    boundary_element_id = 1;
  else if constexpr (edge_dim == 2)
    boundary_element_id = 3;
  else if constexpr (edge_dim == 3)
    boundary_element_id = 8;

  output << "    <Piece NumberOfPoints=\"" << n_plot_points << "\" NumberOfCells= \""
         << n_plot_cells << "\">" << std::endl;
  output << "      <Points>" << std::endl;
  output << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">"
         << std::endl;
  if (plot_options.plot_edges)
  {
    // For each edge, output the corners of the subdivided cells in lexicographic order.
    for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
    {
      auto edge = hyper_graph.hyEdge_geometry(he_number);

      for (unsigned int pt_number = 0; pt_number < points_per_edge; ++pt_number)
      {
        output << "        ";
        const Point<space_dim> point =
          (Point<space_dim>)edge.template lexicographic<n_subpoints>(pt_number, sub_points);
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          output << "  " << std::fixed << std::scientific << std::setprecision(3) << point[dim];
        for (unsigned int dim = space_dim; dim < 3; ++dim)
          output << "  0.0";
        output << std::endl;
      }
    }
  }
  if (plot_options.plot_edge_boundaries)
  {
    // For each edge accumulate edge boundary coordinates
    for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
    {
      auto edge = hyper_graph.hyEdge_geometry(he_number);
      for (hyEdge_index_t boundary = 0; boundary < n_boundaries_per_edge; ++boundary)
      {
        for (unsigned int pt_number = 0; pt_number < points_per_boundary; ++pt_number)
        {
          output << "        ";
          const Point<space_dim> point =
            (Point<space_dim>)edge.template boundary_lexicographic<n_subpoints>(
              pt_number, boundary, plot_options.boundary_scale, boundary_sub_points);

          if (he_number == 0)
            for (unsigned int dim = 0; dim < space_dim; ++dim)
            {
              output << "  " << std::fixed << std::scientific << std::setprecision(3) << point[dim];
            }
          else
            for (unsigned int dim = 0; dim < space_dim; ++dim)
            {
              output << "  " << std::fixed << std::scientific << std::setprecision(3) << point[dim];
            }
          for (unsigned int dim = space_dim; dim < 3; ++dim)
            output << "  0.0";
          output << std::endl;
        }
      }
    }
  }

  output << "        </DataArray>" << std::endl;
  output << "      </Points>" << std::endl;
  output << "      <Cells>" << std::endl;
  output << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">"
         << std::endl;

  pt_index_t offset = 0;
  if (plot_options.plot_edges)
    for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
    {
      vtu_sub_cube_connectivity<edge_dim>(output, n_subpoints, offset);
      offset += points_per_edge;
    }
  if (plot_options.plot_edge_boundaries)
    for (pt_index_t i = 0; i < n_edge_boundaries; ++i)
    {
      vtu_sub_cube_connectivity<edge_dim - 1>(output, n_subpoints, offset);
      offset += points_per_boundary;
    }

  hy_assert(offset == n_plot_points, "We did not write the right number of connectivity data");

  output << "        </DataArray>" << std::endl;
  output << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << std::endl;
  output << "        ";
  int start_number = 0;
  if (plot_options.plot_edges)
  {
    for (hyEdge_index_t he_number = 1; he_number <= n_plot_edges; ++he_number)
    {
      output << "  " << Hypercube<edge_dim>::n_vertices() * he_number;
      start_number = Hypercube<edge_dim>::n_vertices() * he_number;
    }
  }
  if (plot_options.plot_edge_boundaries)
  {
    for (hyEdge_index_t bdr_number = 1; bdr_number <= n_plot_boundaries; ++bdr_number)
      output << "  " << Hypercube<edge_dim - 1>::n_vertices() * bdr_number + start_number;
  }
  output << std::endl;

  output << "        </DataArray>" << std::endl;
  output << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << std::endl;
  output << "        ";
  if (plot_options.plot_edges)
  {
    for (hyEdge_index_t he_number = 0; he_number < n_plot_edges; ++he_number)
      output << "  " << element_id;
  }
  if (plot_options.plot_edge_boundaries)
  {
    for (hyEdge_index_t bdr_number = 1; bdr_number <= n_plot_boundaries; ++bdr_number)
      output << "  " << boundary_element_id;
  }
  output << std::endl;

  output << "        </DataArray>" << std::endl;
  output << "      </Cells>" << std::endl;
}  // end of void plot_vtu_unstructured_geometry

/*!*************************************************************************************************
 * \brief   Check if some file exists and can be opened
 **************************************************************************************************/
void check_file_opened(const std::ofstream& output_file, const std::string filename)
{
  if (!output_file.is_open())
  {
    throw std::ios_base::failure("File  " + filename + " could not be opened");
  }
}
/*!*************************************************************************************************
 * \brief   Check if an outputstream has opened a file
 **************************************************************************************************/
void create_directory_if_needed(
#ifndef NOFILEOUT
  std::ofstream& output_file,
  const std::string filename,
  const PlotOptions& plot_options)
#else
  std::ofstream&,
  const std::string,
  const PlotOptions&)
#endif
{
#ifndef NOFILEOUT
  try
  {
    check_file_opened(output_file, filename);
  }
  catch (std::ios_base::failure& e)
  {
    std::cerr << e.what() << std::endl;
    std::cout << "Trying to create output directory" << std::endl;
    std::filesystem::create_directory(plot_options.outputDir);
  }
#endif
}
}  // end of namespace PlotFunctions

/*!*************************************************************************************************
 * \brief   Auxiliary function to get dof values of edge.
 **************************************************************************************************/
template <unsigned int edge_dim, class HyperGraphT, typename hyEdge_index_t, typename LargeVecT>
std::array<std::array<typename LargeVecT::value_type, HyperGraphT::n_dofs_per_node()>, 2 * edge_dim>
get_edge_dof_values(HyperGraphT& hyper_graph, hyEdge_index_t edge_index, const LargeVecT& lambda)
{
  std::array<std::array<typename LargeVecT::value_type, HyperGraphT::n_dofs_per_node()>,
             2 * edge_dim>
    hyEdge_dofs;
  SmallVec<2 * edge_dim, unsigned int> hyEdge_hyNodes;
  hyEdge_hyNodes = hyper_graph.hyEdge_topology(edge_index).get_hyNode_indices();
  for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
  {
    hyEdge_dofs[hyNode] = hyper_graph.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode],
                                                                      lambda, hyEdge_dofs[hyNode]);
  }
  return hyEdge_dofs;
}

/*!*************************************************************************************************
 * \brief   Auxiliary function to plot solution values on edge.
 **************************************************************************************************/
template <class HyperGraphT,
          class LocalSolverT,
          typename LargeVecT,
          typename floatT,
          unsigned int n_subdivisions = 1,
          typename hyEdge_index_t = unsigned int>
void plot_edge_values(HyperGraphT& hyper_graph,
                      const LocalSolverT& local_solver,
                      const LargeVecT& lambda,
                      std::ofstream& myfile,
                      SmallVec<n_subdivisions + 1, float> abscissas,
                      const floatT time)
{
  using dof_value_t = typename LargeVecT::value_type;
  constexpr unsigned int edge_dim = HyperGraphT::hyEdge_dim();

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  std::array<std::array<dof_value_t, HyperGraphT::n_dofs_per_node()>, 2 * edge_dim> hyEdge_dofs;

  for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
  {
    hyEdge_dofs = get_edge_dof_values<edge_dim, HyperGraphT, hyEdge_index_t, LargeVecT>(
      hyper_graph, he_number, lambda);
    std::array<
      std::array<dof_value_t, Hypercube<HyperGraphT::hyEdge_dim()>::pow(n_subdivisions + 1)>,
      LocalSolverT::system_dimension()>
      local_values;
    if constexpr (PlotFunctions::has_bulk_values<LocalSolverT,
                                                 decltype(local_values)(decltype(abscissas.data())&,
                                                                        decltype(hyEdge_dofs)&,
                                                                        decltype(time))>::value)
      local_values = local_solver.bulk_values(abscissas.data(), hyEdge_dofs, time);
    else if constexpr (PlotFunctions::has_bulk_values<
                         LocalSolverT, decltype(local_values)(
                                         decltype(abscissas.data())&, decltype(hyEdge_dofs)&,
                                         decltype(hyper_graph[he_number])&, decltype(time))>::value)
    {
      auto geometry = hyper_graph[he_number];
      local_values = local_solver.bulk_values(abscissas.data(), hyEdge_dofs, geometry, time);
    }
    else
      hy_assert(false, "Function seems not to be implemented!");

    myfile << "      ";
    for (unsigned int corner = 0; corner < Hypercube<edge_dim>::n_vertices(); ++corner)
    {
      myfile << "  ";
      for (unsigned int d = 0; d < LocalSolverT::system_dimension(); ++d)
        myfile << "  " << local_values[d][corner];  // AR: I switched d and corner!?
      for (unsigned int d = LocalSolverT::system_dimension();
           d < LocalSolverT::node_system_dimension(); ++d)
        myfile << "  " << 0;  // AR: I switched d and corner!?
    }
    myfile << std::endl;
  }
}
/*!*************************************************************************************************
 * \brief   Auxiliary function to plot solution values on edge boundary.
 **************************************************************************************************/
template <class HyperGraphT,
          class LocalSolverT,
          typename LargeVecT,
          unsigned int n_subdivisions = 1,
          typename hyEdge_index_t = unsigned int>
void plot_boundary_values(HyperGraphT& hyper_graph,
                          const LocalSolverT& local_solver,
                          const LargeVecT& lambda,
                          std::ofstream& myfile,
                          SmallVec<n_subdivisions + 1, float> abscissas)
{
  using dof_value_t = typename LargeVecT::value_type;
  constexpr unsigned int edge_dim = HyperGraphT::hyEdge_dim();

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  std::array<std::array<dof_value_t, HyperGraphT::n_dofs_per_node()>, 2 * edge_dim> hyEdge_dofs;
  std::array<
    std::array<dof_value_t, Hypercube<HyperGraphT::hyEdge_dim() - 1>::pow(n_subdivisions + 1)>,
    LocalSolverT::node_system_dimension()>
    local_values;
  for (hyEdge_index_t edge_index = 0; edge_index < n_edges; ++edge_index)
  {
    hyEdge_dofs = get_edge_dof_values<edge_dim, HyperGraphT, hyEdge_index_t, LargeVecT>(
      hyper_graph, edge_index, lambda);
    for (unsigned int bdr_index = 0; bdr_index < edge_dim * 2; ++bdr_index)
    {
      myfile << "      ";

      local_values = local_solver.lambda_values(abscissas.data(), hyEdge_dofs, bdr_index);
      for (unsigned int corner = 0; corner < Hypercube<edge_dim - 1>::n_vertices(); ++corner)
      {
        myfile << "  ";
        for (unsigned int d = 0; d < LocalSolverT::node_system_dimension(); ++d)
          myfile << "  " << local_values[d][corner];  // AR: I switched d and corner!?
      }
      for (unsigned int d = LocalSolverT::node_system_dimension();
           d < LocalSolverT::system_dimension(); ++d)
        myfile << "  " << 0;  // AR: I switched d and corner!?
      myfile << std::endl;
    }
  }
}

/*!*************************************************************************************************
 * \brief   Auxiliary function to plot solution values to vtu file.
 **************************************************************************************************/
template <class HyperGraphT,
          class LocalSolverT,
          typename LargeVecT,
          typename floatT,
          unsigned int n_subdivisions = 1,
          typename hyEdge_index_t = unsigned int>
void plot_vtu(
#ifndef NOFILEOUT
  HyperGraphT& hyper_graph,
  const LocalSolverT& local_solver,
  const LargeVecT& lambda,
  const PlotOptions& plot_options,
  const floatT time = 0.)
#else
  HyperGraphT&,
  const LocalSolverT&,
  const LargeVecT&,
  const PlotOptions&,
  const floatT = 0.)
#endif
{
#ifndef NOFILEOUT
  constexpr unsigned int edge_dim = HyperGraphT::hyEdge_dim();

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  const hyEdge_index_t n_edge_boundaries = n_edges * 2 * edge_dim;
  //  const unsigned int n_points_per_edge = 1 << edge_dim;

  SmallVec<n_subdivisions + 1, float> boundary_abscissas;
  for (unsigned int i = 0; i <= n_subdivisions; ++i)
    boundary_abscissas[i] =
      plot_options.scale * plot_options.boundary_scale * (1. * i / n_subdivisions - 0.5) + 0.5;
  SmallVec<n_subdivisions + 1, float> abscissas;
  for (unsigned int i = 0; i <= n_subdivisions; ++i)
    abscissas[i] = plot_options.scale * (1. * i / n_subdivisions - 0.5) + 0.5;

  static_assert(edge_dim <= 3, "Plotting hyperedges with dimensions larger than 3 is hard.");
  std::ofstream myfile;

  std::string filename = plot_options.outputDir;
  filename.append("/");
  filename.append(plot_options.fileName);
  if (plot_options.printFileNumber)
  {
    filename.append(".");
    filename.append(std::to_string(plot_options.fileNumber));
  }
  filename.append(".vtu");

  myfile.open(filename);
  PlotFunctions::create_directory_if_needed(myfile, filename, plot_options);
  myfile << "<?xml version=\"1.0\"?>" << std::endl;
  myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" "
         << "compressor=\"vtkZLibDataCompressor\">" << std::endl;
  myfile << "  <UnstructuredGrid>" << std::endl;

  PlotFunctions::plot_vtu_unstructured_geometry(myfile, hyper_graph, abscissas, boundary_abscissas,
                                                plot_options);

  myfile << "      <CellData>" << std::endl;
  if (plot_options.numbers)
  {
    myfile << "        <DataArray type=\"Int32\" Name=\"index\" NumberOfComponents=\"1\" "
           << "format=\"ascii\">\n";
    if (plot_options.plot_edges)
    {
      for (hyEdge_index_t he_number = 0; he_number < n_edges; ++he_number)
        for (unsigned int i = 0; i < Hypercube<edge_dim>::pow(n_subdivisions); ++i)
          myfile << ' ' << he_number;
    }
    if (plot_options.plot_edge_boundaries)
    {
      for (hyEdge_index_t bdr_number = 1; bdr_number <= n_edge_boundaries; ++bdr_number)
        for (unsigned int i = 0; i < Hypercube<edge_dim>::pow(n_subdivisions); ++i)
          myfile << "  " << bdr_number;
    }

    myfile << "        </DataArray>";
  }
  myfile << "      </CellData>" << std::endl;

  myfile << "      <PointData>" << std::endl;
  if (LocalSolverT::system_dimension() != 0)
  {
    myfile << "        <DataArray type=\"Float32\" Name=\"values"
           << "\" NumberOfComponents=\""
           << std::max(LocalSolverT::system_dimension(), LocalSolverT::node_system_dimension())
           << "\" format=\"ascii\">" << std::endl;
    if (plot_options.plot_edges)
    {
      plot_edge_values<HyperGraphT, LocalSolverT, LargeVecT, floatT, n_subdivisions,
                       hyEdge_index_t>(hyper_graph, local_solver, lambda, myfile, abscissas, time);
    }
    if (plot_options.plot_edge_boundaries)
    {
      plot_boundary_values<HyperGraphT, LocalSolverT, LargeVecT, n_subdivisions, hyEdge_index_t>(
        hyper_graph, local_solver, lambda, myfile, boundary_abscissas);
    }
    myfile << "        </DataArray>" << std::endl;
  }
  myfile << "      </PointData>" << std::endl;
  myfile << "    </Piece>" << std::endl;
  myfile << "  </UnstructuredGrid>" << std::endl;
  myfile << "</VTKFile>" << std::endl;
  std::cout << plot_options.fileName << " was written\n";
  myfile.close();
#endif
}  // end of void plot_vtu

// -------------------------------------------------------------------------------------------------
// plot()
// -------------------------------------------------------------------------------------------------

template <typename HyperGraphT, typename LocalSolverT, typename LargeVecT, typename floatT>
void plot(HyperGraphT& hyper_graph,
          const LocalSolverT& local_solver,
          const LargeVecT& lambda,
          PlotOptions& plot_options,
          const floatT time)
{
  hy_assert(plot_options.fileEnding == PlotOptions::vtu,
            "Only file ending vtu is supported at the moment. Your choice has been "
              << plot_options.fileEnding << ", which is invalid.");
  hy_assert(!plot_options.fileName.empty(), "File name must not be empty!");
  hy_assert(!plot_options.outputDir.empty(), "Ouput directory must not be empty!");
  plot_vtu(hyper_graph, local_solver, lambda, plot_options, time);
  if (plot_options.incrementFileNumber)
    ++plot_options.fileNumber;
}  // end of void plot
