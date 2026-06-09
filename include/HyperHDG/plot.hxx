#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hypercube.hxx>

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <cstring>

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
    vtu,
    vtkhdf,
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
  /*!***********************************************************************************************
   * \brief   Append per-edge energy components to a CellData/energies dataset (vtkhdf only).
   *
   * Requires the local solver to expose \c n_energy_components() and \c energy(). Defaults to false.
   ************************************************************************************************/
  bool energy = false;
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
    if (value == "vtkhdf")
      plot_options.fileEnding = PlotOptions::vtkhdf;
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
  else if (option == "energy")
    plot_options.energy = (value == "true" || value == "1");
  // else if (option == "n_subintervals")
  //   plot_options.n_subintervals = stoi(value);
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
  else if (option == "energy")
    return_value = std::to_string(plot_options.energy);
  // else if (option == "n_subintervals")
  //   return_value = std::to_string(plot_options.n_subintervals);
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
 * \brief   Prepare struct to check for energy() to exist (cf. compile_time_tricks.hxx).
 **************************************************************************************************/
HAS_MEMBER_FUNCTION(energy, has_energy);
/*!*************************************************************************************************
 * \brief   Prepare struct to check for n_energy_components() to exist.
 **************************************************************************************************/
HAS_MEMBER_FUNCTION(n_energy_components, has_n_energy_components);
/*!*************************************************************************************************
 * \brief   Turn fileType enum into string.
 **************************************************************************************************/
std::string fileType_to_string(const PlotOptions::fileType& type)
{
  switch (type)
  {
    case PlotOptions::fileType::vtu:
      return "vtu";
    case PlotOptions::fileType::vtkhdf:
      return "vtkhdf";
  }
  hy_assert(false, "File type seems to be invalid.");
  return "";
}
/*!*************************************************************************************************
 * \brief   Open stream to file.
 **************************************************************************************************/
std::ofstream open_ofstream(const PlotOptions& plot_options, const bool append = false)
{
  std::ofstream myfile;

  std::string filename = plot_options.outputDir + "/" + plot_options.fileName;
  if (plot_options.printFileNumber)
    filename += "." + std::to_string(plot_options.fileNumber);
  filename += "." + fileType_to_string(plot_options.fileEnding);
  if (std::filesystem::create_directory(plot_options.outputDir))
    std::cout << "Directory \"" << plot_options.outputDir << "\" has been created." << std::endl;

  if (append)
    myfile.open(filename, std::ios_base::app);
  else
    myfile.open(filename, std::ios_base::out);
  hy_assert(myfile.is_open(),
            "The file has not been created. Most likely, the filesystem could not"
              << " create the output directory, since std::filesystem has not been available."
              << std::endl
              << "Please, try to create the output directoy manually and run the code again.");

  return myfile;
}
/*!*************************************************************************************************
 * \brief   Close stream to file.
 **************************************************************************************************/
void close_ofstream(std::ofstream& myfile)
{
  myfile.close();
}

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
template <unsigned int index, typename functions>
static constexpr unsigned int first_dof()
{
  static_assert(index <= std::tuple_size<functions>(), "Index is too large!");
  if constexpr (index == 0)
    return 0;
  else
    return std::tuple_element<index - 1, functions>::n_fun() + first_dof<index - 1>();
}
/*!*************************************************************************************************
 * \brief   Auxiliary function to plot solution values on edge boundary.
 **************************************************************************************************/
template <unsigned int component,
          typename LocalSolverT,
          typename dof_value_t,
          typename lv_t,
          typename hd_t,
          typename pt_t>
void fancy_recursion(__attribute__((unused)) lv_t& local_values,
                     __attribute__((unused)) const hd_t& hyEdge_dofs,
                     __attribute__((unused)) const pt_t& point,
                     __attribute__((unused)) const unsigned int k,
                     __attribute__((unused)) const unsigned int bdr_index)
{
  // CODE THROWS COMPILETIME ERRORS!
  // if constexpr (component == LocalSolverT::node_system_dimension())
  //   return;
  // else
  // {
  //   std::array<
  //     dof_value_t,
  //     std::tuple_element<component, typename
  //     LocalSolverT::node_element::functions>::type::n_fun()> helper_arr;
  //   for (unsigned int k = 0; k < helper_arr.size(); ++k)
  //     helper_arr[k] =
  //       hyEdge_dofs[bdr_index]
  //                  [first_dof<component, typename LocalSolverT::node_element::functions>() + k];

  //   local_values[component][k] =
  //     std::tuple_element<component, typename LocalSolverT::node_element::functions>::type::
  //       template lin_comb_fct_val<float>(SmallVec<helper_arr.size(), dof_value_t>(helper_arr),
  //                                        point);
  //   fancy_recursion<component + 1, LocalSolverT, dof_value_t>(local_values, hyEdge_dofs, point,
  //   k,
  //                                                             bdr_index);
  // }
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
                          const LargeVecT& lambda,
                          std::ofstream& myfile,
                          SmallVec<n_subdivisions + 1, float> abscissas)
{
  using dof_value_t = typename LargeVecT::value_type;
  constexpr unsigned int hyEdge_dim = HyperGraphT::hyEdge_dim();

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  std::array<std::array<dof_value_t, HyperGraphT::n_dofs_per_node()>, 2 * hyEdge_dim> hyEdge_dofs;
  std::array<
    std::array<dof_value_t, Hypercube<HyperGraphT::hyEdge_dim() - 1>::pow(n_subdivisions + 1)>,
    LocalSolverT::node_system_dimension()>
    local_values;
  for (unsigned int i = 0; i < local_values.size(); ++i)
    local_values[i].fill(0.);

  for (hyEdge_index_t edge_index = 0; edge_index < n_edges; ++edge_index)
  {
    hyEdge_dofs = get_edge_dof_values<hyEdge_dim, HyperGraphT, hyEdge_index_t, LargeVecT>(
      hyper_graph, edge_index, lambda);
    for (unsigned int bdr_index = 0; bdr_index < hyEdge_dim * 2; ++bdr_index)
    {
      myfile << "      ";

      for (unsigned int lval = 0; lval < local_values.size(); ++lval)
        fancy_recursion<0, LocalSolverT, dof_value_t>(
          local_values, hyEdge_dofs,
          Hypercube<hyEdge_dim - 1>::template tensorial_pt<Point<hyEdge_dim - 1> >(lval, abscissas),
          lval, bdr_index);
      for (unsigned int corner = 0; corner < Hypercube<hyEdge_dim - 1>::n_vertices(); ++corner)
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
void plot_vtu(HyperGraphT& hyper_graph,
              const LocalSolverT& local_solver,
              const LargeVecT& lambda,
              const PlotOptions& plot_options,
              const floatT time = 0.)
{
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

  std::ofstream myfile = PlotFunctions::open_ofstream(plot_options, false);

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
    myfile << "        <DataArray type=\"Float32\" Name=\"values" << "\" NumberOfComponents=\""
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
        hyper_graph, lambda, myfile, boundary_abscissas);
    }
    myfile << "        </DataArray>" << std::endl;
  }
  myfile << "      </PointData>" << std::endl;
  myfile << "    </Piece>" << std::endl;
  myfile << "  </UnstructuredGrid>" << std::endl;
  myfile << "</VTKFile>" << std::endl;
  PlotFunctions::close_ofstream(myfile);
}  // end of void plot_vtu

#ifdef HYPERHDG_PETSC
#include <hdf5.h>
#include <mpi.h>

// Append n_new_rows to axis 0 of an existing extendable dataset, then write the data.
// For 1D datasets pass n_cols = 1 (ignored). For 2D, axis-1 must match the dataset.
inline void h5_append(hid_t loc, const char* name,
                      hid_t mem_type,
                      hsize_t n_new_rows, hsize_t n_cols,
                      const void* data)
{
  hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
  hy_check(dset >= 0, "h5_append: cannot open '" << name << "'");

  hid_t fspace = H5Dget_space(dset);
  int rank = H5Sget_simple_extent_ndims(fspace);
  hsize_t cur[2] = {0, 0};
  H5Sget_simple_extent_dims(fspace, cur, nullptr);
  H5Sclose(fspace);

  hsize_t new_dims[2] = {cur[0] + n_new_rows, (rank == 2 ? n_cols : 0)};
  H5Dset_extent(dset, new_dims);

  fspace = H5Dget_space(dset);
  hsize_t start[2] = {cur[0], 0};
  hsize_t count[2] = {n_new_rows, (rank == 2 ? n_cols : 0)};
  H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, nullptr, count, nullptr);

  hsize_t mem_dims[2] = {n_new_rows, (rank == 2 ? n_cols : 0)};
  hid_t mspace = H5Screate_simple(rank, mem_dims, nullptr);

  H5Dwrite(dset, mem_type, mspace, fspace, H5P_DEFAULT, data);

  H5Sclose(mspace); H5Sclose(fspace); H5Dclose(dset);
}

// Overwrite element 0 of an existing 1D dataset (used to keep the single-part NumberOf* counts as
// running totals while ranks append their pieces in a token ring).
inline void h5_write_i64_elem0(hid_t loc, const char* name, int64_t value)
{
  hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
  hy_check(dset >= 0, "h5_write_i64_elem0: cannot open '" << name << "'");
  hid_t fspace = H5Dget_space(dset);
  hsize_t start = 0, count = 1;
  H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &start, nullptr, &count, nullptr);
  hid_t mspace = H5Screate_simple(1, &count, nullptr);
  H5Dwrite(dset, H5T_NATIVE_INT64, mspace, fspace, H5P_DEFAULT, &value);
  H5Sclose(mspace); H5Sclose(fspace); H5Dclose(dset);
}

// Update a scalar int64 attribute on a group (delete + recreate; HDF5 attrs are not extendable).
inline void h5_set_attr_i64(hid_t loc, const char* name, int64_t value)
{
  if (H5Aexists(loc, name) > 0)
    H5Adelete(loc, name);
  hid_t space = H5Screate(H5S_SCALAR);
  hid_t attr  = H5Acreate2(loc, name, H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(attr, H5T_NATIVE_INT64, &value);
  H5Aclose(attr); H5Sclose(space);
}

// =================================================================================================
// plot_vtkhdf_mesh: write the static mesh (geometry, topology, NumberOf* arrays, types_points).
// Truncates the file.
// =================================================================================================
template <class HyperGraphT,
          unsigned int n_subdivisions = 1,
          typename hyEdge_index_t = unsigned int,
          typename pt_index_t = unsigned int>
void plot_vtkhdf_mesh(HyperGraphT& hyper_graph,
                      const PlotOptions& plot_options, unsigned int n_components,
                      unsigned int n_energy_components = 0,
                      bool append = false,
                      pt_index_t point_base = 0,
                      pt_index_t cell_base = 0,
                      pt_index_t conn_base = 0)
{
  constexpr unsigned int edge_dim  = HyperGraphT::hyEdge_dim();
  constexpr unsigned int space_dim = HyperGraphT::space_dim();
  static_assert(edge_dim <= 3, "Plotting hyperedges with dim > 3 is hard.");

  constexpr unsigned int n_subpoints     = n_subdivisions + 1;
  constexpr unsigned int points_per_edge = Hypercube<edge_dim>::pow(n_subpoints);
  constexpr unsigned int cells_per_edge  = Hypercube<edge_dim>::pow(n_subdivisions);
  constexpr unsigned int verts_per_cell  = Hypercube<edge_dim>::n_vertices();

  // VTK cell type id (same logic as plot_vtu)
  uint8_t element_id = 0;
  if constexpr (edge_dim == 1) element_id = 3;   // VTK_LINE
  else if constexpr (edge_dim == 2) element_id = 8;
  else if constexpr (edge_dim == 3) element_id = 11;

  // Abscissas (same as plot_vtu, no boundary scaling since we don't do boundaries)
  SmallVec<n_subpoints, float> abscissas;
  for (unsigned int i = 0; i < n_subpoints; ++i)
    abscissas[i] = plot_options.scale * (1.f * i / n_subdivisions - 0.5f) + 0.5f;

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  const pt_index_t n_points    = static_cast<pt_index_t>(n_edges) * points_per_edge;
  const pt_index_t n_cells     = static_cast<pt_index_t>(n_edges) * cells_per_edge;
  const pt_index_t n_conn      = n_cells * verts_per_cell;

  // -----------------------------------------------------------------------
  // Build buffers
  // -----------------------------------------------------------------------

  // Points: (n_points, 3), row-major, padded to 3D
  std::vector<float> points(3 * n_points, 0.f);
  for (hyEdge_index_t he = 0; he < n_edges; ++he) {
    auto edge = hyper_graph.hyEdge_geometry(he);
    for (unsigned int p = 0; p < points_per_edge; ++p) {
      const Point<space_dim> pt =
        (Point<space_dim>)edge.template lexicographic<n_subpoints>(p, abscissas);
      const pt_index_t row = he * points_per_edge + p;
      for (unsigned int d = 0; d < space_dim; ++d)
        points[3 * row + d] = pt[d];
      // dims [space_dim, 3) already zero
    }
  }

  // Connectivity: lexicographic sub-cube vertices per edge.
  // We replicate vtu_sub_cube_connectivity inline here to fill a buffer
  // instead of streaming text. Only handle 1D for now; 2D/3D can follow
  // the same pattern as plot_vtu if needed.
  std::vector<int64_t> connectivity;
  connectivity.reserve(n_conn);
  for (hyEdge_index_t he = 0; he < n_edges; ++he) {
    const pt_index_t offset = he * points_per_edge;
    if constexpr (edge_dim == 1)
      for (unsigned int i = 0; i < n_subdivisions; ++i) {
        connectivity.push_back(offset + i);
        connectivity.push_back(offset + i + 1);
      }
    else
      hy_check(edge_dim == 1, "expected edge_dim == 1 found edge_dim == " << edge_dim);
  }

  // Offsets: length n_cells + 1, starts at 0, ends at n_conn (VTKHDF convention)
  std::vector<int64_t> offsets(n_cells + 1);
  for (pt_index_t i = 0; i <= n_cells; ++i)
    offsets[i] = static_cast<int64_t>(i * verts_per_cell);

  // Types: one per cell
  std::vector<uint8_t> types(n_cells, element_id);

  // Per-piece counts (single piece = length-1 arrays)
  const int64_t np  = n_points;
  const int64_t nc  = n_cells;
  const int64_t nci = n_conn;

  // -----------------------------------------------------------------------
  // Distributed append: this rank's edges are concatenated into the single VTKHDF part written by
  // rank 0. Connectivity is shifted into this rank's point range, offsets into the running
  // connectivity base (dropping the leading 0), and the NumberOf* counts bumped to running totals.
  // Called in a token ring so only one rank touches the file at a time.
  // -----------------------------------------------------------------------
  if (append)
  {
    for (auto& c : connectivity)
      c += static_cast<int64_t>(point_base);

    hid_t file = H5Fopen(plot_options.fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hy_check(file >= 0, "failed to open HDF5 file '" << plot_options.fileName << "'");
    hid_t root = H5Gopen2(file, "VTKHDF", H5P_DEFAULT);

    h5_append(root, "Points",       H5T_NATIVE_FLOAT, (hsize_t)n_points, 3, points.data());
    h5_append(root, "Connectivity", H5T_NATIVE_INT64, (hsize_t)n_conn, 1, connectivity.data());
    {
      std::vector<int64_t> off_shifted(n_cells);
      for (pt_index_t i = 0; i < n_cells; ++i)
        off_shifted[i] = offsets[i + 1] + static_cast<int64_t>(conn_base);
      h5_append(root, "Offsets", H5T_NATIVE_INT64, (hsize_t)n_cells, 1, off_shifted.data());
    }
    h5_append(root, "Types", H5T_NATIVE_UINT8, (hsize_t)n_cells, 1, types.data());

    h5_write_i64_elem0(root, "NumberOfPoints",          static_cast<int64_t>(point_base) + np);
    h5_write_i64_elem0(root, "NumberOfCells",           static_cast<int64_t>(cell_base) + nc);
    h5_write_i64_elem0(root, "NumberOfConnectivityIds", static_cast<int64_t>(conn_base) + nci);

    if constexpr (n_subdivisions == 1 && edge_dim == 1)
    {
      std::vector<int32_t> node_types(n_points);
      for (hyEdge_index_t he = 0; he < n_edges; ++he)
      {
        auto edge = hyper_graph[he];
        const pt_index_t base = he * points_per_edge;
        node_types[base + 0] = static_cast<int32_t>(edge.node_descriptor[0]);
        node_types[base + 1] = static_cast<int32_t>(edge.node_descriptor[1]);
      }
      hid_t pdata = H5Gopen2(root, "PointData", H5P_DEFAULT);
      h5_append(pdata, "types_points", H5T_NATIVE_INT32, (hsize_t)n_points, 1, node_types.data());
      H5Gclose(pdata);
    }

    if (n_edges > 0 && hyper_graph.hyEdge_geometry(0).has_extra_data())
    {
      const unsigned int n_properties =
        static_cast<unsigned int>(hyper_graph.hyEdge_geometry(0).extra_data().size());
      std::vector<float> props_buf(static_cast<size_t>(n_cells) * n_properties, 0.f);
      for (hyEdge_index_t he = 0; he < n_edges; ++he)
      {
        const auto& geom = hyper_graph.hyEdge_geometry(he);
        const auto& props = geom.extra_data();
        for (unsigned int c = 0; c < cells_per_edge; ++c)
        {
          const size_t row = (static_cast<size_t>(he) * cells_per_edge + c) * n_properties;
          for (unsigned int d = 0; d < n_properties; ++d)
            props_buf[row + d] = static_cast<float>(props[d]);
        }
      }
      hid_t cdata = H5Gopen2(root, "CellData", H5P_DEFAULT);
      h5_append(cdata, "properties", H5T_NATIVE_FLOAT, (hsize_t)n_cells, n_properties,
                props_buf.data());
      H5Gclose(cdata);
    }

    H5Gclose(root);
    H5Fclose(file);
    return;
  }

  // -----------------------------------------------------------------------
  // Write file
  // -----------------------------------------------------------------------

  // H5F_ACC_TRUNC - truncate+R/W, i.e. discard all file contents before writing
  hid_t file = H5Fcreate(plot_options.fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  hy_check(file >= 0, "failed to create HDF5 file '" << plot_options.fileName << "'");

  hid_t root = H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // --- attribute: Version = [2, 0]
  {
    hsize_t dim = 2;
    hid_t space = H5Screate_simple(1, &dim, nullptr);
    hid_t attr  = H5Acreate2(root, "Version", H5T_STD_I64LE, space,
                             H5P_DEFAULT, H5P_DEFAULT);
    int64_t version[2] = {2, 0};
    H5Awrite(attr, H5T_NATIVE_INT64, version);
    H5Aclose(attr); H5Sclose(space);
  }
  // --- attribute: Type = "UnstructuredGrid"
  {
    const char* s = "UnstructuredGrid";
    hid_t stype = H5Tcopy(H5T_C_S1);
    H5Tset_size(stype, std::strlen(s));
    H5Tset_strpad(stype, H5T_STR_NULLPAD);
    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr  = H5Acreate2(root, "Type", stype, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, stype, s);
    H5Aclose(attr); H5Sclose(space); H5Tclose(stype);
  }

  // --- helper to write a chunked + unlimited dataset (whole buffer)
  auto write_dset = [&](hid_t loc, const char* name,
                        hid_t file_type, hid_t mem_type,
                        int rank, const hsize_t* dims, const void* data)
  {
    std::vector<hsize_t> maxdims(rank), chunk(rank);
    for (int i = 0; i < rank; ++i)
    {
      maxdims[i] = (i == 0) ? H5S_UNLIMITED : dims[i];
      chunk[i]   = (i == 0) ? std::max<hsize_t>(dims[i], 1) : dims[i];
    }
    hid_t space = H5Screate_simple(rank, dims, maxdims.data());
    hid_t dcpl  = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, rank, chunk.data());
    hid_t dset = H5Dcreate2(loc, name, file_type, space,
                            H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dwrite(dset, mem_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Dclose(dset); H5Pclose(dcpl); H5Sclose(space);
  };

  // --- geometry
  { hsize_t d[2] = {(hsize_t)n_points, 3};
    write_dset(root, "Points", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, 2, d, points.data()); }
  { hsize_t d = (hsize_t)n_conn;
    write_dset(root, "Connectivity", H5T_STD_I64LE, H5T_NATIVE_INT64, 1, &d, connectivity.data()); }
  { hsize_t d = (hsize_t)offsets.size();
    write_dset(root, "Offsets", H5T_STD_I64LE, H5T_NATIVE_INT64, 1, &d, offsets.data()); }
  { hsize_t d = (hsize_t)n_cells;
    write_dset(root, "Types", H5T_STD_U8LE, H5T_NATIVE_UINT8, 1, &d, types.data()); }

  // --- per-piece counts (single piece → length 1)
  { hsize_t d = 1;
    write_dset(root, "NumberOfPoints",          H5T_STD_I64LE, H5T_NATIVE_INT64, 1, &d, &np); }
  { hsize_t d = 1;
    write_dset(root, "NumberOfCells",           H5T_STD_I64LE, H5T_NATIVE_INT64, 1, &d, &nc); }
  { hsize_t d = 1;
    write_dset(root, "NumberOfConnectivityIds", H5T_STD_I64LE, H5T_NATIVE_INT64, 1, &d, &nci); }

  // --- point data  group
  hid_t pdata = H5Gcreate2(root, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // --- types points
  if constexpr (n_subdivisions == 1 && edge_dim == 1) {
    std::vector<int32_t> node_types(n_points);
    for (hyEdge_index_t he = 0; he < n_edges; ++he) {
      auto edge = hyper_graph[he];
      const pt_index_t base = he * points_per_edge;
      node_types[base+0] = static_cast<int32_t>(edge.node_descriptor[0]);
      node_types[base+1] = static_cast<int32_t>(edge.node_descriptor[1]);
    }
    hsize_t d = (hsize_t)n_points;
    write_dset(pdata, "types_points", H5T_STD_I32LE, H5T_NATIVE_INT32, 1, &d, node_types.data());
  }

  // --- empty extendable PointData/values: (0, n_components), unlimited axis 0
  {
    hsize_t dims[2]    = {0, n_components};
    hsize_t maxdims[2] = {H5S_UNLIMITED, n_components};
    hsize_t chunk[2]   = {std::max<hsize_t>(n_points, 1), n_components};
    hid_t space = H5Screate_simple(2, dims, maxdims);
    hid_t dcpl  = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 2, chunk);
    hid_t dset = H5Dcreate2(pdata, "values", H5T_IEEE_F32LE, space,
                            H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dclose(dset); H5Pclose(dcpl); H5Sclose(space);
  }

  // --- static CellData/properties from hyper_edge.geometry.extra_data()
  hid_t cdata = H5Gcreate2(root, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (n_edges > 0 && hyper_graph.hyEdge_geometry(0).has_extra_data()) {
    const unsigned int n_properties =
      static_cast<unsigned int>(hyper_graph.hyEdge_geometry(0).extra_data().size());

    std::vector<float> props_buf(static_cast<size_t>(n_cells) * n_properties, 0.f);
    for (hyEdge_index_t he = 0; he < n_edges; ++he) {
      auto edge = hyper_graph.hyEdge_geometry(he);
      const auto& props = edge.extra_data();
      hy_check(props.size() == n_properties,
               "all hyperedges must have the same number of properties; "
               "edge " << he << " has " << props.size() << ", expected " << n_properties);
      for (unsigned int c = 0; c < cells_per_edge; ++c) {
        const size_t row = (static_cast<size_t>(he) * cells_per_edge + c) * n_properties;
        for (unsigned int d = 0; d < n_properties; ++d)
          props_buf[row + d] = static_cast<float>(props[d]);
      }
    }

    hsize_t d[2] = {(hsize_t)n_cells, n_properties};
    write_dset(cdata, "properties", H5T_IEEE_F32LE, H5T_NATIVE_FLOAT, 2, d, props_buf.data());
  }

  // --- empty extendable CellData/energies: (0, n_energy_components), unlimited axis 0
  if (n_energy_components > 0) {
    hsize_t dims[2]    = {0, n_energy_components};
    hsize_t maxdims[2] = {H5S_UNLIMITED, n_energy_components};
    hsize_t chunk[2]   = {std::max<hsize_t>(n_cells, 1), n_energy_components};
    hid_t space = H5Screate_simple(2, dims, maxdims);
    hid_t dcpl  = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 2, chunk);
    hid_t dset = H5Dcreate2(cdata, "energies", H5T_IEEE_F32LE, space,
                            H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dclose(dset); H5Pclose(dcpl); H5Sclose(space);
  }
  H5Gclose(cdata);

  // --- Steps group + NSteps=0 attribute
  hid_t steps = H5Gcreate2(root, "Steps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5_set_attr_i64(steps, "NSteps", 0);

  // helper for empty 1D extendable datasets
  auto make_empty_1d = [&](hid_t loc, const char* name, hid_t file_type) {
    hsize_t dims    = 0;
    hsize_t maxdims = H5S_UNLIMITED;
    hsize_t chunk   = 1;
    hid_t space = H5Screate_simple(1, &dims, &maxdims);
    hid_t dcpl  = H5Pcreate(H5P_DATASET_CREATE);
    H5Pset_chunk(dcpl, 1, &chunk);
    hid_t dset = H5Dcreate2(loc, name, file_type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
    H5Dclose(dset); H5Pclose(dcpl); H5Sclose(space);
  };

  make_empty_1d(steps, "Values",                H5T_IEEE_F64LE);
  make_empty_1d(steps, "PartOffsets",           H5T_STD_I64LE);
  make_empty_1d(steps, "PointOffsets",          H5T_STD_I64LE);
  make_empty_1d(steps, "CellOffsets",           H5T_STD_I64LE);
  make_empty_1d(steps, "ConnectivityIdOffsets", H5T_STD_I64LE);
  make_empty_1d(steps, "NumberOfParts",         H5T_STD_I64LE);

  hid_t pdo = H5Gcreate2(steps, "PointDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  make_empty_1d(pdo, "values", H5T_STD_I64LE);
  H5Gclose(pdo);

  if (n_energy_components > 0) {
    hid_t cdo = H5Gcreate2(steps, "CellDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    make_empty_1d(cdo, "energies", H5T_STD_I64LE);
    H5Gclose(cdo);
  }

  H5Gclose(steps);
  H5Gclose(pdata);
  H5Gclose(root);
  H5Fclose(file);
}

// =================================================================================================
// plot_vtkhdf_bulk: open existing file R/W, write PointData/values for this step.
// (Currently overwrites — temporal append comes next.)
// =================================================================================================
template <class HyperGraphT,
          class LocalSolverT,
          typename LargeVecT,
          typename floatT,
          unsigned int n_subdivisions = 1,
          typename hyEdge_index_t = unsigned int,
          typename pt_index_t = unsigned int>
void plot_vtkhdf_bulk(HyperGraphT& hyper_graph,
                      const LocalSolverT& local_solver,
                      const LargeVecT& lambda,
                      const PlotOptions& plot_options,
                      const floatT time = 0.,
                      bool append = false)
{
  constexpr unsigned int edge_dim  = HyperGraphT::hyEdge_dim();
  static_assert(edge_dim <= 3, "Plotting hyperedges with dim > 3 is hard.");

  constexpr unsigned int n_subpoints     = n_subdivisions + 1;
  constexpr unsigned int points_per_edge = Hypercube<edge_dim>::pow(n_subpoints);
  constexpr unsigned int cells_per_edge  = Hypercube<edge_dim>::pow(n_subdivisions);

  SmallVec<n_subpoints, float> abscissas;
  for (unsigned int i = 0; i < n_subpoints; ++i)
    abscissas[i] = plot_options.scale * (1.f * i / n_subdivisions - 0.5f) + 0.5f;

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  const pt_index_t n_points    = static_cast<pt_index_t>(n_edges) * points_per_edge;
  const pt_index_t n_cells     = static_cast<pt_index_t>(n_edges) * cells_per_edge;

  if constexpr (LocalSolverT::system_dimension() == 0) return;

  using dof_value_t = typename LargeVecT::value_type;
  constexpr unsigned int n_components = LocalSolverT::system_dimension();

  std::vector<float> values(static_cast<size_t>(n_points) * n_components, 0.f);

  // Energy buffer: only populated when plot_options.energy and LocalSolverT supplies energy().
  std::vector<float> energies_buf;
  constexpr bool has_energy_api =
    PlotFunctions::has_n_energy_components<LocalSolverT, unsigned int()>::value;

  std::array<std::array<dof_value_t, HyperGraphT::n_dofs_per_node()>, 2 * edge_dim>
    hyEdge_dofs;

  if (plot_options.energy) {
    if constexpr (has_energy_api)
      energies_buf.assign(static_cast<size_t>(n_cells) * LocalSolverT::n_energy_components(), 0.f);
    else
      hy_check(false, "plot_options.energy=true but LocalSolverT lacks n_energy_components()");
  }

  for (hyEdge_index_t he = 0; he < n_edges; ++he) {
    hyEdge_dofs = get_edge_dof_values<edge_dim, HyperGraphT, hyEdge_index_t, LargeVecT>(
        hyper_graph, he, lambda);

    std::array<std::array<dof_value_t, points_per_edge>, n_components> local_values;

    using bulk_fn = decltype(local_values)(
      decltype(abscissas.data())&, decltype(hyEdge_dofs)&, decltype(time));
    using bulk_fn_geom = decltype(local_values)(
      decltype(abscissas.data())&, decltype(hyEdge_dofs)&,
      decltype(hyper_graph[he])&, decltype(time));

    if constexpr (PlotFunctions::has_bulk_values<LocalSolverT, bulk_fn>::value) {
      local_values = local_solver.bulk_values(abscissas.data(), hyEdge_dofs, time);
    }
    else if constexpr (PlotFunctions::has_bulk_values<LocalSolverT, bulk_fn_geom>::value) {
      auto geometry = hyper_graph[he];
      local_values = local_solver.bulk_values(abscissas.data(), hyEdge_dofs, geometry, time);
    }
    else {
      hy_check(false, "bulk_values overload not found on LocalSolverT");
    }

    for (unsigned int p = 0; p < points_per_edge; ++p) {
      const size_t row = (static_cast<size_t>(he) * points_per_edge + p) * n_components;
      for (unsigned int d = 0; d < n_components; ++d)
        values[row + d] = static_cast<float>(local_values[d][p]);
    }

    if (plot_options.energy) {
      if constexpr (has_energy_api) {
        constexpr unsigned int n_e_comp = LocalSolverT::n_energy_components();
        auto geometry = hyper_graph[he];
        const auto local_e = local_solver.energy(hyEdge_dofs, geometry, time);
        for (unsigned int c = 0; c < cells_per_edge; ++c) {
          const size_t row = (static_cast<size_t>(he) * cells_per_edge + c) * n_e_comp;
          for (unsigned int d = 0; d < n_e_comp; ++d)
            energies_buf[row + d] = static_cast<float>(local_e[d]);
        }
      } else {
        hy_check(false, "plot_options.energy=true but LocalSolverT lacks energy()/n_energy_components()");
      }
    }
  }

  // --- distributed append: this rank concatenates its point values into the single VTKHDF part.
  // Only the data rows are appended; the per-step Steps bookkeeping is owned by rank 0 (the part
  // that created the step), so it is left untouched here.
  if (append)
  {
    hid_t file = H5Fopen(plot_options.fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    hy_check(file >= 0, "failed to open HDF5 file '" << plot_options.fileName << "'");
    hid_t root  = H5Gopen2(file, "VTKHDF", H5P_DEFAULT);
    hid_t pdata = H5Gopen2(root, "PointData", H5P_DEFAULT);
    h5_append(pdata, "values", H5T_NATIVE_FLOAT, n_points, n_components, values.data());
    H5Gclose(pdata);
    if (plot_options.energy)
    {
      if constexpr (has_energy_api)
      {
        constexpr unsigned int n_e_comp = LocalSolverT::n_energy_components();
        hid_t cdata = H5Gopen2(root, "CellData", H5P_DEFAULT);
        h5_append(cdata, "energies", H5T_NATIVE_FLOAT, n_cells, n_e_comp, energies_buf.data());
        H5Gclose(cdata);
      }
    }
    H5Gclose(root);
    H5Fclose(file);
    return;
  }

  // --- open file R/W
  hid_t file = H5Fopen(plot_options.fileName.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
  hy_check(file >= 0, "failed to open HDF5 file '" << plot_options.fileName << "'");

  hid_t root  = H5Gopen2(file, "VTKHDF",   H5P_DEFAULT);
  hid_t pdata = H5Gopen2(root, "PointData", H5P_DEFAULT);
  hid_t steps = H5Gopen2(root, "Steps",     H5P_DEFAULT);
  hid_t pdo   = H5Gopen2(steps, "PointDataOffsets", H5P_DEFAULT);

  const int64_t step_index = static_cast<int64_t>(plot_options.fileNumber);

  // --- append PointData/values: n_points new rows
  h5_append(pdata, "values", H5T_NATIVE_FLOAT,
                           n_points, n_components, values.data());

  // --- append Steps/Values: one timestamp
  double t = static_cast<double>(time);
  h5_append(steps, "Values", H5T_NATIVE_DOUBLE, 1, 1, &t);

  // --- append Steps offset entries (all zeros for single-part static mesh)
  int64_t zero = 0, one = 1;
  h5_append(steps, "PartOffsets",           H5T_NATIVE_INT64, 1, 1, &zero);
  h5_append(steps, "PointOffsets",          H5T_NATIVE_INT64, 1, 1, &zero);
  h5_append(steps, "CellOffsets",           H5T_NATIVE_INT64, 1, 1, &zero);
  h5_append(steps, "ConnectivityIdOffsets", H5T_NATIVE_INT64, 1, 1, &zero);
  h5_append(steps, "NumberOfParts",         H5T_NATIVE_INT64, 1, 1, &one);

  // --- append PointDataOffsets entries
  int64_t off_values = step_index * static_cast<int64_t>(n_points);
  h5_append(pdo, "values", H5T_NATIVE_INT64, 1, 1, &off_values);

  // --- append CellData/energies and CellDataOffsets/energies
  if (plot_options.energy) {
    if constexpr (has_energy_api) {
      constexpr unsigned int n_e_comp = LocalSolverT::n_energy_components();
      hid_t cdata = H5Gopen2(root, "CellData", H5P_DEFAULT);
      hid_t cdo   = H5Gopen2(steps, "CellDataOffsets", H5P_DEFAULT);
      hy_check(cdata >= 0 && cdo >= 0,
               "energies dataset missing: was plot_options.energy=true at fileNumber=0?");
      h5_append(cdata, "energies", H5T_NATIVE_FLOAT,
                n_cells, n_e_comp, energies_buf.data());
      int64_t off_energies = step_index * static_cast<int64_t>(n_cells);
      h5_append(cdo, "energies", H5T_NATIVE_INT64, 1, 1, &off_energies);
      H5Gclose(cdo);
      H5Gclose(cdata);
    }
  }

  // --- update NSteps
  h5_set_attr_i64(steps, "NSteps", step_index + 1);

  H5Gclose(pdo);
  H5Gclose(steps);
  H5Gclose(pdata);
  H5Gclose(root);
  H5Fclose(file);
}

// =================================================================================================
// plot_vtkhdf: dispatcher.
// =================================================================================================
template <class HyperGraphT,
          class LocalSolverT,
          typename LargeVecT,
          typename floatT,
          unsigned int n_subdivisions = 1,
          typename hyEdge_index_t = unsigned int,
          typename pt_index_t = unsigned int>
void plot_vtkhdf(HyperGraphT& hyper_graph,
                 const LocalSolverT& local_solver,
                 const LargeVecT& lambda,
                 const PlotOptions& plot_options,
                 const floatT time = 0.)
{
  // Distributed plot: each rank holds its owned hyperedges. The plot is per-edge (no shared-node
  // deduplication), so the ranks' pieces concatenate into one VTKHDF part. A token ring serializes
  // file access and carries the running point/cell/connectivity bases used to shift this rank's
  // connectivity/offsets into the concatenated part. Single rank => append=false (unchanged).
  int rank = 0, size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  hy_check(size == 1 || plot_options.fileNumber == 0,
           "distributed (multi-rank) vtkhdf plotting currently supports a single static step "
           "(fileNumber == 0) only.");

  constexpr unsigned int edge_dim        = HyperGraphT::hyEdge_dim();
  constexpr unsigned int n_subpoints     = n_subdivisions + 1;
  constexpr unsigned int points_per_edge = Hypercube<edge_dim>::pow(n_subpoints);
  constexpr unsigned int cells_per_edge  = Hypercube<edge_dim>::pow(n_subdivisions);
  constexpr unsigned int verts_per_cell  = Hypercube<edge_dim>::n_vertices();
  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  const int64_t my_np  = static_cast<int64_t>(n_edges) * points_per_edge;
  const int64_t my_nc  = static_cast<int64_t>(n_edges) * cells_per_edge;
  const int64_t my_nci = my_nc * verts_per_cell;

  const int tag = 7;
  int64_t bases[3] = {0, 0, 0};  // running totals: point_base, cell_base, conn_base
  if (rank > 0)
    MPI_Recv(bases, 3, MPI_LONG_LONG, rank - 1, tag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  const bool append = (rank > 0);

  if (plot_options.fileNumber == 0) {
    unsigned int n_e_comp = 0;
    if (plot_options.energy) {
      if constexpr (PlotFunctions::has_n_energy_components<LocalSolverT, unsigned int()>::value)
        n_e_comp = LocalSolverT::n_energy_components();
      else
        hy_check(false, "plot_options.energy=true but LocalSolverT lacks n_energy_components()");
    }
    plot_vtkhdf_mesh<HyperGraphT, n_subdivisions, hyEdge_index_t, pt_index_t>(
      hyper_graph, plot_options, LocalSolverT::system_dimension(), n_e_comp, append,
      static_cast<pt_index_t>(bases[0]), static_cast<pt_index_t>(bases[1]),
      static_cast<pt_index_t>(bases[2]));
  }
  plot_vtkhdf_bulk<HyperGraphT, LocalSolverT, LargeVecT, floatT,
                   n_subdivisions, hyEdge_index_t, pt_index_t>(
    hyper_graph, local_solver, lambda, plot_options, time, append);

  int64_t next[3] = {bases[0] + my_np, bases[1] + my_nc, bases[2] + my_nci};
  if (rank < size - 1)
    MPI_Send(next, 3, MPI_LONG_LONG, rank + 1, tag, MPI_COMM_WORLD);
}
#endif

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
  hy_check(!plot_options.fileName.empty(), "expected non-empty file name");
  hy_check(!plot_options.outputDir.empty(), "output directory must not be empty!");

  if (plot_options.fileEnding == PlotOptions::vtu)
    plot_vtu(hyper_graph, local_solver, lambda, plot_options, time);
#ifdef HYPERHDG_PETSC
  else if (plot_options.fileEnding == PlotOptions::vtkhdf)
    plot_vtkhdf(hyper_graph, local_solver, lambda, plot_options, time);
#endif
  else
    hy_assert(false, "Unsupported file ending.");

  if (plot_options.incrementFileNumber)
    ++plot_options.fileNumber;
}
