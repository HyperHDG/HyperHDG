#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hypercube.hxx>
#ifdef HYPERHDG_HAVE_HDF5
#include <HyperHDG/read_domain_hdf5.hxx>
#endif

#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <vector>

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
    vtkhdf
  };
  /*!***********************************************************************************************
   * \brief   File ending and also way of plotting.
   *
   * "vtu" writes one XML file per call. "vtkhdf" (needs HyperHDG configured with HDF5) writes
   * one HDF5 file in the VTKHDF format: the mesh once, and every call with \c fileNumber > 0
   * appends a time step, so a transient run becomes one file that ParaView reads as a time
   * series (cf. plot_vtkhdf). Default is "vtu".
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
   * Requires the local solver to expose \c n_energy_components() and \c energy(). Defaults to
   * false.
   ************************************************************************************************/
  bool energy = false;
  /*!***********************************************************************************************
   * \brief   Component selection for PointData/values (vtkhdf only).
   *
   * "all", "none", "i,j,k" (keep the listed components, in that order) or "mag:i,j,k" (a single
   * column holding the euclidean norm of the listed components). The per-step values block is the
   * dominant payload for large networks, so restricting it (e.g. to the displacement components)
   * scales the file size by n_selected / n_components. Defaults to "all".
   ************************************************************************************************/
  std::string values_select = "all";
  /*!***********************************************************************************************
   * \brief   Column selection for the static CellData/properties (vtkhdf only).
   *
   * "all", "none" or "i,j,k". Columns are written in the listed order. Defaults to "all".
   ************************************************************************************************/
  std::string properties_select = "all";
};  // end of class PlotOptions

/*!*************************************************************************************************
 * \brief   Parsed component selection ("all" | "none" | "i,j,k" | "mag:i,j,k"), see PlotOptions.
 **************************************************************************************************/
struct ComponentSelect
{
  bool all = true;                  ///< keep every component
  bool mag = false;                 ///< single euclidean-norm column over #comps
  std::vector<unsigned int> comps;  ///< selected source components (empty for "all"/"none")

  static ComponentSelect parse(const std::string& spec)
  {
    ComponentSelect sel;
    if (spec.empty() || spec == "all")
      return sel;
    sel.all = false;
    if (spec == "none")
      return sel;
    std::string list = spec;
    if (spec.rfind("mag:", 0) == 0)
    {
      sel.mag = true;
      list = spec.substr(4);
    }
    size_t pos = 0;
    while (pos < list.size())
    {
      size_t next = list.find(',', pos);
      if (next == std::string::npos)
        next = list.size();
      sel.comps.push_back(std::stoul(list.substr(pos, next - pos)));
      pos = next + 1;
    }
    hy_check(!sel.comps.empty(), "empty component list in select spec '" << spec << "'");
    return sel;
  }

  unsigned int n_out(unsigned int n_full) const
  {
    if (all)
      return n_full;
    if (mag)
      return 1;
    return comps.size();
  }
};  // end of struct ComponentSelect

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
    else if (value == "vtkhdf")
      plot_options.fileEnding = PlotOptions::vtkhdf;
    else
      hy_assert(false, "You have chosen an invalid file type!");
  }
  else if (option == "energy")
    plot_options.energy = (value == "true" || value == "1");
  else if (option == "valuesSelect")
    plot_options.values_select = value;
  else if (option == "propertiesSelect")
    plot_options.properties_select = value;
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
    return_value = plot_options.fileEnding == PlotOptions::vtu ? "vtu" : "vtkhdf";
  else if (option == "energy")
    return_value = std::to_string(plot_options.energy);
  else if (option == "valuesSelect")
    return_value = plot_options.values_select;
  else if (option == "propertiesSelect")
    return_value = plot_options.properties_select;
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
    if constexpr (requires {
                    local_solver.bulk_values(abscissas.data(), hyEdge_dofs, prvalue_of(time));
                  })
      local_values = local_solver.bulk_values(abscissas.data(), hyEdge_dofs, time);
    else if constexpr (requires(decltype(hyper_graph[he_number])& geometry) {
                         local_solver.bulk_values(abscissas.data(), hyEdge_dofs, geometry,
                                                  prvalue_of(time));
                       })
    {
      auto geometry = hyper_graph[he_number];
      local_values = local_solver.bulk_values(abscissas.data(), hyEdge_dofs, geometry, time);
    }
    else
      static_assert(always_false_v<LocalSolverT>,
                    "LocalSolverT implements no usable overload of bulk_values!");

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

#ifdef HYPERHDG_HAVE_HDF5
/*!*************************************************************************************************
 * \brief   The VTKHDF writer (PlotOptions::vtkhdf): one HDF5 file per run, read by ParaView.
 *
 * Layout (VTKHDF 2.0, UnstructuredGrid, one part, transient):
 *
 *     /VTKHDF/Points, Connectivity, Offsets, Types, NumberOf*   the mesh, written at fileNumber 0
 *     /VTKHDF/PointData/values   (n_steps * n_points, n_values)  float32, one block per time step
 *     /VTKHDF/PointData/types_points                            node types (networks)
 *     /VTKHDF/CellData/properties (n_cells, n_properties)       geometry's extra data, if any
 *     /VTKHDF/CellData/energies  (n_steps * n_cells, n_energies) if PlotOptions::energy
 *     /VTKHDF/Steps/...                                         time values and block offsets
 *
 * Every hyperedge is plotted on its own (no sharing of points between hyperedges, like the vtu
 * writer with plot_edges). The mesh is written once, when \c fileNumber is 0; every call appends
 * one time step (\c fileNumber is the step, \c time its value), so a transient run gives one file
 * that ParaView reads as a time series. The datasets are chunked with an unlimited first axis;
 * the same layout, written collectively with one slab per rank, is what a distributed run will
 * produce later.
 **************************************************************************************************/
namespace PlotVTKHDF
{
using HDF5Domain::Id;

/*!*************************************************************************************************
 * \brief   The file name of the plot: outputDir/fileName.vtkhdf (the directory is created).
 **************************************************************************************************/
inline std::string file_path(const PlotOptions& plot_options)
{
  if (std::filesystem::create_directory(plot_options.outputDir))
    std::cout << "Directory \"" << plot_options.outputDir << "\" has been created." << std::endl;
  return plot_options.outputDir + "/" + plot_options.fileName + ".vtkhdf";
}

/*!*************************************************************************************************
 * \brief   Create a chunked dataset whose first axis is unlimited (initial extents \c dims).
 **************************************************************************************************/
inline hid_t create_extendable(const hid_t loc,
                               const char* name,
                               const hid_t file_type,
                               const int ndim,
                               const hsize_t* dims,
                               const hsize_t chunk_rows)
{
  hsize_t maxdims[2] = {H5S_UNLIMITED, ndim == 2 ? dims[1] : 0};
  hsize_t chunk[2] = {std::max<hsize_t>(chunk_rows, 1),
                      ndim == 2 ? std::max<hsize_t>(dims[1], 1) : 0};
  Id space(H5Screate_simple(ndim, dims, maxdims), H5Sclose, "create dataspace");
  Id dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, "create dataset property list");
  H5Pset_chunk(dcpl, ndim, chunk);
  const hid_t dset = H5Dcreate2(loc, name, file_type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  hy_check(dset >= 0, "HDF5: creating dataset " << name << " failed.");
  return dset;
}

/*!*************************************************************************************************
 * \brief   Write \c n_rows rows (of \c n_cols columns for 2d datasets) starting at \c row_start.
 **************************************************************************************************/
inline void write_rows(const hid_t dset,
                       const hid_t mem_type,
                       const hsize_t row_start,
                       const hsize_t n_rows,
                       const hsize_t n_cols,
                       const void* data)
{
  if (n_rows == 0)
    return;
  Id fspace(H5Dget_space(dset), H5Sclose, "get dataspace");
  const int ndim = H5Sget_simple_extent_ndims(fspace);
  const hsize_t start[2] = {row_start, 0}, count[2] = {n_rows, n_cols};
  Id mspace(H5Screate_simple(ndim, count, nullptr), H5Sclose, "create memory dataspace");
  H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, nullptr, count, nullptr);
  hy_check(H5Dwrite(dset, mem_type, mspace, fspace, H5P_DEFAULT, data) >= 0,
           "HDF5: writing rows failed.");
}

/*!*************************************************************************************************
 * \brief   Grow the first axis of an open dataset to \c new_rows.
 **************************************************************************************************/
inline void set_rows(const hid_t dset, const hsize_t new_rows)
{
  Id fspace(H5Dget_space(dset), H5Sclose, "get dataspace");
  hsize_t dims[2] = {0, 0};
  H5Sget_simple_extent_dims(fspace, dims, nullptr);
  dims[0] = new_rows;
  hy_check(H5Dset_extent(dset, dims) >= 0, "HDF5: extending dataset failed.");
}

/*!*************************************************************************************************
 * \brief   Store one scalar at index \c step of an extendable 1d dataset (grown to step + 1).
 **************************************************************************************************/
inline void append_scalar(const hid_t loc,
                          const char* name,
                          const hid_t mem_type,
                          const hsize_t step,
                          const void* value)
{
  Id dset(H5Dopen2(loc, name, H5P_DEFAULT), H5Dclose, std::string("open dataset ") + name);
  set_rows(dset, step + 1);
  write_rows(dset, mem_type, step, 1, 1, value);
}

/*!*************************************************************************************************
 * \brief   Set an integer attribute (attributes cannot be resized, so delete and recreate).
 **************************************************************************************************/
inline void set_attribute(const hid_t loc, const char* name, const long long value)
{
  if (H5Aexists(loc, name) > 0)
    H5Adelete(loc, name);
  HDF5Domain::write_attribute(loc, name, value);
}

/*!*************************************************************************************************
 * \brief   Write the mesh and the empty, extendable data arrays (truncates the file).
 **************************************************************************************************/
template <class HyperGraphT,
          unsigned int n_subdivisions = 1,
          typename hyEdge_index_t = unsigned int,
          typename pt_index_t = unsigned int>
void write_mesh(HyperGraphT& hyper_graph,
                const PlotOptions& plot_options,
                const unsigned int n_components,
                const unsigned int n_energy_components)
{
  constexpr unsigned int edge_dim = HyperGraphT::hyEdge_dim();
  constexpr unsigned int space_dim = HyperGraphT::space_dim();
  static_assert(edge_dim >= 1 && edge_dim <= 3 && space_dim <= 3,
                "VTKHDF plots exist for lines, quads and hexahedra in up to 3 dimensions.");
  constexpr unsigned int n_subpoints = n_subdivisions + 1;
  constexpr unsigned int points_per_edge = Hypercube<edge_dim>::pow(n_subpoints);
  constexpr unsigned int cells_per_edge = Hypercube<edge_dim>::pow(n_subdivisions);
  constexpr unsigned int verts_per_cell = Hypercube<edge_dim>::n_vertices();
  constexpr unsigned char vtk_type = edge_dim == 1   ? 3
                                     : edge_dim == 2 ? 8
                                                     : 11;  // line/pixel/voxel

  SmallVec<n_subpoints, float> abscissas;
  for (unsigned int i = 0; i < n_subpoints; ++i)
    abscissas[i] = plot_options.scale * (1.f * i / n_subdivisions - 0.5f) + 0.5f;

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  const hsize_t n_points = (hsize_t)n_edges * points_per_edge;
  const hsize_t n_cells = (hsize_t)n_edges * cells_per_edge;
  const hsize_t n_conn = n_cells * verts_per_cell;

  std::vector<float> points(3 * n_points, 0.f);
  for (hyEdge_index_t he = 0; he < n_edges; ++he)
  {
    auto edge = hyper_graph.hyEdge_geometry(he);
    for (unsigned int p = 0; p < points_per_edge; ++p)
    {
      const Point<space_dim> pt =
        (Point<space_dim>)edge.template lexicographic<n_subpoints>(p, abscissas);
      for (unsigned int d = 0; d < space_dim; ++d)
        points[3 * (he * points_per_edge + p) + d] = pt[d];
    }
  }

  // sub-cubes of every hyperedge in lexicographic order; their vertices in VTK's pixel/voxel
  // order, which is the lexicographic one
  std::vector<long long> connectivity;
  connectivity.reserve(n_conn);
  for (hyEdge_index_t he = 0; he < n_edges; ++he)
    for (unsigned int cell = 0; cell < cells_per_edge; ++cell)
      for (unsigned int vertex = 0; vertex < verts_per_cell; ++vertex)
      {
        unsigned int index = 0, stride = 1, rest = cell;
        for (unsigned int k = 0; k < edge_dim; ++k, rest /= n_subdivisions, stride *= n_subpoints)
          index += (rest % n_subdivisions + ((vertex >> k) & 1)) * stride;
        connectivity.push_back((long long)he * points_per_edge + index);
      }
  std::vector<long long> offsets(n_cells + 1);
  for (hsize_t cell = 0; cell <= n_cells; ++cell)
    offsets[cell] = cell * verts_per_cell;
  const std::vector<unsigned char> types(n_cells, vtk_type);

  Id file(H5Fcreate(file_path(plot_options).c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT),
          H5Fclose, "create file " + file_path(plot_options));
  Id root(H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
          "create group VTKHDF");
  {
    const long long version[2] = {2, 0};
    const hsize_t two = 2;
    Id space(H5Screate_simple(1, &two, nullptr), H5Sclose, "create dataspace");
    Id attr(H5Acreate2(root, "Version", H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose,
            "create attribute Version");
    H5Awrite(attr, H5T_NATIVE_LLONG, version);
  }
  HDF5Domain::write_attribute(root, "Type", std::string("UnstructuredGrid"));

  {
    const hsize_t dims[2] = {n_points, 3};
    Id ds(create_extendable(root, "Points", H5T_IEEE_F32LE, 2, dims, n_points), H5Dclose, "");
    write_rows(ds, H5T_NATIVE_FLOAT, 0, n_points, 3, points.data());
  }
  {
    Id ds(create_extendable(root, "Connectivity", H5T_STD_I64LE, 1, &n_conn, n_conn), H5Dclose, "");
    write_rows(ds, H5T_NATIVE_LLONG, 0, n_conn, 1, connectivity.data());
  }
  {
    const hsize_t n = n_cells + 1;
    Id ds(create_extendable(root, "Offsets", H5T_STD_I64LE, 1, &n, n), H5Dclose, "");
    write_rows(ds, H5T_NATIVE_LLONG, 0, n, 1, offsets.data());
  }
  {
    Id ds(create_extendable(root, "Types", H5T_STD_U8LE, 1, &n_cells, n_cells), H5Dclose, "");
    write_rows(ds, H5T_NATIVE_UCHAR, 0, n_cells, 1, types.data());
  }
  const long long counts[3] = {(long long)n_points, (long long)n_cells, (long long)n_conn};
  const char* count_names[3] = {"NumberOfPoints", "NumberOfCells", "NumberOfConnectivityIds"};
  for (unsigned int k = 0; k < 3; ++k)
  {
    const hsize_t one = 1;
    Id ds(create_extendable(root, count_names[k], H5T_STD_I64LE, 1, &one, 1), H5Dclose, "");
    write_rows(ds, H5T_NATIVE_LLONG, 0, 1, 1, counts + k);
  }

  Id pdata(H5Gcreate2(root, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
           "create group PointData");
  if constexpr (n_subdivisions == 1 && edge_dim == 1)
  {
    std::vector<int> node_types(n_points);
    for (hyEdge_index_t he = 0; he < n_edges; ++he)
    {
      auto edge = hyper_graph[he];
      node_types[2 * he] = edge.node_descriptor[0];
      node_types[2 * he + 1] = edge.node_descriptor[1];
    }
    Id ds(create_extendable(pdata, "types_points", H5T_STD_I32LE, 1, &n_points, n_points), H5Dclose,
          "");
    write_rows(ds, H5T_NATIVE_INT, 0, n_points, 1, node_types.data());
  }
  const unsigned int n_values_out =
    ComponentSelect::parse(plot_options.values_select).n_out(n_components);
  if (n_values_out > 0)
  {
    const hsize_t dims[2] = {0, n_values_out};
    Id(create_extendable(pdata, "values", H5T_IEEE_F32LE, 2, dims, n_points), H5Dclose, "");
  }

  Id cdata(H5Gcreate2(root, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
           "create group CellData");
  if constexpr (requires(decltype(hyper_graph.hyEdge_geometry(0)) g) { g.has_extra_data(); })
  {
    const ComponentSelect props_sel = ComponentSelect::parse(plot_options.properties_select);
    hy_check(!props_sel.mag, "'mag:' is not meaningful for propertiesSelect");
    auto first = hyper_graph.hyEdge_geometry(0);
    const unsigned int n_properties =
      n_edges > 0 && first.has_extra_data() ? first.extra_data().size() : 0;
    const unsigned int n_props_out = props_sel.n_out(n_properties);
    if (n_properties > 0 && n_props_out > 0)
    {
      for (unsigned int c : props_sel.comps)
        hy_check(c < n_properties,
                 "propertiesSelect column " << c << " out of range (" << n_properties << ")");
      std::vector<float> props(n_cells * n_props_out);
      for (hyEdge_index_t he = 0; he < n_edges; ++he)
      {
        auto geometry = hyper_graph.hyEdge_geometry(he);
        const auto& data = geometry.extra_data();
        hy_check(data.size() == n_properties,
                 "all hyperedges must have the same number of properties; edge "
                   << he << " has " << data.size() << ", expected " << n_properties);
        for (unsigned int cell = 0; cell < cells_per_edge; ++cell)
          for (unsigned int d = 0; d < n_props_out; ++d)
            props[((hsize_t)he * cells_per_edge + cell) * n_props_out + d] =
              data[props_sel.all ? d : props_sel.comps[d]];
      }
      const hsize_t dims[2] = {n_cells, n_props_out};
      Id ds(create_extendable(cdata, "properties", H5T_IEEE_F32LE, 2, dims, n_cells), H5Dclose, "");
      write_rows(ds, H5T_NATIVE_FLOAT, 0, n_cells, n_props_out, props.data());
    }
  }
  if (n_energy_components > 0)
  {
    const hsize_t dims[2] = {0, n_energy_components};
    Id(create_extendable(cdata, "energies", H5T_IEEE_F32LE, 2, dims, n_cells), H5Dclose, "");
  }

  Id steps(H5Gcreate2(root, "Steps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
           "create group Steps");
  set_attribute(steps, "NSteps", 0);
  auto empty_1d = [](const hid_t loc, const char* name, const hid_t type)
  {
    const hsize_t zero = 0;
    Id(create_extendable(loc, name, type, 1, &zero, 1), H5Dclose, "");
  };
  empty_1d(steps, "Values", H5T_IEEE_F64LE);
  for (const char* name :
       {"PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets", "NumberOfParts"})
    empty_1d(steps, name, H5T_STD_I64LE);
  Id pdo(H5Gcreate2(steps, "PointDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
         "create group PointDataOffsets");
  if (n_values_out > 0)
    empty_1d(pdo, "values", H5T_STD_I64LE);
  if (n_energy_components > 0)
  {
    Id cdo(H5Gcreate2(steps, "CellDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
           "create group CellDataOffsets");
    empty_1d(cdo, "energies", H5T_STD_I64LE);
  }
}

/*!*************************************************************************************************
 * \brief   Append one time step: PointData/values (and CellData/energies) and the Steps entries.
 **************************************************************************************************/
template <class HyperGraphT,
          class LocalSolverT,
          typename LargeVecT,
          typename floatT,
          unsigned int n_subdivisions = 1,
          typename hyEdge_index_t = unsigned int>
void write_step(HyperGraphT& hyper_graph,
                const LocalSolverT& local_solver,
                const LargeVecT& lambda,
                const PlotOptions& plot_options,
                const floatT time)
{
  constexpr unsigned int edge_dim = HyperGraphT::hyEdge_dim();
  constexpr unsigned int n_subpoints = n_subdivisions + 1;
  constexpr unsigned int points_per_edge = Hypercube<edge_dim>::pow(n_subpoints);
  constexpr unsigned int cells_per_edge = Hypercube<edge_dim>::pow(n_subdivisions);
  constexpr unsigned int n_components = LocalSolverT::system_dimension();
  constexpr bool has_energy_api = requires { LocalSolverT::n_energy_components(); };
  using dof_value_t = typename LargeVecT::value_type;
  using plot_value_t = decltype(std::real(std::declval<dof_value_t>()));

  SmallVec<n_subpoints, float> abscissas;
  for (unsigned int i = 0; i < n_subpoints; ++i)
    abscissas[i] = plot_options.scale * (1.f * i / n_subdivisions - 0.5f) + 0.5f;

  const hyEdge_index_t n_edges = hyper_graph.n_hyEdges();
  const hsize_t n_points = (hsize_t)n_edges * points_per_edge;
  const hsize_t n_cells = (hsize_t)n_edges * cells_per_edge;
  const hsize_t step = plot_options.fileNumber;

  const ComponentSelect values_sel = ComponentSelect::parse(plot_options.values_select);
  const unsigned int n_values_out = values_sel.n_out(n_components);
  for (unsigned int c : values_sel.comps)
    hy_check(c < n_components,
             "valuesSelect component " << c << " out of range (" << n_components << ")");
  if (plot_options.energy)
    hy_check(has_energy_api,
             "energy = true, but the local solver has no n_energy_components() / energy().");

  std::vector<float> values(n_points * n_values_out), energies;
  if constexpr (has_energy_api)
    if (plot_options.energy)
      energies.resize(n_cells * LocalSolverT::n_energy_components());

  if (n_values_out > 0 || plot_options.energy)
    for (hyEdge_index_t he = 0; he < n_edges; ++he)
    {
      const auto hyEdge_dofs =
        get_edge_dof_values<edge_dim, HyperGraphT, hyEdge_index_t, LargeVecT>(hyper_graph, he,
                                                                              lambda);
      if (n_values_out > 0)
      {
        std::array<std::array<plot_value_t, points_per_edge>, n_components> local_values;
        if constexpr (requires {
                        local_solver.bulk_values(abscissas.data(), hyEdge_dofs, prvalue_of(time));
                      })
          local_values = local_solver.bulk_values(abscissas.data(), hyEdge_dofs, time);
        else if constexpr (requires(decltype(hyper_graph[he])& geometry) {
                             local_solver.bulk_values(abscissas.data(), hyEdge_dofs, geometry,
                                                      time);
                           })
        {
          auto geometry = hyper_graph[he];
          local_values = local_solver.bulk_values(abscissas.data(), hyEdge_dofs, geometry, time);
        }
        else
          static_assert(always_false_v<LocalSolverT>,
                        "LocalSolverT implements no usable overload of bulk_values!");

        for (unsigned int p = 0; p < points_per_edge; ++p)
        {
          float* row = values.data() + ((hsize_t)he * points_per_edge + p) * n_values_out;
          if (values_sel.all)
            for (unsigned int d = 0; d < n_components; ++d)
              row[d] = local_values[d][p];
          else if (values_sel.mag)
          {
            plot_value_t sq = 0;
            for (unsigned int c : values_sel.comps)
              sq += local_values[c][p] * local_values[c][p];
            row[0] = std::sqrt(sq);
          }
          else
            for (unsigned int d = 0; d < n_values_out; ++d)
              row[d] = local_values[values_sel.comps[d]][p];
        }
      }
      if constexpr (has_energy_api)
        if (plot_options.energy)
        {
          constexpr unsigned int n_e = LocalSolverT::n_energy_components();
          auto geometry = hyper_graph[he];
          const auto local_e = local_solver.energy(hyEdge_dofs, geometry, time);
          for (unsigned int cell = 0; cell < cells_per_edge; ++cell)
            for (unsigned int d = 0; d < n_e; ++d)
              energies[((hsize_t)he * cells_per_edge + cell) * n_e + d] = local_e[d];
        }
    }

  Id file(H5Fopen(file_path(plot_options).c_str(), H5F_ACC_RDWR, H5P_DEFAULT), H5Fclose,
          "open file " + file_path(plot_options));
  Id root(H5Gopen2(file, "VTKHDF", H5P_DEFAULT), H5Gclose, "open group VTKHDF");
  Id steps(H5Gopen2(root, "Steps", H5P_DEFAULT), H5Gclose, "open group Steps");

  if (n_values_out > 0)
  {
    Id pdata(H5Gopen2(root, "PointData", H5P_DEFAULT), H5Gclose, "open group PointData");
    Id ds(H5Dopen2(pdata, "values", H5P_DEFAULT), H5Dclose, "open dataset values");
    set_rows(ds, (step + 1) * n_points);
    write_rows(ds, H5T_NATIVE_FLOAT, step * n_points, n_points, n_values_out, values.data());
    Id pdo(H5Gopen2(steps, "PointDataOffsets", H5P_DEFAULT), H5Gclose,
           "open group PointDataOffsets");
    const long long offset = step * n_points;
    append_scalar(pdo, "values", H5T_NATIVE_LLONG, step, &offset);
  }
  if constexpr (has_energy_api)
    if (plot_options.energy)
    {
      Id cdata(H5Gopen2(root, "CellData", H5P_DEFAULT), H5Gclose, "open group CellData");
      Id ds(H5Dopen2(cdata, "energies", H5P_DEFAULT), H5Dclose,
            "open dataset energies (was energy = true at fileNumber 0?)");
      set_rows(ds, (step + 1) * n_cells);
      write_rows(ds, H5T_NATIVE_FLOAT, step * n_cells, n_cells, LocalSolverT::n_energy_components(),
                 energies.data());
      Id cdo(H5Gopen2(steps, "CellDataOffsets", H5P_DEFAULT), H5Gclose,
             "open group CellDataOffsets");
      const long long offset = step * n_cells;
      append_scalar(cdo, "energies", H5T_NATIVE_LLONG, step, &offset);
    }

  const double t = time;
  const long long zero = 0, one = 1;
  append_scalar(steps, "Values", H5T_NATIVE_DOUBLE, step, &t);
  for (const char* name : {"PartOffsets", "PointOffsets", "CellOffsets", "ConnectivityIdOffsets"})
    append_scalar(steps, name, H5T_NATIVE_LLONG, step, &zero);
  append_scalar(steps, "NumberOfParts", H5T_NATIVE_LLONG, step, &one);
  set_attribute(steps, "NSteps", step + 1);
}

}  // end of namespace PlotVTKHDF

/*!*************************************************************************************************
 * \brief   Plot in the VTKHDF format: the mesh at fileNumber 0, then one time step per call.
 **************************************************************************************************/
template <class HyperGraphT, class LocalSolverT, typename LargeVecT, typename floatT>
void plot_vtkhdf(HyperGraphT& hyper_graph,
                 const LocalSolverT& local_solver,
                 const LargeVecT& lambda,
                 const PlotOptions& plot_options,
                 const floatT time)
{
  if (plot_options.fileNumber == 0)
  {
    unsigned int n_energy_components = 0;
    if constexpr (requires { LocalSolverT::n_energy_components(); })
      if (plot_options.energy)
        n_energy_components = LocalSolverT::n_energy_components();
    PlotVTKHDF::write_mesh(hyper_graph, plot_options, LocalSolverT::system_dimension(),
                           n_energy_components);
  }
  PlotVTKHDF::write_step(hyper_graph, local_solver, lambda, plot_options, time);
}
#endif  // HYPERHDG_HAVE_HDF5

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
  std::cout << plot_options.fileName;
  if (plot_options.printFileNumber)
    std::cout << "." << plot_options.fileNumber << " was written\n";
  PlotFunctions::close_ofstream(myfile);
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
  hy_assert(!plot_options.fileName.empty(), "File name must not be empty!");
  hy_assert(!plot_options.outputDir.empty(), "Ouput directory must not be empty!");
  if (plot_options.fileEnding == PlotOptions::vtu)
    plot_vtu(hyper_graph, local_solver, lambda, plot_options, time);
  else
  {
#ifdef HYPERHDG_HAVE_HDF5
    plot_vtkhdf(hyper_graph, local_solver, lambda, plot_options, time);
#else
    hy_check(false, "Plotting as vtkhdf needs HyperHDG configured with HDF5 (HYPERHDG_HAVE_HDF5): "
                      << "either CMake did not find HDF5, or this is a Cython module.");
#endif
  }
  if (plot_options.incrementFileNumber)
    ++plot_options.fileNumber;
}  // end of void plot
