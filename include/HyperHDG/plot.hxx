#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hypercube.hxx>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <tuple>
#include <cstring>
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
   * "all", "none" or "i,j,k". Columns are written in the listed order, so readers must address
   * positions within the filtered dataset (e.g. netvis --beams-cols). Defaults to "all".
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
  else if (option == "valuesSelect")
    plot_options.values_select = value;
  else if (option == "propertiesSelect")
    plot_options.properties_select = value;
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
  else if (option == "valuesSelect")
    return_value = plot_options.values_select;
  else if (option == "propertiesSelect")
    return_value = plot_options.properties_select;
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

// --- Parallel (MPI-IO) HDF5 helpers ------------------------------------------------------------
// The vtkhdf plot is written collectively: every rank contributes its owned hyperedges as one
// contiguous slab of the single VTKHDF part, placed at a global offset obtained from an exclusive
// prefix sum (MPI_Exscan). All ranks open the file with the MPI-IO driver and all dataset
// create/extent/close calls are collective; only the data hyperslabs differ per rank. This removes
// the serial token-ring entirely and supports an arbitrary number of time steps.

// File-access property list bound to MPI-IO over MPI_COMM_WORLD (caller closes).
inline hid_t h5p_fapl_mpio()
{
  hid_t fapl = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(fapl, MPI_COMM_WORLD, MPI_INFO_NULL);
  return fapl;
}

// Collectively create a chunked dataset whose axis 0 is unlimited (extendable). gdims holds the
// GLOBAL extents; chunk0 is the axis-0 chunk length. Returns the open dataset (caller closes).
inline hid_t h5p_create(hid_t loc, const char* name, hid_t file_type,
                        int ndim, const hsize_t* gdims, hsize_t chunk0)
{
  hsize_t maxdims[2], chunk[2];
  for (int i = 0; i < ndim; ++i) maxdims[i] = gdims[i];
  maxdims[0] = H5S_UNLIMITED;
  chunk[0]   = std::max<hsize_t>(chunk0, 1);
  for (int i = 1; i < ndim; ++i) chunk[i] = std::max<hsize_t>(gdims[i], 1);
  hid_t space = H5Screate_simple(ndim, gdims, maxdims);
  hid_t dcpl  = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(dcpl, ndim, chunk);
  hid_t dset = H5Dcreate2(loc, name, file_type, space, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  H5Pclose(dcpl); H5Sclose(space);
  return dset;
}

// Collective hyperslab write of n_rows rows (n_cols columns, ignored for 1D datasets) starting at
// file row row_start. Ranks with n_rows == 0 participate with an empty selection.
inline void h5p_write_rows(hid_t dset, hid_t mem_type, hsize_t row_start,
                           hsize_t n_rows, hsize_t n_cols, const void* data)
{
  hid_t fspace = H5Dget_space(dset);
  int ndim = H5Sget_simple_extent_ndims(fspace);
  hsize_t start[2] = {row_start, 0};
  hsize_t count[2] = {n_rows, n_cols};
  hsize_t mdims[2] = {n_rows, n_cols};
  hid_t mspace = H5Screate_simple(ndim, mdims, nullptr);
  if (n_rows == 0) { H5Sselect_none(fspace); H5Sselect_none(mspace); }
  else H5Sselect_hyperslab(fspace, H5S_SELECT_SET, start, nullptr, count, nullptr);
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset, mem_type, mspace, fspace, dxpl, data);
  H5Pclose(dxpl); H5Sclose(mspace); H5Sclose(fspace);
}

// Collectively grow an open dataset's axis-0 extent to new_rows (other extents unchanged).
inline void h5p_set_rows(hid_t dset, hsize_t new_rows)
{
  hid_t fspace = H5Dget_space(dset);
  hsize_t dims[2] = {0, 0};
  H5Sget_simple_extent_dims(fspace, dims, nullptr);
  H5Sclose(fspace);
  dims[0] = new_rows;
  H5Dset_extent(dset, dims);
}

// Collectively append one scalar at index `step` to an extendable 1D dataset; all ranks perform the
// (collective) extent, only rank 0 writes the value.
inline void h5p_append1(hid_t loc, const char* name, hid_t mem_type,
                        int64_t step, const void* val, int rank)
{
  hid_t dset = H5Dopen2(loc, name, H5P_DEFAULT);
  hy_check(dset >= 0, "h5p_append1: cannot open '" << name << "'");
  hsize_t newsize = static_cast<hsize_t>(step + 1);
  H5Dset_extent(dset, &newsize);
  hid_t fspace = H5Dget_space(dset);
  hsize_t start = static_cast<hsize_t>(step), count = 1;
  hid_t mspace = H5Screate_simple(1, &count, nullptr);
  if (rank == 0) H5Sselect_hyperslab(fspace, H5S_SELECT_SET, &start, nullptr, &count, nullptr);
  else { H5Sselect_none(fspace); H5Sselect_none(mspace); }
  hid_t dxpl = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(dxpl, H5FD_MPIO_COLLECTIVE);
  H5Dwrite(dset, mem_type, mspace, fspace, dxpl, val);
  H5Pclose(dxpl); H5Sclose(mspace); H5Sclose(fspace); H5Dclose(dset);
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
                      unsigned int n_energy_components = 0)
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

  // Types: one per cell. Offsets (VTKHDF convention, length n_cells+1) are written analytically in
  // the collective section below, so no buffer is built here.
  std::vector<uint8_t> types(n_cells, element_id);

  // This rank's owned counts (concatenated into the single VTKHDF part further below).
  const int64_t np  = n_points;
  const int64_t nc  = n_cells;
  const int64_t nci = n_conn;

  // -----------------------------------------------------------------------
  // Collective layout: every rank contributes its owned edges as one slab of the single VTKHDF part.
  // Global totals (Allreduce) size the datasets; exclusive prefix sums (Exscan) give this rank's
  // starting point/cell/connectivity row. Connectivity is shifted into this rank's global point
  // range. Runtime-dependent CellData/properties presence is reduced so all ranks create the same
  // datasets even when some own zero edges.
  // -----------------------------------------------------------------------
  int mpi_rank = 0, mpi_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  long long mine[3] = {np, nc, nci};   // {points, cells, connectivity ids} owned by this rank
  long long off[3]  = {0, 0, 0};       // exclusive prefix sums (this rank's first row)
  long long tot[3]  = {0, 0, 0};       // global totals
  MPI_Exscan(mine, off, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(mine, tot, 3, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (mpi_rank == 0) { off[0] = off[1] = off[2] = 0; }

  for (auto& c : connectivity)
    c += static_cast<int64_t>(off[0]);

  int has_props_l = (n_edges > 0 && hyper_graph.hyEdge_geometry(0).has_extra_data()) ? 1 : 0;
  int n_props_l   = has_props_l
                      ? static_cast<int>(hyper_graph.hyEdge_geometry(0).extra_data().size()) : 0;
  int has_props = 0, n_properties = 0;
  MPI_Allreduce(&has_props_l, &has_props,     1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce(&n_props_l,   &n_properties,  1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

  // -----------------------------------------------------------------------
  // Write file
  // -----------------------------------------------------------------------

  // H5F_ACC_TRUNC - truncate; opened collectively via the MPI-IO driver so every rank writes its
  // slab into the shared file.
  hid_t fapl = h5p_fapl_mpio();
  hid_t file = H5Fcreate(plot_options.fileName.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, fapl);
  H5Pclose(fapl);
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

  // --- geometry: each rank writes its slab into globally-sized datasets
  { hsize_t gd[2] = {(hsize_t)tot[0], 3};
    hid_t ds = h5p_create(root, "Points", H5T_IEEE_F32LE, 2, gd, (hsize_t)tot[0]);
    h5p_write_rows(ds, H5T_NATIVE_FLOAT, (hsize_t)off[0], (hsize_t)np, 3, points.data());
    H5Dclose(ds); }
  { hsize_t gd = (hsize_t)tot[2];
    hid_t ds = h5p_create(root, "Connectivity", H5T_STD_I64LE, 1, &gd, gd);
    h5p_write_rows(ds, H5T_NATIVE_INT64, (hsize_t)off[2], (hsize_t)nci, 1, connectivity.data());
    H5Dclose(ds); }
  // Offsets has tot_cells+1 entries (VTKHDF convention, leading 0). Offsets[i] = i*verts_per_cell.
  // Rank 0 writes its cells plus the leading 0; other ranks write only their cells' running ends.
  { hsize_t gd = (hsize_t)(tot[1] + 1);
    hid_t ds = h5p_create(root, "Offsets", H5T_STD_I64LE, 1, &gd, gd);
    const hsize_t ostart = (mpi_rank == 0) ? 0 : (hsize_t)(off[1] + 1);
    const hsize_t ocount = (mpi_rank == 0) ? (hsize_t)(nc + 1) : (hsize_t)nc;
    std::vector<int64_t> obuf(ocount);
    for (hsize_t i = 0; i < ocount; ++i)
      obuf[i] = static_cast<int64_t>((ostart + i) * verts_per_cell);
    h5p_write_rows(ds, H5T_NATIVE_INT64, ostart, ocount, 1, obuf.data());
    H5Dclose(ds); }
  { hsize_t gd = (hsize_t)tot[1];
    hid_t ds = h5p_create(root, "Types", H5T_STD_U8LE, 1, &gd, gd);
    h5p_write_rows(ds, H5T_NATIVE_UINT8, (hsize_t)off[1], (hsize_t)nc, 1, types.data());
    H5Dclose(ds); }

  // --- per-piece counts (single concatenated piece → length 1, written by rank 0)
  auto write_count = [&](const char* name, int64_t total) {
    hsize_t gd = 1;
    hid_t ds = h5p_create(root, name, H5T_STD_I64LE, 1, &gd, 1);
    h5p_write_rows(ds, H5T_NATIVE_INT64, 0, (mpi_rank == 0 ? 1 : 0), 1, &total);
    H5Dclose(ds);
  };
  write_count("NumberOfPoints",          tot[0]);
  write_count("NumberOfCells",           tot[1]);
  write_count("NumberOfConnectivityIds", tot[2]);

  // --- point data group
  hid_t pdata = H5Gcreate2(root, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  // --- types points
  if constexpr (n_subdivisions == 1 && edge_dim == 1) {
    std::vector<int32_t> node_types(n_points);
    for (hyEdge_index_t he = 0; he < n_edges; ++he) {
      auto edge = hyper_graph[he];
      const pt_index_t b = he * points_per_edge;
      node_types[b+0] = static_cast<int32_t>(edge.node_descriptor[0]);
      node_types[b+1] = static_cast<int32_t>(edge.node_descriptor[1]);
    }
    hsize_t gd = (hsize_t)tot[0];
    hid_t ds = h5p_create(pdata, "types_points", H5T_STD_I32LE, 1, &gd, gd);
    h5p_write_rows(ds, H5T_NATIVE_INT32, (hsize_t)off[0], (hsize_t)n_points, 1, node_types.data());
    H5Dclose(ds);
  }

  // --- empty extendable PointData/values: (0, n_values_out); one time step occupies tot_points
  // rows. valuesSelect == "none" drops the dataset (and its Steps offsets) entirely.
  const unsigned int n_values_out =
    ComponentSelect::parse(plot_options.values_select).n_out(n_components);
  if (n_values_out > 0)
  {
    hsize_t gd[2] = {0, n_values_out};
    hid_t ds = h5p_create(pdata, "values", H5T_IEEE_F32LE, 2, gd, (hsize_t)tot[0]);
    H5Dclose(ds);
  }

  // --- static CellData/properties from hyper_edge.geometry.extra_data() (collective if any rank
  // has). propertiesSelect filters/reorders columns; "none" skips the dataset.
  hid_t cdata = H5Gcreate2(root, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  const ComponentSelect props_sel = ComponentSelect::parse(plot_options.properties_select);
  hy_check(!props_sel.mag, "'mag:' is not meaningful for propertiesSelect");
  const unsigned int n_props_out = has_props ? props_sel.n_out(n_properties) : 0;
  if (has_props && n_props_out > 0) {
    for (unsigned int c : props_sel.comps)
      hy_check((int)c < n_properties,
               "propertiesSelect column " << c << " out of range (" << n_properties << ")");
    std::vector<float> props_buf;
    hsize_t myrows = 0;
    if (n_edges > 0 && hyper_graph.hyEdge_geometry(0).has_extra_data()) {
      props_buf.assign(static_cast<size_t>(n_cells) * n_props_out, 0.f);
      myrows = (hsize_t)n_cells;
      for (hyEdge_index_t he = 0; he < n_edges; ++he) {
        auto geom = hyper_graph.hyEdge_geometry(he);
        const auto& props = geom.extra_data();
        hy_check((int)props.size() == n_properties,
                 "all hyperedges must have the same number of properties; "
                 "edge " << he << " has " << props.size() << ", expected " << n_properties);
        for (unsigned int c = 0; c < cells_per_edge; ++c) {
          const size_t row = (static_cast<size_t>(he) * cells_per_edge + c) * n_props_out;
          for (unsigned int d = 0; d < n_props_out; ++d)
            props_buf[row + d] = static_cast<float>(props[props_sel.all ? d : props_sel.comps[d]]);
        }
      }
    }
    hsize_t gd[2] = {(hsize_t)tot[1], (hsize_t)n_props_out};
    hid_t ds = h5p_create(cdata, "properties", H5T_IEEE_F32LE, 2, gd, (hsize_t)tot[1]);
    h5p_write_rows(ds, H5T_NATIVE_FLOAT, (hsize_t)off[1], myrows, (hsize_t)n_props_out,
                   props_buf.data());
    H5Dclose(ds);
  }

  // --- empty extendable CellData/energies: (0, n_energy_components); one step = tot_cells rows
  if (n_energy_components > 0) {
    hsize_t gd[2] = {0, n_energy_components};
    hid_t ds = h5p_create(cdata, "energies", H5T_IEEE_F32LE, 2, gd, (hsize_t)tot[1]);
    H5Dclose(ds);
  }
  H5Gclose(cdata);

  // --- Steps group + NSteps=0 attribute
  hid_t steps = H5Gcreate2(root, "Steps", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h5_set_attr_i64(steps, "NSteps", 0);

  auto make_empty_1d = [&](hid_t loc, const char* name, hid_t file_type) {
    hsize_t gd = 0;
    hid_t ds = h5p_create(loc, name, file_type, 1, &gd, 1);
    H5Dclose(ds);
  };

  make_empty_1d(steps, "Values",                H5T_IEEE_F64LE);
  make_empty_1d(steps, "PartOffsets",           H5T_STD_I64LE);
  make_empty_1d(steps, "PointOffsets",          H5T_STD_I64LE);
  make_empty_1d(steps, "CellOffsets",           H5T_STD_I64LE);
  make_empty_1d(steps, "ConnectivityIdOffsets", H5T_STD_I64LE);
  make_empty_1d(steps, "NumberOfParts",         H5T_STD_I64LE);

  hid_t pdo = H5Gcreate2(steps, "PointDataOffsets", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (n_values_out > 0)
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
                      const floatT time = 0.)
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

  const ComponentSelect values_sel = ComponentSelect::parse(plot_options.values_select);
  const unsigned int n_values_out = values_sel.n_out(n_components);
  for (unsigned int c : values_sel.comps)
    hy_check(c < n_components,
             "valuesSelect component " << c << " out of range (" << n_components << ")");

  std::vector<float> values(static_cast<size_t>(n_points) * n_values_out, 0.f);

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
    if (n_values_out == 0 && !plot_options.energy)
      break;  // only Steps bookkeeping left to write

    hyEdge_dofs = get_edge_dof_values<edge_dim, HyperGraphT, hyEdge_index_t, LargeVecT>(
        hyper_graph, he, lambda);

    std::array<std::array<dof_value_t, points_per_edge>, n_components> local_values;

    using bulk_fn = decltype(local_values)(
      decltype(abscissas.data())&, decltype(hyEdge_dofs)&, decltype(time));
    using bulk_fn_geom = decltype(local_values)(
      decltype(abscissas.data())&, decltype(hyEdge_dofs)&,
      decltype(hyper_graph[he])&, decltype(time));

    if (n_values_out > 0) {
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
        const size_t row = (static_cast<size_t>(he) * points_per_edge + p) * n_values_out;
        if (values_sel.all)
          for (unsigned int d = 0; d < n_components; ++d)
            values[row + d] = static_cast<float>(local_values[d][p]);
        else if (values_sel.mag) {
          dof_value_t sq = 0;
          for (unsigned int c : values_sel.comps)
            sq += local_values[c][p] * local_values[c][p];
          values[row] = static_cast<float>(std::sqrt(sq));
        }
        else
          for (unsigned int d = 0; d < n_values_out; ++d)
            values[row + d] = static_cast<float>(local_values[values_sel.comps[d]][p]);
      }
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

  // --- collective temporal append. Every rank adds its slab of this step's point values to the end
  // of the (static-mesh) PointData/values dataset; the per-step Steps bookkeeping is written once by
  // rank 0. Global totals/offsets come from the same Allreduce/Exscan as the mesh layout.
  int mpi_rank = 0, mpi_size = 1;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

  long long mine[2] = {(long long)n_points, (long long)n_cells};
  long long off[2]  = {0, 0};   // exclusive prefix (this rank's first row within the step)
  long long tot[2]  = {0, 0};   // global rows per step
  MPI_Exscan(mine, off, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(mine, tot, 2, MPI_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
  if (mpi_rank == 0) { off[0] = off[1] = 0; }

  const int64_t step = static_cast<int64_t>(plot_options.fileNumber);

  hid_t fapl = h5p_fapl_mpio();
  hid_t file = H5Fopen(plot_options.fileName.c_str(), H5F_ACC_RDWR, fapl);
  H5Pclose(fapl);
  hy_check(file >= 0, "failed to open HDF5 file '" << plot_options.fileName << "'");

  hid_t root  = H5Gopen2(file, "VTKHDF",   H5P_DEFAULT);
  hid_t pdata = H5Gopen2(root, "PointData", H5P_DEFAULT);
  hid_t steps = H5Gopen2(root, "Steps",     H5P_DEFAULT);
  hid_t pdo   = H5Gopen2(steps, "PointDataOffsets", H5P_DEFAULT);

  // --- PointData/values: grow to (step+1)*tot_points rows, write this rank's slab for this step
  if (n_values_out > 0)
  {
    hid_t ds = H5Dopen2(pdata, "values", H5P_DEFAULT);
    h5p_set_rows(ds, static_cast<hsize_t>((step + 1) * tot[0]));
    h5p_write_rows(ds, H5T_NATIVE_FLOAT, static_cast<hsize_t>(step * tot[0] + off[0]),
                   (hsize_t)n_points, n_values_out, values.data());
    H5Dclose(ds);
  }

  // --- CellData/energies (same scheme) + CellDataOffsets/energies
  if (plot_options.energy) {
    if constexpr (has_energy_api) {
      constexpr unsigned int n_e_comp = LocalSolverT::n_energy_components();
      hid_t cdata = H5Gopen2(root, "CellData", H5P_DEFAULT);
      hid_t ds = H5Dopen2(cdata, "energies", H5P_DEFAULT);
      hy_check(ds >= 0, "energies dataset missing: was plot_options.energy=true at fileNumber=0?");
      h5p_set_rows(ds, static_cast<hsize_t>((step + 1) * tot[1]));
      h5p_write_rows(ds, H5T_NATIVE_FLOAT, static_cast<hsize_t>(step * tot[1] + off[1]),
                     (hsize_t)n_cells, n_e_comp, energies_buf.data());
      H5Dclose(ds);
      hid_t cdo = H5Gopen2(steps, "CellDataOffsets", H5P_DEFAULT);
      int64_t off_energies = step * tot[1];
      h5p_append1(cdo, "energies", H5T_NATIVE_INT64, step, &off_energies, mpi_rank);
      H5Gclose(cdo);
      H5Gclose(cdata);
    }
  }

  // --- Steps bookkeeping (single concatenated part): rank 0 writes the values, all ranks extend
  double t = static_cast<double>(time);
  int64_t zero = 0, one = 1, off_values = step * tot[0];
  h5p_append1(steps, "Values",                H5T_NATIVE_DOUBLE, step, &t,    mpi_rank);
  h5p_append1(steps, "PartOffsets",           H5T_NATIVE_INT64,  step, &zero, mpi_rank);
  h5p_append1(steps, "PointOffsets",          H5T_NATIVE_INT64,  step, &zero, mpi_rank);
  h5p_append1(steps, "CellOffsets",           H5T_NATIVE_INT64,  step, &zero, mpi_rank);
  h5p_append1(steps, "ConnectivityIdOffsets", H5T_NATIVE_INT64,  step, &zero, mpi_rank);
  h5p_append1(steps, "NumberOfParts",         H5T_NATIVE_INT64,  step, &one,  mpi_rank);
  if (n_values_out > 0)
    h5p_append1(pdo,   "values",              H5T_NATIVE_INT64,  step, &off_values, mpi_rank);

  // --- update NSteps
  h5_set_attr_i64(steps, "NSteps", step + 1);

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
  // deduplication), so the ranks' slabs concatenate into one VTKHDF part. The mesh (static) is
  // written once at fileNumber == 0; each step appends a block of point/cell data. The writers open
  // the file collectively (MPI-IO) and place each rank's slab at a global offset, so any number of
  // ranks and time steps is supported with no serial token ring.
  if (plot_options.fileNumber == 0) {
    unsigned int n_e_comp = 0;
    if (plot_options.energy) {
      if constexpr (PlotFunctions::has_n_energy_components<LocalSolverT, unsigned int()>::value)
        n_e_comp = LocalSolverT::n_energy_components();
      else
        hy_check(false, "plot_options.energy=true but LocalSolverT lacks n_energy_components()");
    }
    plot_vtkhdf_mesh<HyperGraphT, n_subdivisions, hyEdge_index_t, pt_index_t>(
      hyper_graph, plot_options, LocalSolverT::system_dimension(), n_e_comp);
  }
  plot_vtkhdf_bulk<HyperGraphT, LocalSolverT, LargeVecT, floatT,
                   n_subdivisions, hyEdge_index_t, pt_index_t>(
    hyper_graph, local_solver, lambda, plot_options, time);
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
