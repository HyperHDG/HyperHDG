#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/epsilon_neighborhood_graph.hxx>
#include <HyperHDG/hy_assert.hxx>

#include <algorithm>
#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>

/*!*************************************************************************************************
 * \brief   Check whether a \c std::vector does not contain duplicate entries.
 *
 * This function returns true if the \c std::vector \c vec does not contain duplicates, and false if
 * it contains duplicates.
 *
 * \tparam  vectorT         Class template for vector.
 * \param   vec             General vector to be checked for duplicates.
 * \retval  is_unique       True if vector does not contain duplicates and false otherwise.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename vectorT>
bool is_unique(const vectorT& vec)
{
  vectorT x(vec);                                              // Copy vector enabling reorder.
  std::sort(x.begin(), x.end());                               // Sort vector in O(N log N).
  return (std::adjacent_find(x.begin(), x.end()) == x.end());  // Check for duplicates in O(n).
}

/*!*************************************************************************************************
 * \brief   Hold all topological and geometrical information of a hypergraph.
 *
 * This struct holds all topological and geometrical information of a hypergraph that has been read
 * from a file. Here, the hyperedges are supposed to be hypercuboids, i.e., we do not allow for
 * simplicial structures here. These have not yet been implemented.
 *
 * \tparam  hyEdge_dim      The local dimension of a hyperedge.
 * \tparam  space_dim       The dimension of the surrounding space.
 * \tparam  vectorT         The typename of a large vector holding e.g. all points.
 * \tparam  pointT          The typename of a point.
 * \tparam  hyEdge_index_t  The index type for hyperedges. Default is \c unsigned \c int.
 * \tparam  hyNode_index_t  The index type for hypernodes. Default is hyEdge_index_t.
 * \tparam  pt_index_t      The index type for points. Default is hyNode_index_t.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          template <typename...>
          typename vectorT,
          typename pointT,
          typename hyEdge_index_t = unsigned int,
          typename hyNode_index_t = hyEdge_index_t,
          typename pt_index_t = hyNode_index_t>
struct DomainInfo
{
  /*!***********************************************************************************************
   * \brief   Total amount of hyperedges.
   ************************************************************************************************/
  hyEdge_index_t n_hyEdges;
  /*!***********************************************************************************************
   * \brief   Total amount of hypernodes.
   ************************************************************************************************/
  hyNode_index_t n_hyNodes;
  /*!***********************************************************************************************
   * \brief   Total amount of points.
   ************************************************************************************************/
  pt_index_t n_points;

  unsigned int n_properties;
  /*!***********************************************************************************************
   * \brief   Vector containing points.
   ************************************************************************************************/
  vectorT<pointT> points;
  /*!***********************************************************************************************
   * \brief   Vector containing hypernode indices per hyperedge.
   ************************************************************************************************/
  vectorT<std::array<hyNode_index_t, 2 * hyEdge_dim> > hyNodes_hyEdge, hyFaces_hyEdge;
  /*!***********************************************************************************************
   * \brief   Vector containing vertex indices per hyperedge.
   ************************************************************************************************/
  vectorT<std::array<pt_index_t, 1 << hyEdge_dim> > points_hyEdge;  // 2 ^ hyEdge_dim
  /*!***********************************************************************************************
   * \brief   Vector containing vertex indices per hyperedge.
   ************************************************************************************************/
  vectorT<vectorT<typename pointT::value_type> > hyEdge_properties;
  /*!***********************************************************************************************
   * \brief   Constructor of DomainInfo struct.
   ************************************************************************************************/
  DomainInfo(const pt_index_t n_points,
             const hyEdge_index_t n_hyEdge,
             const hyNode_index_t n_hyNode,
             const pt_index_t n_point)
  : n_hyEdges(n_hyEdge),
    n_hyNodes(n_hyNode),
    n_points(n_point),
    n_properties(0),
    points(n_points),
    hyNodes_hyEdge(n_hyEdges),
    hyFaces_hyEdge(n_hyEdges),
    points_hyEdge(n_hyEdges)
  {
  }
  /*!***********************************************************************************************
   * \brief   Check, whether DomainInfo is in a consistent state.
   ************************************************************************************************/
  bool check_consistency()
  {
    bool consistent = true;

    std::for_each(hyNodes_hyEdge.begin(), hyNodes_hyEdge.end(),
                  [&](std::array<hyNode_index_t, 2 * hyEdge_dim> hyEdge)
                  {
                    hy_assert(is_unique(hyEdge), "Duplicate points in hyperedge.");
                    for (unsigned int i = 0; i < hyEdge.size(); ++i)
                    {
                      consistent = (hyEdge[i] < n_hyNodes && hyEdge[i] >= 0);
                      hy_assert(consistent, "At least one hypernode index is invalid!");
                    }
                  });

    std::for_each(points_hyEdge.begin(), points_hyEdge.end(),
                  [&](std::array<pt_index_t, 1 << hyEdge_dim> hyEdge)
                  {
                    hy_assert(is_unique(hyEdge), "Duplicate points in hyperedge.");
                    for (unsigned int i = 0; i < hyEdge.size(); ++i)
                    {
                      consistent = (hyEdge[i] < n_points && hyEdge[i] >= 0);
                      hy_assert(consistent, "At least one point index is invalid!");
                    }
                  });

    consistent = is_unique(points);
    hy_assert(consistent, "DomainInfo.points contains duplicate points!");
    if (!consistent)
      return false;

    consistent = is_unique(hyNodes_hyEdge);
    hy_assert(consistent, "DomainInfo.hyNodes_hyEdge contains duplicate hypernode!");
    if (!consistent)
      return false;

    consistent = is_unique(points_hyEdge);
    hy_assert(consistent, "DomainInfo.points_hyEdge contains duplicate points!");
    if (!consistent)
      return false;

    return true;
  }  // end of check_consistency
};  // end of struct DomainInfo

/*!*************************************************************************************************
 * \brief   Function to read geo file.
 *
 * Function to read self-defined .geo files are supported. These files are supposed to comprise
 * hypergraphs consisting of hypercuboids only. Other versions have not yet been implemented.
 *
 * \tparam  hyEdge_dim      The local dimension of a hyperedge.
 * \tparam  space_dim       The dimension of the surrounding space.
 * \tparam  vectorT         The typename of a large vector holding e.g. all points.
 * \tparam  pointT          The typename of a point.
 * \tparam  hyEdge_index_t  The index type for hyperedges. Default is \c unsigned \c int.
 * \tparam  hyNode_index_t  The index type for hypernodes. Default is hyEdge_index_t.
 * \tparam  pt_index_t      The index type for points. Default is hyNode_index_t.
 *
 * \param   filename        Name of the .geo file to be read.
 * \retval  domain_info     Topological and geometrical information of hypergraph.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          template <typename...> typename vectorT = std::vector,
          typename pointT = Point<space_dim, double>,
          typename hyEdge_index_t = unsigned int,
          typename hyNode_index_t = hyEdge_index_t,
          typename pt_index_t = hyNode_index_t>
DomainInfo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
read_domain_geo(const std::string& filename)
{
  hy_assert(filename.substr(filename.size() - 4, filename.size()) == ".geo",
            "The given file needs to be a .geo file for this function to be applicable!");

  std::ifstream infile(filename);
  std::istringstream linestream;
  std::string line, keyword, equal_sign;

  unsigned int Space_Dim, HyperEdge_Dim;
  pt_index_t N_Points = 0, pt_iter;
  hyEdge_index_t N_HyperEdges = 0, hyEdge_iter;
  hyNode_index_t N_HyperNodes = 0;
  pointT pt;

  while (keyword != "Space_Dim" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> Space_Dim;

  hy_assert(keyword == "Space_Dim",
            "The keyword Space_Dim has not been found in the file " << filename << "!");
  hy_assert(equal_sign == "=", "The keyword " << keyword << " has not been followd by = symbol!");
  hy_assert(Space_Dim == space_dim, "Space_Dim in " << filename << " is " << Space_Dim
                                                    << ", but should be " << space_dim << "!");

  while (keyword != "HyperEdge_Dim" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> HyperEdge_Dim;

  hy_assert(keyword == "HyperEdge_Dim",
            "The keyword HyperEdge_Dim has not been found in the file " << filename << "!");
  hy_assert(equal_sign == "=", "The keyword " << keyword << " has not been followd by = symbol!");
  hy_assert(HyperEdge_Dim == hyEdge_dim, "HyperEdge_Dim in " << filename << " is " << Space_Dim
                                                             << ", but should be " << hyEdge_dim
                                                             << "!");

  while (keyword != "N_Points" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> N_Points;

  hy_assert(keyword == "N_Points",
            "The keyword N_Points has not been found in the file " << filename << "!");
  hy_assert(equal_sign == "=", "The keyword " << keyword << " has not been followd by = symbol!");
  hy_assert(N_Points != 0, "The value of N_Points has not been set correctly!");

  while (keyword != "N_HyperNodes" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> N_HyperNodes;

  hy_assert(keyword == "N_HyperNodes",
            "The keyword N_Points has not been found in the file " << filename << "!");
  hy_assert(equal_sign == "=", "The keyword " << keyword << " has not been followd by = symbol!");
  hy_assert(N_HyperNodes != 0, "The value of N_HyperNodes has not been set correctly!");

  while (keyword != "N_HyperEdges" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> N_HyperEdges;

  hy_assert(keyword == "N_HyperEdges",
            "The keyword N_Points has not been found in the file " << filename << "!");
  hy_assert(equal_sign == "=", "The keyword " << keyword << " has not been followd by = symbol!");
  hy_assert(N_HyperEdges != 0, "The value of N_HyperEdges has not been set correctly!");

  DomainInfo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
    domain_info(N_Points, N_HyperEdges, N_HyperNodes, N_Points);

  while (keyword != "POINTS:" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }

  hy_assert(keyword == "POINTS:",
            "The keyword 'POINTS:' has not been found in the file " << filename << "!");

  for (pt_iter = 0; pt_iter < N_Points && std::getline(infile, line); ++pt_iter)
  {
    linestream = std::istringstream(line);
    for (unsigned int dim = 0; dim < space_dim; ++dim)
      linestream >> pt[dim];
    domain_info.points[pt_iter] = pt;
  }

  hy_assert(pt_iter == N_Points, "Not all points have been added to the list!");

  while (keyword != "HYPERNODES_OF_HYPEREDGES:" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }

  hy_assert(
    keyword == "HYPERNODES_OF_HYPEREDGES:",
    "The keyword 'HYPERNODES_OF_HYPEREDGES:' has not been found in the file " << filename << "!");

  for (hyEdge_iter = 0; hyEdge_iter < N_HyperEdges && std::getline(infile, line); ++hyEdge_iter)
  {
    linestream = std::istringstream(line);
    for (unsigned int i = 0; i < domain_info.hyNodes_hyEdge[hyEdge_iter].size(); ++i)
      linestream >> domain_info.hyNodes_hyEdge[hyEdge_iter][i];
  }

  hy_assert(hyEdge_iter == N_HyperEdges, "Not all hyperedges have been added to the list!");

  while (keyword != "TYPES_OF_HYPERFACES:" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }

  hy_assert(
    keyword == "TYPES_OF_HYPERFACES:",
    "The keyword 'TYPES_OF_HYPERFACES:' has not been found in the file " << filename << "!");

  for (hyEdge_iter = 0; hyEdge_iter < N_HyperEdges && std::getline(infile, line); ++hyEdge_iter)
  {
    linestream = std::istringstream(line);
    for (unsigned int i = 0; i < domain_info.hyFaces_hyEdge[hyEdge_iter].size(); ++i)
      linestream >> domain_info.hyFaces_hyEdge[hyEdge_iter][i];
  }

  hy_assert(hyEdge_iter == N_HyperEdges, "Not all hyperedges have been added to the list!");

  while (keyword != "POINTS_OF_HYPEREDGES:" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }

  hy_assert(
    keyword == "POINTS_OF_HYPEREDGES:",
    "The keyword 'POINTS_OF_HYPEREDGES:' has not been found in the file " << filename << "!");

  for (hyEdge_iter = 0; hyEdge_iter < N_HyperEdges && std::getline(infile, line); ++hyEdge_iter)
  {
    linestream = std::istringstream(line);
    for (unsigned int i = 0; i < domain_info.points_hyEdge[hyEdge_iter].size(); ++i)
      linestream >> domain_info.points_hyEdge[hyEdge_iter][i];
  }

  hy_assert(hyEdge_iter == N_HyperEdges, "Not all hyperedges have been added to the list!");

  while (keyword != "HYPEREDGE_PROPERTIES:" && std::getline(infile, line))
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }

  if (keyword != "HYPEREDGE_PROPERTIES:")
  {
    domain_info.hyEdge_properties = vectorT<vectorT<typename pointT::value_type> >();
    infile.close();
    return domain_info;
  }

  if (keyword == "HYPEREDGE_PROPERTIES:")
  {
    domain_info.hyEdge_properties = vectorT<vectorT<typename pointT::value_type> >(N_HyperEdges);
    linestream >> domain_info.n_properties;
  }

  for (hyEdge_iter = 0; hyEdge_iter < N_HyperEdges && std::getline(infile, line); ++hyEdge_iter)
  {
    linestream = std::istringstream(line);
    domain_info.hyEdge_properties[hyEdge_iter].resize(domain_info.n_properties);
    for (unsigned int i = 0; i < domain_info.hyEdge_properties[hyEdge_iter].size(); ++i)
      linestream >> domain_info.hyEdge_properties[hyEdge_iter][i];
  }

  hy_assert(hyEdge_iter == N_HyperEdges, "Not all hyperedges have been added to the list!");

  infile.close();
  return domain_info;
}  // end of read_domain_geo

/*!*************************************************************************************************
 * \brief   General Function to read domain from input file.
 *
 * At the moment, only self-defined .geo files are supported. These files are supposed to comprise
 * hypergraphs consisting of hypercuboids only. Other versions have not yet been implemented.
 *
 * \tparam  hyEdge_dim      The local dimension of a hyperedge.
 * \tparam  space_dim       The dimension of the surrounding space.
 * \tparam  vectorT         The typename of a large vector holding e.g. all points.
 * \tparam  pointT          The typename of a point.
 * \tparam  hyEdge_index_t  The index type for hyperedges. Default is \c unsigned \c int.
 * \tparam  hyNode_index_t  The index type for hypernodes. Default is hyEdge_index_t.
 * \tparam  pt_index_t      The index type for points. Default is hyNode_index_t.
 *
 * \param   filename        Name of the .geo file to be read.
 * \retval  domain_info     Topological and geometrical information of hypergraph.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          template <typename...> typename vectorT = std::vector,
          typename pointT = Point<space_dim, double>,
          typename hyEdge_index_t = unsigned int,
          typename hyNode_index_t = hyEdge_index_t,
          typename pt_index_t = hyNode_index_t>
DomainInfo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
read_domain(std::string filename)
{
  hy_assert(std::filesystem::exists(filename), "File does not exist.");

  if (filename.substr(filename.size() - 4, filename.size()) == ".pts")
  {
    hy_assert(hyEdge_dim == 1, "This only works for graphs, so far!");
    make_epsilon_neighborhood_graph<space_dim, vectorT, pointT, hyEdge_index_t>(filename);
  }

  hy_assert(filename.substr(filename.size() - 4, filename.size()) == ".geo",
            "The given file needs to be a .geo file, since no other input file types are currently"
              << " implemented.");

  DomainInfo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
    domain_info = read_domain_geo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t,
                                  hyNode_index_t, pt_index_t>(filename);

  hy_assert(domain_info.check_consistency(),
            "Domain info appears to be inconsistent!"
              << std::endl
              << "This assertion is never to be thrown since it can only be caused by internal "
              << "assertions of DomainInfo.check_consistency()!");

  return domain_info;
}  // end of read_domain
