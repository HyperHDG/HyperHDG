#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/read_domain.hxx>

#include <string>

namespace Topology
{
/*!*************************************************************************************************
 * \brief   Hypergraph topology based on an input file.
 *
 * The topology class File is a set of hyperedges. Each of these tensorial hyperedges is represented
 * by its hypernodes (given within the file). For consistency, it is assumed that the vertices and
 * the hypernodes are assumed to be given in lexicographical order to ensure that geometry and
 * topology of all hyperedges fit.
 *
 * \tparam  hyEdge_dimT     Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT      The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyEdge_dimT.
 * \tparam  hyEdge_index_t  The index type for hyperedges. Default is \c unsigned \c int.
 * \tparam  hyNode_index_t  The index type for hypernodes. Default is \c hyEdge_index_t.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dimT,
          template <typename> typename vectorT = std::vector,
          typename pointT = Point<space_dimT, float>,
          typename hyEdge_index_t = unsigned int,
          typename hyNode_index_t = hyEdge_index_t,
          typename pt_index_t = hyNode_index_t>
class File
{
  /*!***********************************************************************************************
   * \brief   Definition of the topology of a hypergraph's edges.
   ************************************************************************************************/
  class hyEdge
  {
   public:
    static constexpr unsigned int n_hyNodes() { return 2 * hyEdge_dimT; }

   private:
    /*!*********************************************************************************************
     * \brief   Reference to parent hypergraph.
     **********************************************************************************************/
    const File& hyGraph_topology_;
    /*!*********************************************************************************************
     * \brief   Index of the hyperedge within the hypergraph
     **********************************************************************************************/
    const hyEdge_index_t index_;

   public:
    /*!*********************************************************************************************
     * \brief   Construct hyperedge from hypergraph and index.
     **********************************************************************************************/
    hyEdge(const File& hyGraph_topology, const hyEdge_index_t index)
    : hyGraph_topology_(hyGraph_topology), index_(index)
    {
    }
    /*!*********************************************************************************************
     * \brief   Return hypernodes of a hyperedge.
     **********************************************************************************************/
    const auto& get_hyNode_indices() const
    {
      return hyGraph_topology_.domain_info_.hyNodes_hyEdge[index_];
    }
  };  // end of class hyEdge

 public:
  /*!***********************************************************************************************
   * \brief   Return local dimension of the hypergraph's hyperedges.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Return dimension of the surrounding space.
   ************************************************************************************************/
  static constexpr unsigned int space_dim() { return space_dimT; }

 private:
  /*!***********************************************************************************************
   * \brief   Domain Info containing all the information of the hypergraph (cf. ReadDomain.hxx).
   ************************************************************************************************/
  const DomainInfo<hyEdge_dimT,
                   space_dimT,
                   vectorT,
                   pointT,
                   hyEdge_index_t,
                   hyNode_index_t,
                   pt_index_t>
    domain_info_;

 public:
  /*!***********************************************************************************************
   * \brief   Defines the return value of the class.
   *
   * The \c class \c HyperGraph_Cubic defines the topology of the hypergraph. It "contains" the
   * different hyperedges (that actually are constructed everytime access is needed from e.g. the
   * solver class). Thus, its main purpose is to provide a structure that administrates the
   * hyperedges that are the return value of this structure.
   ************************************************************************************************/
  typedef hyEdge value_type;
  /*!***********************************************************************************************
   * \brief   Defines the value type of input argument for standard constructor.
   *
   * To receive a very general \c AbstractProblem, constructors need to account for the fact that
   * the specific topology / geometry of a hypergraph influences the way in which the hypergraph
   * needs to be constructed. The \c typedef implements the aspect, that a cubic hypergraph
   * topology is by default constructed by a std::vector that contains amounts of elements in the
   * different dimensions.
   ************************************************************************************************/
  typedef std::string constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Construct a topology from a given filename.
   *
   * \param   filename    Name of file containing the information.
   ************************************************************************************************/
  File(const constructor_value_type& filename)
  : domain_info_(read_domain<hyEdge_dimT,
                             space_dimT,
                             vectorT,
                             pointT,
                             hyEdge_index_t,
                             hyNode_index_t,
                             pt_index_t>(filename))
  {
  }
  /*!***********************************************************************************************
   * \brief   Copy constructor.
   ************************************************************************************************/
  File(const File<hyEdge_dimT, space_dimT>& other) : domain_info_(other.domain_info) {}

  /*!***********************************************************************************************
   * \brief   Get topological hyperedge of given index.
   *
   * \todo  Here, we repeatedly return a large object. This is done since the object could be
   *        locally created in regular topologies/geometries! Return shared-pointer?
   *        -> This is not really a large object, is it? I mean, it only consists of a reference
   *        and an index. I do not really see the advantage of returning a shared pointer. Does
   *        it make any difference, here?
   *
   * This is equivalent to \c get_hyEdge.
   *
   * \param   index           The index of the hyperedge to be returned.
   * \retval  hyperedge       Topological information on the hyperedge (cf. \c value_type).
   ************************************************************************************************/
  const value_type operator[](const hyEdge_index_t index) const { return get_hyEdge(index); }
  /*!***********************************************************************************************
   * \brief   Get topological hyperedge of given index.
   *
   * This is equivalent to \c operator[].
   *
   * \param   index           The index of the hyperedge to be returned.
   * \retval  hyperedge       Topological information on the hyperedge (cf. \c value_type).
   ************************************************************************************************/
  value_type get_hyEdge(const hyEdge_index_t index) const
  {
    hy_assert(index < domain_info_.n_hyEdges && index >= 0,
              "Index must be non-negative and smaller than "
                << domain_info_.n_hyEdges << " (which is the amount of hyperedges). It was "
                << index << "!");
    return hyEdge(*this, index);
  }
  /*!***********************************************************************************************
   * \brief   Returns the number of hyperedges making up the hypergraph.
   *
   * \retval  n_hyperedges    The total amount of hyperedges of a hypergraph.
   ************************************************************************************************/
  const hyEdge_index_t n_hyEdges() const { return domain_info_.n_hyEdges; }
  /*!***********************************************************************************************
   * \brief   Returns the number of hypernodes making up the hypergraph.
   *
   * \retval  n_hypernodes    The total amount of hypernodes of a hypergraph.
   ************************************************************************************************/
  const hyNode_index_t n_hyNodes() const { return domain_info_.n_hyNodes; }
  /*!***********************************************************************************************
   * \brief   Returns the whole domain info related to a hypergraph.
   *
   * \retval  domain_indo     Const reference to domain info.
   ************************************************************************************************/
  const DomainInfo<hyEdge_dimT,
                   space_dimT,
                   vectorT,
                   pointT,
                   hyEdge_index_t,
                   hyNode_index_t,
                   pt_index_t>&
  domain_info() const
  {
    return domain_info_;
  }
};  // end of class File

}  // end of namespace Topology
