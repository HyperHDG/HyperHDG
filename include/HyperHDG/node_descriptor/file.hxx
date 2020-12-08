#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/topology/file.hxx>

#include <string>

namespace NodeDescriptor
{
/*!*************************************************************************************************
 * \brief   Hypergraph topology based on an input file.
 *
 * The node descriptor class of a hypergraph which is read from a file. Thus, this class uses the
 * file's information on hypernode types and returns them.
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
   * \brief   Definition of the node types of a hypergraph's edges.
   ************************************************************************************************/
  class hyEdge
  {
   public:
    static constexpr unsigned int n_hyNodes() { return 2 * hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Return dimension of the hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Return dimension of the surrounding space.
     **********************************************************************************************/
    static constexpr unsigned int space_dim() { return space_dimT; }

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
     * \brief   Return hypernode types of the \c hyEdge.
     **********************************************************************************************/
    const auto& get_hyFaces_types() const
    {
      return hyGraph_topology_.domain_info_.hyFaces_hyEdge[index_];
    }
    /*!*********************************************************************************************
     * \brief   Return hypernode type of given index.
     **********************************************************************************************/
    unsigned int operator[](const unsigned int index) const
    {
      return hyGraph_topology_.domain_info_.hyFaces_hyEdge[index_][index];
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
                   pt_index_t>& domain_info_;

 public:
  /*!***********************************************************************************************
   * \brief   Define the return value of the class.
   ************************************************************************************************/
  typedef hyEdge value_type;
  /*!***********************************************************************************************
   * \brief   Defines the value type of input argument for standard constructor.
   ************************************************************************************************/
  typedef Topology::
    File<hyEdge_dimT, space_dimT, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
      constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Construct a node descriptor from a file topology.
   ************************************************************************************************/
  File(const constructor_value_type& topology) : domain_info_(topology.domain_info()) {}
  /*!***********************************************************************************************
   * \brief   Copy constructor.
   ************************************************************************************************/
  File(const File<hyEdge_dimT,
                  space_dimT,
                  vectorT,
                  pointT,
                  hyEdge_index_t,
                  hyNode_index_t,
                  pt_index_t>& other)
  : domain_info_(other.domain_info)
  {
  }

  /*!***********************************************************************************************
   * \brief   Return hyperegde of given index.
   ************************************************************************************************/
  value_type operator[](const hyEdge_index_t index) const { return get_hyEdge(index); }
  /*!***********************************************************************************************
   * \brief   Return hyperedge of given index.
   ************************************************************************************************/
  value_type get_hyEdge(const hyEdge_index_t index) const
  {
    hy_assert(index < domain_info_.n_hyEdges && index >= 0,
              "Index must be non-negative and smaller than "
                << domain_info_.n_hyEdges << " (which is the amount of hyperedges). It was "
                << index << "!");
    return hyEdge(*this, index);
  }
};  // end of class File

}  // namespace NodeDescriptor
