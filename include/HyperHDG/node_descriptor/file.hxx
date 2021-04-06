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
          template <typename...> typename vectorT = std::vector,
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

    /*!*********************************************************************************************
     * \brief   Indices of the hypernodes adjacent to the hyperedge.
     **********************************************************************************************/
    std::array<hyNode_index_t, 2 * hyEdge_dimT> hyFace_types_;

   public:
    /*!*********************************************************************************************
     * \brief   Construct hyperedge from hypergraph and index.
     **********************************************************************************************/
    hyEdge(const File& node_desc, const hyEdge_index_t index)
    : hyGraph_topology_(node_desc), index_(index)
    {
      hyFace_types_ = node_desc.domain_info_
                        .hyFaces_hyEdge[index / Wrapper::n_elements(node_desc.tpcc_ref_elem_)];

      Wrapper::tpcc_elem_t<hyEdge_dimT, hyEdge_dimT> ref_elem = Wrapper::get_element(
        node_desc.tpcc_ref_elem_, index % Wrapper::n_elements(node_desc.tpcc_ref_elem_));
      for (unsigned int i = 0; i < hyFace_types_.size(); ++i)
      {
        Wrapper::tpcc_elem_t<hyEdge_dimT - 1, hyEdge_dimT> face = Wrapper::get_face(ref_elem, i);
        if (Wrapper::exterior_coordinate(face, 0) != 0 &&
            Wrapper::exterior_coordinate(face, 0) != node_desc.n_subintervals_)
          hyFace_types_[i] = 0;
      }
    }
    /*!*********************************************************************************************
     * \brief   Return hypernode types of the \c hyEdge.
     **********************************************************************************************/
    const auto& get_hyFaces_types() const { return hyFace_types_; }
    /*!*********************************************************************************************
     * \brief   Return hypernode type of given index.
     **********************************************************************************************/
    unsigned int operator[](const unsigned int index) const { return hyFace_types_[index]; }
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

  /*!***********************************************************************************************
   * \brief   Refinment level corresponds to number of subintervals per dimension.
   ************************************************************************************************/
  unsigned int n_subintervals_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for refined elements of hypergedge.
   ************************************************************************************************/
  Wrapper::tpcc_t<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyNode_index_t> tpcc_ref_elem_;
  /*!***********************************************************************************************
   * \brief   Total amount of hyperedges in hypergraph.
   ************************************************************************************************/
  hyEdge_index_t n_hyEdges_;

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
  File(const constructor_value_type& topology)
  : domain_info_(topology.domain_info()),
    n_subintervals_(topology.get_refinement()),
    tpcc_ref_elem_(topology.tpcc_ref_elem()),
    n_hyEdges_(topology.n_hyEdges())
  {
  }
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
  : domain_info_(other.domain_info),
    n_subintervals_(1),
    tpcc_ref_elem_(
      Wrapper::create_tpcc<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        SmallVec<hyEdge_dimT, unsigned int>(n_subintervals_))),
    n_hyEdges_(domain_info_.n_hyEdges)
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
    hy_assert(index < domain_info_.n_hyEdges * Wrapper::n_elements(tpcc_ref_elem_) && index >= 0,
              "Index must be non-negative and smaller than "
                << domain_info_.n_hyEdges << " (which is the amount of hyperedges). It was "
                << index << "!");
    return hyEdge(*this, index);
  }
  /*!***********************************************************************************************
   * \brief   Return the refinement level (equal to number of subintervals).
   ************************************************************************************************/
  unsigned int get_refinement() const { return n_subintervals_; }
  /*!***********************************************************************************************
   * \brief   Set the refinement level (equal to number of subintervals).
   ************************************************************************************************/
  void set_refinement(unsigned int level)
  {
    n_subintervals_ = level;
    tpcc_ref_elem_ =
      Wrapper::create_tpcc<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        SmallVec<hyEdge_dimT, unsigned int>(n_subintervals_));
    n_hyEdges_ = domain_info_.n_hyEdges * Wrapper::n_elements(tpcc_ref_elem_);
  }
};  // end of class File

}  // namespace NodeDescriptor
