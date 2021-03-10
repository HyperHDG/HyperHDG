#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/read_domain.hxx>
#include <HyperHDG/wrapper/tpcc.hxx>

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
 * \tparam  vectorT         The typename of the large vector type. Defaults to std::vector.
 * \tparam  pointT          The typename of a Point class. Defaults to \c Point<space_dimT, float>.
 * \tparam  hyEdge_index_t  The index type for hyperedges. Default is \c unsigned \c int.
 * \tparam  hyNode_index_t  The index type for hypernodes. Default is \c hyNode_index_t.
 * \tparam  pt_index_t      The index type of points. Default is \c hyNode_index_t.
 * \tparam NodeOrientationT The class type that encodes the orientation of the hypernodes with
 *                          respect to a given hyperedge.
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
          typename pt_index_t = hyNode_index_t,
          typename NodeOrientationT = SmallVec<2 * hyEdge_dimT - 2, unsigned int> >
class File
{
  /*!***********************************************************************************************
   * \brief   Definition of the topology of a hypergraph's edges.
   ************************************************************************************************/
  class hyEdge
  {
   public:
    /*!*********************************************************************************************
     * \brief   Number of hypernodes per hyperedge.
     **********************************************************************************************/
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

    /*!*********************************************************************************************
     * \brief   Indices of the hypernodes adjacent to the hyperedge.
     **********************************************************************************************/
    std::array<hyNode_index_t, 2 * hyEdge_dimT> hyNode_indices_;

   public:
    /*!*********************************************************************************************
     * \brief   Construct hyperedge from hypergraph and index.
     **********************************************************************************************/
    hyEdge(const File& topology, const hyEdge_index_t index)
    : hyGraph_topology_(topology), index_(index)
    {
      hyNode_indices_ =
        hyGraph_topology_.domain_info_.hyNodes_hyEdge[index / topology.n_elem_per_elem];

      Wrapper::tpcc_elem_t<hyEdge_dimT, hyEdge_dimT> ref_elem =
        Wrapper::get_element(topology.tpcc_ref_elem_, index % topology.n_elem_per_elem);
      for (unsigned int i = 0; i < hyNode_indices_.size(); ++i)
      {
        Wrapper::tpcc_elem_t<hyEdge_dimT - 1, hyEdge_dimT> face = Wrapper::get_face(ref_elem, i);
        if (Wrapper::exterior_coordinate(face, 0) == 0 ||
            Wrapper::exterior_coordinate(face, 0) == topology.n_subintervals_)
          hyNode_indices_[i] = topology.n_face_per_face * hyNode_indices_[i] +
                               Wrapper::get_index_in_slice(topology.tpcc_ref_faces_, face);
        else
          hyNode_indices_[i] = topology.n_coarse_face * topology.n_face_per_face +
                               topology.n_face_per_elem * (index / topology.n_elem_per_elem) +
                               Wrapper::get_index(topology.tpcc_ref_faces_, face);
      }
    }
    /*!*********************************************************************************************
     * \brief   Return hypernodes of a hyperedge.
     **********************************************************************************************/
    const auto& get_hyNode_indices() const { return hyNode_indices_; }
    /*!*********************************************************************************************
     * \brief   Return orienation of hypernode.
     **********************************************************************************************/
    const NodeOrientationT get_hyNode_oriantation(unsigned int) const { return NodeOrientationT(); }
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

  /*!***********************************************************************************************
   * \brief   Refinment level corresponds to number of subintervals per dimension.
   ************************************************************************************************/
  unsigned int n_subintervals_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for refining a hyperedge into elements.
   ************************************************************************************************/
  Wrapper::tpcc_t<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyNode_index_t> tpcc_ref_elem_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for refining a hyperedge into faces.
   ************************************************************************************************/
  Wrapper::tpcc_t<hyEdge_dimT - 1, hyEdge_dimT, TPCC::boundaries::none, hyNode_index_t>
    tpcc_ref_faces_;
  /*!***********************************************************************************************
   * \brief   Number of refined elements per corase element.
   ************************************************************************************************/
  unsigned int n_elem_per_elem;
  /*!***********************************************************************************************
   * \brief   Number of refined faces per corase face.
   ************************************************************************************************/
  unsigned int n_face_per_face;
  /*!***********************************************************************************************
   * \brief   Number of refined faces per corase element.
   ************************************************************************************************/
  unsigned int n_face_per_elem;
  /*!***********************************************************************************************
   * \brief   Number of refined corase elements.
   ************************************************************************************************/
  unsigned int n_coarse_elem;
  /*!***********************************************************************************************
   * \brief   Number of refined corase faces.
   ************************************************************************************************/
  unsigned int n_coarse_face;

  /*!***********************************************************************************************
   * \brief   Total amount of hyperedges.
   ************************************************************************************************/
  hyEdge_index_t n_hyEdges_;
  /*!***********************************************************************************************
   * \brief   Total amount of hypernodes.
   ************************************************************************************************/
  hyNode_index_t n_hyNodes_;

 public:
  /*!***********************************************************************************************
   * \brief   Defines the return value of the class.
   ************************************************************************************************/
  typedef hyEdge value_type;
  /*!***********************************************************************************************
   * \brief   Defines the value type of input argument for standard constructor.
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
                             pt_index_t>(filename)),
    n_subintervals_(1),
    tpcc_ref_elem_(
      Wrapper::create_tpcc<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        SmallVec<hyEdge_dimT, unsigned int>(n_subintervals_))),
    tpcc_ref_faces_(Wrapper::tpcc_faces<TPCC::boundaries::none>(tpcc_ref_elem_)),
    n_elem_per_elem(Wrapper::n_elements(tpcc_ref_elem_)),
    n_face_per_face(Hypercube<hyEdge_dimT - 1>::pow(n_subintervals_)),
    n_face_per_elem(Wrapper::n_elements(tpcc_ref_faces_)),
    n_coarse_elem(domain_info_.n_hyEdges),
    n_coarse_face(domain_info_.n_hyNodes),
    n_hyEdges_(n_coarse_elem * n_elem_per_elem),
    n_hyNodes_(n_coarse_face * n_face_per_face + n_coarse_elem * n_face_per_elem)
  {
  }
  /*!***********************************************************************************************
   * \brief   Copy constructor.
   ************************************************************************************************/
  File(const File<hyEdge_dimT, space_dimT>& other)
  : domain_info_(other.domain_info()),
    n_subintervals_(other.get_refinement()),
    tpcc_ref_elem_(other.tpcc_ref_elem_),
    tpcc_ref_faces_(other.tpcc_ref_faces_),
    n_elem_per_elem(other.n_elem_per_elem),
    n_face_per_face(other.n_face_per_face),
    n_face_per_elem(other.n_face_per_elem),
    n_coarse_elem(other.n_coarse_elem),
    n_coarse_face(other.n_coarse_face),
    n_hyEdges_(other.n_hyEdges_),
    n_hyNodes_(other.n_hyNodes_)
  {
  }

  /*!***********************************************************************************************
   * \brief   Get topological hyperedge of given index.
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
  const value_type get_hyEdge(const hyEdge_index_t index) const
  {
    hy_assert(index < n_hyEdges_ && index >= 0, "Index must be non-negative and smaller than "
                                                  << domain_info_.n_hyEdges
                                                  << " (which is the amount of hyperedges). It was "
                                                  << index << "!");
    return hyEdge(*this, index);
  }
  /*!***********************************************************************************************
   * \brief   Return the number of hyperedges making up the hypergraph.
   *
   * \retval  n_hyperedges    The total amount of hyperedges of a hypergraph.
   ************************************************************************************************/
  const hyEdge_index_t n_hyEdges() const { return n_hyEdges_; }
  /*!***********************************************************************************************
   * \brief   Return the number of hypernodes making up the hypergraph.
   *
   * \retval  n_hypernodes    The total amount of hypernodes of a hypergraph.
   ************************************************************************************************/
  const hyNode_index_t n_hyNodes() const { return n_hyNodes_; }
  /*!***********************************************************************************************
   * \brief   Return the whole domain info related to a hypergraph.
   *
   * \retval  domain_info     Const reference to domain info.
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

  const Wrapper::tpcc_t<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyNode_index_t>&
  tpcc_ref_elem() const
  {
    return tpcc_ref_elem_;
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
    tpcc_ref_faces_ = Wrapper::tpcc_faces<TPCC::boundaries::none>(tpcc_ref_elem_);
    n_elem_per_elem = Wrapper::n_elements(tpcc_ref_elem_);
    n_face_per_face = Hypercube<hyEdge_dimT - 1>::pow(n_subintervals_);
    n_face_per_elem = Wrapper::n_elements(tpcc_ref_faces_);
    n_hyEdges_ = n_coarse_elem * n_elem_per_elem;
    n_hyNodes_ = n_coarse_face * n_face_per_face + n_coarse_elem * n_face_per_elem;
  }
};  // end of class File

}  // end of namespace Topology
