#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/wrapper/tpcc.hxx>

namespace Topology
{
/*!*************************************************************************************************
 * \brief   Definition of the topology of a hypergraph --- Cubic HyperGraphs.
 *
 * One of the advantages of this software package is the strict discrimination between the topology
 * and the geometry of the domain \f$\Omega\f$. Thus, one can exemplarily define a single topology
 * (the one of a cube) to approximate PDEs that live on the cube's boundary and PDEs that live on a
 * sphere, since their topology is the same. However, different geometries have to be defined, since
 * these obviously are not equal. Thus, all parts of the code that involve communication and/or
 * solving systems of equations are reusable in a much more general (than the standard) sense.
 * Beyond that, absurd (on first sight) domains can be defined easily. This also covers variously
 * periodic domains, for example.
 *
 * \tparam  hyEdge_dimT     Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is
 *                          for PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dimT      The dimension of the space, the object is located in. This number should
 *                          be larger than or equal to hyEdge_dimT.
 * \tparam  NodeIncesVecT   The vector type of an array containing the node indices of an hyperedge.
 * \tparam  ConstructorVecT The vector type of the constructor.
 * \tparam  hyEdge_index_t  The index type of an hyperedge.
 * \tparam NodeOrientationT The class type that encodes the orientation of the hypernodes with
 *                          respect to a given hyperedge.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dimT,
          unsigned int n_subintervalsT = 1,
          typename NodeIndexVecT = SmallVec<2 * hyEdge_dimT, unsigned int>,
          typename ConstructorVecT = SmallVec<space_dimT, unsigned int>,
          typename hyEdge_index_t = typename NodeIndexVecT::value_type,
          typename NodeOrientationT = SmallVec<2 * hyEdge_dimT - 2, unsigned int> >
class CubicRefined
{
  using hyNode_index_t = typename NodeIndexVecT::value_type;
  /*!***********************************************************************************************
   * \brief   Definition of the topology of a hypergraph's edges --- Cubic HyperGraph's edges.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
   * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
   ************************************************************************************************/
  class hyEdge
  {
   private:
    /*!*********************************************************************************************
     * \brief   Indices of the hypernodes adjacent to the hyperedge.
     **********************************************************************************************/
    NodeIndexVecT hyNode_indices_;

   public:
    /*!*********************************************************************************************
     * \brief   Construct a cubic hyperedge from its index and a the global topology class.
     *
     * \param   index               The index of the hyperedge to be created.
     * \param   topology            A cubic topology.
     **********************************************************************************************/
    hyEdge(const hyEdge_index_t index, const Cubic& topology)
    {
      Wrapper::tpcc_elem_t<hyEdge_dimT, space_dimT> elem =
        Wrapper::get_element(topology.tpcc_elements_, index / topology.n_elem_per_elem);
      for (unsigned int i = 0; i < hyNode_indices_.size(); ++i)
      {
        Wrapper::tpcc_elem_t<hyEdge_dimT - 1, space_dimT> face = Wrapper::get_face(elem, i);
        hyNode_indices_[i] = Wrapper::get_index(topology.tpcc_faces_, face);
      }

      Wrapper::tpcc_elem_t<hyEdge_dimT, hyEdge_dimT> ref_elem =
        Wrapper::get_element(topology.tpcc_ref_elem_, index % topology.n_elem_per_elem);
      for (unsigned int i = 0; i < hyNode_indices_.size(); ++i)
      {
        Wrapper::tpcc_elem_t<hyEdge_dimT - 1, hyEdge_dimT> face = Wrapper::get_face(ref_elem, i);
        if (Wrapper::exterior_coordinate(face, 0) == 0 ||
            Wrapper::exterior_coordinate(face, 0) == n_subintervalsT + 1)
          hyNode_indices_[i] = n_face_per_face * hyNode_indices_[i] +
                               Wraper::get_index_in_slice(topology.tpcc_ref_faces_, face);
        else
          hyNode_indices_[i] = n_face_per_elem * index(topology.tpcc_ref_elem_, ref_elem) +
                               Wrapper::get_index(topology.tpcc_ref_faces_, face);
      }
    }
    /*!*********************************************************************************************
     * \brief   Return indices of hypernodes adjacent to the hyperedge.
     *
     * \retval  hypernode_indices   Topological information on the hyperedge (cf. \c value_type).
     **********************************************************************************************/
    const NodeIndexVecT& get_hyNode_indices() const { return hyNode_indices_; }
    /*!*********************************************************************************************
     * \brief   Return orienation of hypernode.
     **********************************************************************************************/
    const NodeOrientationT get_hyNode_oriantation(unsigned int) const { return NodeOrientationT(); }
  };  // end of class hyEdge

 public:
  /*!***********************************************************************************************
   * \brief   Return the template parameter representing the dimension of a hyperedge.
   *
   * \retval  hyEdge_dimT     The dimension of a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; };
  /*!***********************************************************************************************
   * \brief   Return the template parameter representing the dimension of the space.
   *
   * \retval  space_dimT      The dimension of the space.
   ************************************************************************************************/
  static constexpr unsigned int space_dim() { return space_dimT; };

 private:
  /*!***********************************************************************************************
   * \brief   Number of elements per spatial dimension.
   ************************************************************************************************/
  const ConstructorVecT n_elements_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for elements.
   ************************************************************************************************/
  const Wrapper::tpcc_t<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyNode_index_t>
    tpcc_elements_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for faces.
   ************************************************************************************************/
  const Wrapper::tpcc_t<hyEdge_dimT - 1, space_dimT, TPCC::boundaries::both, hyNode_index_t>
    tpcc_faces_;

  const Wrapper::tpcc_t<hyEdge_dimT, hyEdge_dimT, hyNode_index_t> tpcc_ref_elem_;
  const Wrapper::tpcc_t<hyEdge_dimT - 1, hyEdge_dimT, hyNode_index_t> tpcc_ref_faces_;

  const unsigned int n_elem_per_elem, n_face_per_face, n_face_per_elem, n_coarse_elem,
    n_coarse_face;

  /*!***********************************************************************************************
   * \brief   Total amount of hyperedges.
   ************************************************************************************************/
  const hyEdge_index_t n_hyEdges_;
  /*!***********************************************************************************************
   * \brief   Total amount of hypernodes.
   ************************************************************************************************/
  const hyNode_index_t n_hyNodes_;

 public:
  /*!***********************************************************************************************
   * \brief   Define the return value of the class.
   ************************************************************************************************/
  typedef hyEdge value_type;
  /*!***********************************************************************************************
   * \brief   Define the value type of input argument for standard constructor.
   ************************************************************************************************/
  typedef ConstructorVecT constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Construct a cubic hypergraph.
   *
   * \param   n_elements      A vector / array containing number of elements per dimension.
   ************************************************************************************************/
  CubicRefined(const constructor_value_type& n_elements)
  : n_elements_(n_elements),
    tpcc_elements_(
      Wrapper::create_tpcc<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        n_elements)),
    tpcc_faces_(Wrapper::tpcc_faces<TPCC::boundaries::both>(tpcc_elements_)),
    tpcc_ref_elem_(
      Wrapper::create_tpcc<hyEdge_dimT, hyEdge_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        SmallVec<hyEdge_dimT, unsigned int>(n_subintervalsT))),
    tpcc_ref_faces_(Wrapper::tpcc_faces<TPCC::boundaries::none>(tpcc_ref_elem_)),

    n_elem_per_elem(Wrapper::n_elements(tpcc_ref_elem_)),
    n_face_per_face(Hypercube<hyEdge_dimT - 1>::pow(n_subintervalsT)),
    n_face_per_elem(Wrapper::n_elements(tpcc_ref_faces_)),
    n_coarse_elem(Wrapper::n_elements(tpcc_elements_)),
    c_coarse_face(Wrapper::n_elements(tpcc_faces_)),
    n_hyEdges_(n_coarse_elem * n_elem_per_elem),
    n_hyNodes_(n_coarse_face * n_face_per_face + n_coarse_elem * n_face_per_elem)
  {
  }
  /*!***********************************************************************************************
   * \brief   Construct a hypergraph topology from another hypergraph topology.
   *
   * Create a (value based) copy of another hypergraph.
   *
   * \param   other           Hypergraph to be copied.
   ************************************************************************************************/
  CubicRefined(const Cubic<hyEdge_dimT, space_dimT>& other)
  : n_elements_(other.n_elements_),
    tpcc_elements_(other.tpcc_elements_),
    tpcc_faces_(other.tpcc_faces_),
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
   * This function returns the hyperedge of the given index, i.e., it returns the topological
   * hyperedge (\b not the geometrical information). The topological informatiom comprises the
   * indices of adjacent hypernodes and information about their respective orientations.
   *
   * \param   index           The index of the hyperedge to be returned.
   * \retval  hyperedge       Topological information on the hyperedge (cf. \c value_type).
   ************************************************************************************************/
  const value_type operator[](const hyEdge_index_t index) const
  {
    hy_assert(index >= 0 && index < n_hyEdges_,
              "The index of an hyperedge must be non-negative and smaller than the total amount "
                << "of hyperedges, which is " << n_hyEdges_ << ". Nonetheless, the " << index
                << "-th hyperedge is tried to be accessed.");
    return hyEdge(index, *this);
  }
  /*!***********************************************************************************************
   * \brief   Read the array of elements per dimensions.
   *
   * \retval  n_elements      A vector / arary containing the elements in the repective dimension.
   ************************************************************************************************/
  const ConstructorVecT& n_elements() const { return n_elements_; }
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for elements.
   ************************************************************************************************/
  const Wrapper::tpcc_t<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyNode_index_t>&
  tpcc_elem() const
  {
    return tpcc_elements_;
  }
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for faces.
   ************************************************************************************************/
  const Wrapper::tpcc_t<hyEdge_dimT - 1, space_dimT, TPCC::boundaries::both, hyNode_index_t>&
  tpcc_face() const
  {
    return tpcc_faces_;
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
};  // end of class CubicRefined

}  // end of namespace Topology
