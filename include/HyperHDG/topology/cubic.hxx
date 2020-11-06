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
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dimT,
          typename NodeIndexVecT = SmallVec<2 * hyEdge_dimT, unsigned int>,
          typename ConstructorVecT = SmallVec<space_dimT, unsigned int>,
          typename hyEdge_index_t = typename NodeIndexVecT::value_type,
          typename NodeOrientationT = SmallVec<2 * hyEdge_dimT - 2, unsigned int> >
class Cubic
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
     *
     * A \c std::array comprising the indices of the hypernodes adjacent to a hyperedge.
     **********************************************************************************************/
    NodeIndexVecT hyNode_indices_;

   public:
    /*!*********************************************************************************************
     * \brief   Construct a cubic hyperedge from its index and a \c std::array of elements in each
     *          spatial dimension.
     *
     * Constructs a hyperedge from a \c std::array containing the elementens per spatial dimension
     * which is given as input data and the index of the hyperedge to be constructed.
     *
     * \param   index           The index of the hyperedge to be created.
     * \param   topology        A cubic topology.
     **********************************************************************************************/
    hyEdge(const hyEdge_index_t index, const Cubic& topology)
    {
      Wrapper::tpcc_elem_t<hyEdge_dimT, space_dimT> elem =
        Wrapper::get_element<hyEdge_dimT, space_dimT, unsigned int>(topology.tpcc_elements_, index);
      for (unsigned int i = 0; i < hyNode_indices_.size(); ++i)
      {
        Wrapper::tpcc_elem_t<hyEdge_dimT - 1, space_dimT> face =
          Wrapper::get_face<hyEdge_dimT, space_dimT>(elem, i);
        hyNode_indices_[i] =
          Wrapper::get_index<hyEdge_dimT - 1, space_dimT, unsigned int>(topology.tpcc_faces_, face);
      }
    }
    /*!*********************************************************************************************
     * \brief   Return indices of hypernodes adjacent to the hyperedge.
     *
     * Return a \c std::array comprising the indices of the hypernodes adjacent to a hyperedge.
     *
     * \retval  hypernode_indeices  Topological information on the hyperedge (cf. \c value_type).
     **********************************************************************************************/
    const NodeIndexVecT& get_hyNode_indices() const { return hyNode_indices_; }
    /*!*********************************************************************************************
     * \brief   Return orienation of hypernode.
     **********************************************************************************************/
    const NodeOrientationT get_hyNode_oriantation(unsigned int node) const
    {
      return NodeOrientationT();
    }
  };  // end of class hyEdge

 public:
  /*!***********************************************************************************************
   * \brief   Returns the template parameter representing the dimension of a hyperedge.
   *
   * \retval  hyEdge_dimT   The dimension of a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; };
  /*!***********************************************************************************************
   * \brief   Returns the template parameter representing the dimension of the space.
   *
   * \retval  space_dimT       The dimension of the space.
   ************************************************************************************************/
  static constexpr unsigned int space_dim() { return space_dimT; };

 private:
  /*!***********************************************************************************************
   * \brief   Number of elements per spatial dimension.
   *
   * A \c std::array comprising the number of elements in each spatial dimension.
   ************************************************************************************************/
  const ConstructorVecT n_elements_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for elements.
   ************************************************************************************************/
  const Wrapper::tpcc_t<hyEdge_dimT, space_dimT, hyNode_index_t> tpcc_elements_;
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for faces.
   ************************************************************************************************/
  const Wrapper::tpcc_t<hyEdge_dimT - 1, space_dimT, hyNode_index_t> tpcc_faces_;
  /*!***********************************************************************************************
   * \brief   Total amount of hyperedges.
   *
   * The number of hyperedges that form the hypergraph. This information is needed to allow to go
   * through all hyperedges and execute some code. The number of hyperedges can be computed from
   * the \c std::array \c n_elements_.
   ************************************************************************************************/
  const hyEdge_index_t n_hyEdges_;
  /*!***********************************************************************************************
   * \brief   Total amount of hypernodes.
   *
   * The number of hypernodes that make up the hypergraph. This information is needed to have the
   * appropriate version of a \c HyperNodeFactory. It can be vomputed from the \c std::array
   * \c n_elements_.
   ************************************************************************************************/
  const hyNode_index_t n_hyNodes_;

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
  typedef ConstructorVecT constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Construct a cubic hypergraph from a \c std::vector.
   *
   * Constructs a hypergraph from a \c std::vector containing the elementens per spatial dimension
   * which is given as input data. If the input vector is shorter that \c space_dimT, the
   * remaining amounts of elemnts are assumed to be equal to zero. If the vector is longer than
   * \c space_dimT, the first \c space_dimT entries are considered only.
   *
   * \todo    If the vector is too short, an error is thrown in the test program and the behavior
   *          is undefined for Python (most likely an error is thrown, too) at the moment.
   *
   * \param   n_elements    A \c std::vector containing number of elements per dimension.
   ************************************************************************************************/
  Cubic(const constructor_value_type& n_elements)
  : n_elements_(n_elements),
    tpcc_elements_(Wrapper::create_tpcc<hyEdge_dimT, space_dimT, hyEdge_index_t>(n_elements)),
    tpcc_faces_(Wrapper::tpcc_faces<hyEdge_dimT, space_dimT, hyEdge_index_t>(tpcc_elements_)),
    n_hyEdges_(Wrapper::n_elements<hyEdge_dimT, space_dimT, hyEdge_index_t>(tpcc_elements_)),
    n_hyNodes_(Wrapper::n_elements<hyEdge_dimT - 1, space_dimT, hyEdge_index_t>(tpcc_faces_))
  {
  }
  /*!***********************************************************************************************
   * \brief   Construct a hypergraph from another hypergraph.
   *
   * \todo    Guido: If possible, this function computes the amount of hyperedges and hypernodes.
   *
   * Create a (value based) copy of another hypergraph.
   *
   * \param   other           Hypergraph to be copied.
   ************************************************************************************************/
  Cubic(const Cubic<hyEdge_dimT, space_dimT>& other)
  : n_elements_(other.n_elements_),
    tpcc_elements_(other.tpcc_elements_),
    tpcc_faces_(other.tpcc_faces_),
    n_hyEdges_(other.n_hyEdges_),
    n_hyNodes_(other.n_hyNodes_)
  {
  }
  /*!***********************************************************************************************
   * \brief   Get topological hyperedge of given index.
   *
   * \todo  Here, we repeatedly return a large object. This is done since the object could be
   *        locally created in regular topologies/geometries! Return shared-pointer?
   *        -> AR: I do not really see the point, but might just be stupid ;)
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
   * \retval  n_elements    A \c std::array containing the elements in the repective dimension.
   ************************************************************************************************/
  const ConstructorVecT& n_elements() const { return n_elements_; }
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for elements.
   ************************************************************************************************/
  Wrapper::tpcc_t<hyEdge_dimT, space_dimT, hyNode_index_t> tpcc_elem() const
  {
    return tpcc_elements_;
  }
  /*!***********************************************************************************************
   * \brief   Tensor product chain complex for faces.
   ************************************************************************************************/
  Wrapper::tpcc_t<hyEdge_dimT - 1, space_dimT, hyNode_index_t> tpcc_face() const
  {
    return tpcc_faces_;
  }
  /*!***********************************************************************************************
   * \brief   Returns the number of hyperedges making up the hypergraph.
   *
   * \retval  n_hyperedges    The total amount of hyperedges of a hypergraph.
   ************************************************************************************************/
  const hyEdge_index_t n_hyEdges() const { return n_hyEdges_; }
  /*!***********************************************************************************************
   * \brief   Returns the number of hypernodes making up the hypergraph.
   *
   * \retval  n_hypernodes    The total amount of hypernodes of a hypergraph.
   ************************************************************************************************/
  const hyNode_index_t n_hyNodes() const { return n_hyNodes_; }
};  // end of class Cubic

}  // end of namespace Topology
