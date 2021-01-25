#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/topology/cubic.hxx>
#include <cmath>

namespace NodeDescriptor
{
/*!*************************************************************************************************
 * \brief   Definition of the node types of a cubic hypergraph.
 *
 * In a cubic hypergraph \f$\mathcal G\f$, the nodes receive their indices according to their
 * positions. Thus, the node type \f$t\f$ of node \f$\mathcal N\f$ is defined as
 * \f[
 *  t = \sum_{n=1}^\text{space_dim} 3^n \left[ ( \Pi_n \mathcal N == \{ \min \Pi_n \mathcal G \} )
 *      + 2 ( \Pi_n \mathcal N == \{ \max \Pi_n \mathcal G \} ) \right],
 * \f]
 * where \f$\Pi_n\f$ is the orthonormal projection onto the \f$n\f$-th axis of the canonical basis.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dimT,
          typename NodeTypeVecT = SmallVec<2 * hyEdge_dimT, unsigned int>,
          typename ConstructorVecT = SmallVec<space_dimT, unsigned int>,
          typename hyEdge_index_t = unsigned int,
          typename hyNode_index_t = hyEdge_index_t>
class Cubic
{
  /*!***********************************************************************************************
   * \brief   Definition of the node types of a hypergraph's edges.
   *
   * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
   * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
   ************************************************************************************************/
  class hyEdge
  {
   private:
    /*!*********************************************************************************************
     * \brief   Types of the hypernodes adjacent to the hyperedge.
     **********************************************************************************************/
    NodeTypeVecT hyFace_types_;

   public:
    /*!*********************************************************************************************
     * \brief   Construct a \c hyperEdge, i.e. determine its node indices.
     *
     * \param   index           The index of the hyperedge to be created.
     * \param   node_desc       Cubic node descriptor.
     **********************************************************************************************/
    hyEdge(const hyEdge_index_t index, const Cubic& node_desc)
    {
      Wrapper::tpcc_elem_t<hyEdge_dimT, space_dimT> elem =
        Wrapper::get_element<hyEdge_dimT, space_dimT, TPCC::boundaries::both, unsigned int>(
          node_desc.tpcc_elements_, index);
      for (unsigned int i = 0; i < hyFace_types_.size(); ++i)
      {
        Wrapper::tpcc_elem_t<hyEdge_dimT - 1, space_dimT> face =
          Wrapper::get_face<hyEdge_dimT, space_dimT>(elem, i);
        hyFace_types_[i] = 0;
        for (unsigned int dim = 0; dim < space_dimT - hyEdge_dimT + 1; ++dim)
        {
          unsigned int coordinate =
            Wrapper::exterior_coordinate<hyEdge_dimT - 1, space_dimT>(face, dim);
          unsigned int direction =
            Wrapper::exterior_direction<hyEdge_dimT - 1, space_dimT>(face, dim);
          if (coordinate == 0 || coordinate == node_desc.n_elements().operator[](direction))
            hyFace_types_[i] += (1 + (coordinate != 0)) * std::pow(3, direction);
        }
      }
    }
    /*!*********************************************************************************************
     * \brief   Return types of hypernodes adjacent to the hyperedge.
     **********************************************************************************************/
    const NodeTypeVecT& get_hyFace_types() const { return hyFace_types_; }
    /*!*********************************************************************************************
     * \brief   Return types of hypernode of given index.
     **********************************************************************************************/
    unsigned int operator[](const unsigned int index) const { return hyFace_types_[index]; }
  };  // end of class hyEdge

 public:
  /*!***********************************************************************************************
   * \brief   Return the template parameter representing the dimension of a hyperedge.
   *
   * \retval  hyEdge_dimT   The dimension of a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; };
  /*!***********************************************************************************************
   * \brief   Return the template parameter representing the dimension of the space.
   *
   * \retval  space_dimT    The dimension of the space.
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
  /*!***********************************************************************************************
   * \brief   Total amount of hyperedges in hypergraph.
   ************************************************************************************************/
  const hyEdge_index_t n_hyEdges_;

 public:
  /*!***********************************************************************************************
   * \brief   Defines the return value of the class to be a \c hyEdge.
   ************************************************************************************************/
  typedef hyEdge value_type;
  /*!***********************************************************************************************
   * \brief   Defines the value type of input argument for standard constructor.
   ************************************************************************************************/
  typedef ConstructorVecT constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Construct a node descriptor of a cubic hypergraph.
   *
   * \param   n_elements    A vector/array containing number of elements per spatial dimension.
   ************************************************************************************************/
  Cubic(const ConstructorVecT& n_elements)
  : n_elements_(n_elements),
    tpcc_elements_(
      Wrapper::create_tpcc<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        n_elements)),
    tpcc_faces_(
      Wrapper::tpcc_faces<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyEdge_index_t>(
        tpcc_elements_)),
    n_hyEdges_(Wrapper::n_elements<hyEdge_dimT, space_dimT, TPCC::boundaries::both, hyEdge_index_t>(
      tpcc_elements_))
  {
  }
  /*!***********************************************************************************************
   * \brief   Construct a cubic node descruiptor from a cubic topology.
   *
   * \param   other         Topology to which the node descriptor will fit.
   ************************************************************************************************/
  Cubic(const Topology::Cubic<hyEdge_dimT, space_dimT>& other)
  : n_elements_(other.n_elements()),
    tpcc_elements_(other.tpcc_elem()),
    tpcc_faces_(other.tpcc_face()),
    n_hyEdges_(other.n_hyEdges())
  {
  }
  /*!***********************************************************************************************
   * \brief   Get topological hyperedge of given index.
   *
   * This function returns the hyperedge of the given index, i.e., it returns the  hyperedge 's node
   * description info.
   *
   * \param   index         The index of the hyperedge to be returned.
   * \retval  hyperedge     Information on the hyperedge's nodes (cf. \c value_type).
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
   * \brief   Return the array / vector of elements per dimensions.
   *
   * \retval  n_elements    A vector / array containing the elements in the repective dimension.
   ************************************************************************************************/
  const ConstructorVecT& n_elements() const { return n_elements_; }
  /*!***********************************************************************************************
   * \brief   Returns the number of hyperedges making up the hypergraph.
   *
   * \retval  n_hyperedges  The total amount of hyperedges of a hypergraph.
   ************************************************************************************************/
  const hyEdge_index_t n_hyEdges() const { return n_hyEdges_; }
};  // end of class Cubic

}  // end of namespace NodeDescriptor
