#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/HyAssert.hxx>
#include <HyperHDG/HyperNodeFactory.hxx>
#include <HyperHDG/HyDataContainer.hxx>

#include <memory>

class EmptyC{};

/*!*************************************************************************************************
 * \brief   The class template uniting topology and geometry of a hypergraph with the topology of
 *          the skeleton space of the HDG method.
 * 
 * The main class representing a hypergraph. It uses a class \c Topology to represent the collection
 * of nodes and edges as well as a class \c Geometry presenting the physical coordinates of the
 * edges. It behaves like a random access container of hyperedges and has additional access to its
 * nodes.
 *
 * In our abstraction, nodes only carry degrees of freedom. Thus, they can be obtained from one
 * object \c HyperNodeFactory for any graph. Their location, if such a notion is reasonable, must be
 * determined by that of the boundaries of an edge. The meaning of their degrees of freedom is
 * decided by the local solvers of the HDG method applied. The \c Geometry class may use degrees of
 * freedom of the nodes as well.
 *
 * \tparam  n_dofs_per_nodeT The number of degrees of freedom of a single hypernode which is assumed
 *                          to be the same for all hypernodes.
 * \tparam  TopoT           Class that contains the topology of the hypergraph. This class is needs
 *                          to provide a getter function to the topological information of a
 *                          hyperedge of given index and can be arbitrarily implemented.
 * \tparam  GeomT           Class that contains the topology of the hypergraph. This class is needs
 *                          to provide a getter function to the topological information of a
 *                          hyperedge of given index and can be arbitrarily implemented.
 * \tparam  NodeDescT       Class that contains the descriptions of hyperedges' faces.
 * \tparam  hyEdge_index_t  Unsigned integer type specification. Default is unsigned int.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  unsigned int n_dofs_per_nodeT, class TopoT, class GeomT, class NodeT, class DataT = EmptyC,
  typename hyEdge_index_t = unsigned int
>
class HDGHyperGraph
{
  
  /*!***********************************************************************************************
   * \brief   The type for a hyperedge returned by \c operator[].
   *
   * This \c typedef \c struct is returned by the \c operator[] of an \c HDGHyperGraph. It contains
   * topological and geometrical information about a single hyperedge. It is therefore defined as
   * the \c value_type of class \c HDGHyperGraph.
   ************************************************************************************************/
  typedef struct hyEdge
  {
    /*!*********************************************************************************************
     * \brief   Topological information of a hyperedge.
     *
     * A \c TopoT::value_type comprising the topological information of a hyperedge.
     **********************************************************************************************/
    typename TopoT::value_type topology;
    /*!*********************************************************************************************
     * \brief   Geometrical information of a hyperedge.
     *
     * A \c TopoT::value_type comprising the geometrical information of a hyperedge.
     **********************************************************************************************/
    typename GeomT::value_type geometry;
    /*!*********************************************************************************************
     * \brief   Type of nodes information of a hyperedge.
     *
     * A \c NodeDescT::value_type comprising the information about faces of a hyperedge.
     **********************************************************************************************/
    typename NodeT::value_type node_descriptor;
    /*!*********************************************************************************************
     * \brief   Type of nodes information of a hyperedge.
     *
     * A \c NodeDescT::value_type comprising the information about faces of a hyperedge.
     **********************************************************************************************/
    DataT& data;
    /*!*********************************************************************************************
     * \brief   Construct the \c struct that contains geometrical and topological information on a
     *          hyperedge.
     *
     * Construct a \c struct \c HyperEdge which is the value type of a \c HDGHyperGraph and contains
     * topolopgical and geometrical information on a hyperedge.
     *
     * \param   topo    Topological information of a hyperedge.
     * \param   geom    Geometrical information of a hyperedge.
     **********************************************************************************************/
    hyEdge
    (
      const typename TopoT::value_type& topo,
      const typename GeomT::value_type& geom,
      const typename NodeT::value_type& node,
      DataT& dat
    )
    : topology(topo), geometry(geom), node_descriptor(node), data(dat) { }
  } value_type; // end of typedef struct hyEdge
  
  /*!***********************************************************************************************
   * \brief   Iterator for \c struct \c hyEdge returned by \c operator[].
   *
   * Iterator that allows to go through the hyperedges of a hypergraph forwards and backwards. This
   * iterator fulfills the preconditions to allow the use of \c std::for_each on the set of
   * hyperedges that are contained in the \c HDGHyperGraph.
   ************************************************************************************************/
  class iterator
  {
    private:
      /*!*******************************************************************************************
       * \brief   Reference to the \c HDGHyperGraph of the iterator.
       *
       * The \c hyEdge is characterized via its respective \c HDGHypergraph (of which the
       * reference is saved) and its index who need to be members of the \c iterator.
       ********************************************************************************************/
      HDGHyperGraph& hyGraph_;
      /*!*******************************************************************************************
       * \brief   Index of the \c hyEdge of the iterator.
       *
       * The \c hyEdge is characterized via its respective \c HDGHypergraph (of which the
       * reference is saved) and its index who need to be members of the \c iterator.
       ********************************************************************************************/
      hyEdge_index_t index_;
    public:
/*    
      using iterator_category = std::random_access_iterator_tag;
      using value_type = HDGHyperGraph::value_type;
      using difference_type = int;
      using pointer = HDGHyperGraph::value_type*;
      using reference = HDGHyperGraph::value_type&;
*/      
      /*!*******************************************************************************************
       * \brief   Construct an iterator from an \c HDGHyperGraph and an index.
       * 
       * Construct \c HDGHyperGraph::iterator by passing over an \c HDGHyperGraph object and the
       * index the iterator is supposed to dot at.
       *
       * \param   hyGraph    The \c HDGHyperGraph, the iterator refers to.
       * \param   index         Index of the object, the iterator dots at.
       ********************************************************************************************/
      iterator(HDGHyperGraph& hyGraph, const hyEdge_index_t index)
      : hyGraph_(hyGraph), index_(index) { }
      /*!*******************************************************************************************
       * \brief   Copy--construct an iterator from another iterator.
       * 
       * Construct \c HDGHyperGraph::iterator as copy of another one.
       *
       * \param   other         Other \c iterator which is copied.
       ********************************************************************************************/
      iterator(const iterator& other)
      : hyGraph_(other.hyGraph_), index_(other.index_) { }
      /*!*******************************************************************************************
       * \brief   Copy--assign an iterator from another iterator.
       * 
       * Asign a given \c HDGHyperGraph::iterator to be a copy of another one
       *
       * \param   other         Other \c iterator which is copied.
       ********************************************************************************************/
      iterator& operator=(const iterator& other) = default;
      /*!*******************************************************************************************
       * \brief   Increment iterator and return incremented iterator.
       *
       * This function incements the iterator and returns the incremented iterator. Thus, no new
       * iterator needs to be constructed and only a reference needs to be returned. This makes the
       * function more performant compared to \c iterator \c operator++(int).
       * It is executed using \c ++iterator and \b not \c iterator++.
       *
       * \retval  incemented    The incremented iterator.
       ********************************************************************************************/
      iterator& operator++() { ++index_; return *this; }
      /*!*******************************************************************************************
       * \brief   Decrement iterator and return incremented iterator.
       *
       * This function decements the iterator and returns the decremented iterator. Thus, no new
       * iterator needs to be constructed and only a reference needs to be returned. This makes the
       * function more performant compared to \c iterator \c operator--(int).
       * It is executed using \c --iterator and \b not \c iterator--.
       *
       * \retval  decemented    The decremented iterator.
       ********************************************************************************************/
      iterator& operator--() { --index_; return *this; }
      /*!*******************************************************************************************
       * \brief   Increment iterator and return old iterator.
       *
       * This function incements the iterator and returns the old iterator. Thus, a new iterator
       * needs to be constructed and only a reference needs to be returned. This makes the function
       * less performant compared to \c iterator \c operator++().
       * It is executed using \c iterator++ and \b not \c ++iterator.
       *
       * \retval  incemented    The incremented iterator.
       ********************************************************************************************/
      iterator operator++(int) { return iterator(hyGraph_, index_++); }
      /*!*******************************************************************************************
       * \brief   Decrement iterator and return old iterator.
       *
       * This function decements the iterator and returns the old iterator. Thus, a new iterator
       * needs to be constructed and only a reference needs to be returned. This makes the function
       * less performant compared to \c iterator \c operator--().
       * It is executed using \c iterator-- and \b not \c --iterator.
       *
       * \retval  decemented    The decremented iterator.
       ********************************************************************************************/
      iterator operator--(int) { return iterator(hyGraph_, index_--); }
      /*!*******************************************************************************************
       * \brief   Dereference \c iterator to \c hyEdge.
       *
       * This function dereferences the iterator and returns the \c hyEdge this iterator dots at.
       *
       * \retval  hyEdge     The hyperedge described by the iterator.
       ********************************************************************************************/
      HDGHyperGraph::value_type operator*() { return hyGraph_[index_]; }
      /*!*******************************************************************************************
       * \brief   Check for equality with another iterator.
       *
       * This function checks whether the current iterator is equal to anoother \c iterator. In this
       * context equal means that they refer to the same \c HDGHyperGraph and have the same index.
       *
       * \param   other         \c iterator which is checked to be equal.
       * \retval  is_equal      \c boolean which is true if both iterators are equal and false
       *                        otherwise.
       ********************************************************************************************/
      bool operator==(const iterator& other)
      {
        return index_ == other.index_ && std::addressof(hyGraph_) == std::addressof(other.hyGraph_);
      }
      /*!*******************************************************************************************
       * \brief   Check for unequality with another iterator.
       *
       * This function checks whether the current iterator is equal to anoother \c iterator. In this
       * context unequal means that they do not refer to the same \c HDGHyperGraph or do not have
       * the same index.
       *
       * \param   other         \c iterator which is checked to be unequal.
       * \retval  is_equal      \c boolean which is false if both iterators are equal and true
       *                        otherwise.
       ********************************************************************************************/
      bool operator!=(const iterator& other)
      { 
        return index_ != other.index_ || std::addressof(hyGraph_) != std::addressof(other.hyGraph_);
      }
  }; // end of class iterator
  
  public:
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of a hyperedge.
     *
     * \retval  hyEdge_dim            The dimension of a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return TopoT::hyEdge_dim(); }
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of the space.
     *
     * \retval  space_dim             The dimension of the space.
     **********************************************************************************************/
    static constexpr unsigned int space_dim() { return TopoT::space_dim(); }
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the amount of dofs per node.
     *
     * \retval  n_dofs_per_nodeT       The amount of degrees of freedom per node.
     **********************************************************************************************/
    static constexpr unsigned int n_dofs_per_node() { return n_dofs_per_nodeT; }
  private:
    /*!*********************************************************************************************
     * \brief   Topology of the hypergraph.
     *
     * This object contains the topology of the hypergraph, i.e., it encodes which hyperedges
     * connect which hypernodes.
     **********************************************************************************************/
    std::shared_ptr<const TopoT> hyGraph_topology_;
    /*!*********************************************************************************************
     * \brief   Geometry of the hypergraph.
     *
     * This object contains the geometry of the hypergraph, i.e., it encodes the geometry of the
     * various hyperedges.
     **********************************************************************************************/
    std::shared_ptr<const GeomT> hyGraph_geometry_;
    /*!*********************************************************************************************
     * \brief   Node descriptions of the hypergraph.
     *
     * This object contains the nodal information of the hypergraph.
     **********************************************************************************************/
    std::shared_ptr<const NodeT> hyGraph_node_des_;
    /*!*********************************************************************************************
     * \brief   Hypernode factory administrating the access to degrees of freedom.
     *
     * A \c HyperNodeFactory allowing to connect nodes of the hypergraph to degrees of freedom which
     * are located in some \c std::vector. Note that this HyperNodeFactory has the same index type
     * for the hypernodes as this class for the hyperedges.
     **********************************************************************************************/
    const HyperNodeFactory<n_dofs_per_nodeT, hyEdge_index_t> hyNode_factory_;
    /*!*********************************************************************************************
     * \brief   Hypernode factory administrating the access to degrees of freedom.
     *
     * A \c HyperNodeFactory allowing to connect nodes of the hypergraph to degrees of freedom which
     * are located in some \c std::vector. Note that this HyperNodeFactory has the same index type
     * for the hypernodes as this class for the hyperedges.
     **********************************************************************************************/
    HyDataContainer<DataT> hyData_cont_;
  public:
    /*!*********************************************************************************************
     * \brief   Construct \c HDGHyperGraph from \c constructor_value_type.
     *
     * \todo Where does an HDGHyperGraph constructed like this get its geometry from?
     *
     * \todo Do not put "\c" in front of class names
     *
     * This is one of two standard ways of constructing a hypergraph.
     * That is, a hypergraph is constructed by providing the necessary data in form of the
     * respective \c constructor_value_type.
     *
     * \param   constructor           Information needed to deduce topological and geometrical data
     *                                to construct a \c HDGHyperGraph.
     **********************************************************************************************/
    HDGHyperGraph ( const typename TopoT::constructor_value_type& construct_topo )
    : hyGraph_topology_ ( std::make_shared<const TopoT>( construct_topo )),
      hyGraph_geometry_ ( std::make_shared<const GeomT>( *hyGraph_topology_ )),
      hyGraph_node_des_ ( std::make_shared<const NodeT>( *hyGraph_topology_ )),
      hyNode_factory_   ( hyGraph_topology_->n_hyNodes()),
      hyData_cont_      ( hyGraph_topology_->n_hyEdges())
    {
      static_assert( TopoT::hyEdge_dim() == GeomT::hyEdge_dim() ,
                     "The dimension of topology and geometry should be equal!" );
      hy_assert( hyNode_factory_.n_hyNodes() == hyGraph_topology_->n_hyNodes() ,
                 "The amount of hypernodes known to the hypernode factory is " <<
                 hyNode_factory_.n_hyNodes() << ", which is not equal to the amount that the"
                 << " hypergraph assumes, i.e., " << hyGraph_topology_->n_hyNodes() << "." );
      hy_assert( hyNode_factory_.n_hyNodes() >= 2 ,
                 "A hypergraph is assumed to consist of at least two hypernodes. This graph only "
                 << "consists of " << hyNode_factory_.n_hyNodes() << " hypernodes." );
      hy_assert( hyGraph_topology_->n_hyEdges() > 0 ,
                 "A hypergraph is supposed to consist of at least one hyperedge. This graph "
                 << "consists of " << hyGraph_topology_->n_hyEdges() << " hyperedges." );
    }
    /*!*********************************************************************************************
     * \brief   Construct \c HDGHyperGraph from \c constructor_value_type.
     *
     * This is one of two standard ways of constructing a hypergraph.
     * That is, a hypergraph is constructed by providing the necessary data to construct its
     * topology and its geometry seperately in form of the respective \c constructor_value_type
     * (plural, two).
     *
     * \param   construct_topo        Information needed to deduce topological data.
     * \param   construct_geom        Information needed to deduce geometrical data.
     **********************************************************************************************/
    HDGHyperGraph( const typename TopoT::constructor_value_type& construct_topo,
                   const typename GeomT::constructor_value_type& construct_geom )
    : hyGraph_topology_ ( std::make_shared<const TopoT>(construct_topo)),
      hyGraph_geometry_ ( std::make_shared<const GeomT>(construct_geom)),
      hyGraph_node_des_ ( std::make_shared<const NodeT>(*hyGraph_topology_)),
      hyNode_factory_   ( hyGraph_topology_->n_hyNodes()),
      hyData_cont_      ( hyGraph_topology_->n_hyEdges())
    {
      static_assert( TopoT::hyEdge_dim() == GeomT::hyEdge_dim() ,
                     "The dimension of topology and geometry should be equal!" );
      hy_assert( hyNode_factory_.n_hyNodes() == hyGraph_topology_->n_hyNodes() ,
                 "The amount of hypernodes known to the hypernode factory is " <<
                 hyNode_factory_.n_hyNodes() << ", which is not equal to the amount that the"
                 << " hypergraph assumes, i.e., " << hyGraph_topology_->n_hyNodes() << "." );
      hy_assert( hyNode_factory_.n_hyNodes() >= 2 ,
                 "A hypergraph is assumed to consist of at least two hypernodes. This graph only "
                 << "consists of " << hyNode_factory_.n_hyNodes() << " hypernodes." );
      hy_assert( hyGraph_topology_->n_hyEdges() > 0 ,
                 "A hypergraph is supposed to consist of at least one hyperedge. This graph "
                 << "consists of " << hyGraph_topology_->n_hyEdges() << " hyperedges." );
    }
    /*!*********************************************************************************************
     * \brief Construct HDGHyperGraph from existing topology and geometry
     *
     * This is one of two standard ways of constructing a hypergraph.
     * That is, a hypergraph is constructed by providing the necessary data to construct its
     * topology and its geometry seperately in form of the respective \c constructor_value_type
     * (plural, two).
     *
     * \param   construct_topo        Information needed to deduce topological data.
     * \param   construct_geom        Information needed to deduce geometrical data.
     **********************************************************************************************/
    HDGHyperGraph
    ( std::shared_ptr<const TopoT> topo,
      std::shared_ptr<const GeomT> geom,
      std::shared_ptr<const NodeT> node
    )
    : hyGraph_topology_ ( topo ), hyGraph_geometry_ ( geom ), hyGraph_node_des_ ( node ),
      hyNode_factory_   ( hyGraph_topology_->n_hyNodes() ),
      hyData_cont_      ( hyGraph_topology_->n_hyEdges())
    {
      static_assert( TopoT::hyEdge_dim() == GeomT::hyEdge_dim() ,
                     "The dimension of topology and geometry should be equal!" );
      hy_assert( hyNode_factory_.n_hyNodes() == hyGraph_topology_->n_hyNodes() ,
                 "The amount of hypernodes known to the hypernode factory is " <<
                 hyNode_factory_.n_hyNodes() << ", which is not equal to the amount that the"
                 << " hypergraph assumes, i.e., " << hyGraph_topology_->n_hyNodes() << "." );
      hy_assert( hyNode_factory_.n_hyNodes() >= 2 ,
                 "A hypergraph is assumed to consist of at least two hypernodes. This graph only "
                 << "consists of " << hyNode_factory_.n_hyNodes() << " hypernodes." );
      hy_assert( hyGraph_topology_->n_hyEdges() > 0 ,
                 "A hypergraph is supposed to consist of at least one hyperedge. This graph "
                 << "consists of " << hyGraph_topology_->n_hyEdges() << " hyperedges." );
    }
    /*!*********************************************************************************************
     * \brief   Subscript operator of a \c HDGHyperGraph.
     *
     * \todo  Here, we repeatedly return a large object. This is done since the object could be
     *        locally created in regular topologies/geometries! Return shared-pointer?
     *
     * The subscript operator takes an index referring to an hyperedge and returns the respective
     * \c hyEdge containing its topological and geometrical information. Thus, this operator can
     * be bypassed by using the functions \c hyEdge_topology (only returning the topological
     * data) and \c hyEdge_geometry (ony returning the geometrical data).
     *
     * \param   index                 Index of the \c hyEdge to be returned.
     * \retval  hyEdge             The \c hyEdge of the given index.
     **********************************************************************************************/
    value_type operator[] (const hyEdge_index_t index)
    { return value_type(hyEdge_topology(index), hyEdge_geometry(index), hyNode_descriptor(index),
                        hyEdge_data(index)); }
    /*!*********************************************************************************************
     * \brief   Return iterator to first \c hyEdge of \c HDGHyperGraph.
     *
     * This function returns an \c HDGHyperGraph::iterator that refers to the first \c hyEdge of
     * the hypergraph (index = 0). Thus, it can be used to mark the starting point in \c for_each
     * loops.
     *
     * \retval  hyEdge             Iterator referring to first \c hyEdge.
     **********************************************************************************************/
    typename HDGHyperGraph<n_dofs_per_nodeT, TopoT, GeomT,NodeT,DataT, hyEdge_index_t >::iterator
    begin()
    { 
      return HDGHyperGraph< n_dofs_per_nodeT, TopoT, GeomT, NodeT,DataT, hyEdge_index_t >::iterator
               (*this, 0);
    }
    /*!*********************************************************************************************
     * \brief   Return iterator to the end of \c hyEdge list.
     *
     * This function returns an \c HDGHyperGraph::iterator that refers to the position of an (non-
     * existing) \c hyEdge of the hypergraph (index = n_hyEdges), i.e., the position  directly
     * after the last valid entry of the \c HDGHyperGraph. Thus, it can be used to mark the ending
     * point in \c for_each loops.
     *
     * \retval  hyEdge             Iterator referring to position behind last \c hyEdge.
     **********************************************************************************************/
    typename HDGHyperGraph<n_dofs_per_nodeT, TopoT, GeomT, NodeT, DataT, hyEdge_index_t >::iterator
    end()
    { 
      return HDGHyperGraph<n_dofs_per_nodeT,TopoT,GeomT,NodeT,DataT,hyEdge_index_t>::iterator
               (*this, n_hyEdges());
    }
    /*!*********************************************************************************************
     * \brief   Return const reference to HyperNodeFactory.
     *
     * \todo    Why is this public? Why not get_hypernode()?
     *          -> Because a hypernode would have several functions and we decided not to introduce
     *          a hypernode class, but to only have a hypernode factory covering all those aspects.
     *          What we could do is to repeat all functions of the hypernode_factory in this class.
     * 
     * This function returns an \c HyperNodeFactory handling the access to the degrees of freedom
     * encoded in some \c std::vector.
     *
     * \retval  hypernode_factory     The \c HyperNodeFactory belonging the hypergraph.
     **********************************************************************************************/
    const HyperNodeFactory<n_dofs_per_nodeT,hyEdge_index_t>& hyNode_factory() const
    { return hyNode_factory_; }
    /*!*********************************************************************************************
     * \brief   Topological information of prescribed hyperedge.
     *
     * \todo  Here, we repeatedly return a large object. This is done since the object could be
     *        locally created in regular topologies/geometries! Return shared-pointer?
     *        -> I do not see the advanatge of that. Do we not just return a pointer to a large
     *        object then?
     *
     * Return the topological information of a specific hyperedge identified via its index. This
     * function can be used to bypass the subscript operator which returns topological and geometric
     * information about a hyperedge of given index.
     *
     * \param   index                 Index of the hyperedge to be returned.
     * \retval  hyEdge_topology    Topological information about hyperedge.
     **********************************************************************************************/
    const typename TopoT::value_type hyEdge_topology(const hyEdge_index_t index) const
    { return hyGraph_topology_->operator[](index); }
    /*!*********************************************************************************************
     * \brief   Geometrical information of prescribed hyperedge.
     *
     * \todo  Here, we repeatedly return a large object. This is done since the object could be
     *        locally created in regular topologies/geometries! Return shared-pointer?
     *        -> I do not see the advanatge of that. Do we not just return a pointer to a large
     *        object then?
     *
     * Return the geometrical information of a specific hyperedge identified via its index. This
     * function can be used to bypass the subscript operator which returns topological and geometric
     * information about a hyperedge of given index.
     *
     * \param   index                 Index of the hyperedge to be returned.
     * \retval  hyEdge_geometry    Geometrical information about hyperedge.
     **********************************************************************************************/
    const typename GeomT::value_type hyEdge_geometry(const hyEdge_index_t index) const
    { return hyGraph_geometry_->operator[](index); }
    /*!*********************************************************************************************
     * \brief   Nodal information of prescribed hyperedge.
     *
     * \todo  Here, we repeatedly return a large object. This is done since the object could be
     *        locally created in regular topologies/geometries! Return shared-pointer?
     *        -> I do not see the advanatge of that. Do we not just return a pointer to a large
     *        object then?
     *
     * \todo  All
     *
     * \param   index                 Index of the hyperedge to be returned.
     * \retval  hyEdge_geometry    Geometrical information about hyperedge.
     **********************************************************************************************/
    const typename NodeT::value_type hyNode_descriptor(const hyEdge_index_t index) const
    { return hyGraph_node_des_->operator[](index); }
    /*!*********************************************************************************************
     * \brief   Nodal information of prescribed hyperedge.
     *
     * \todo  Here, we repeatedly return a large object. This is done since the object could be
     *        locally created in regular topologies/geometries! Return shared-pointer?
     *        -> I do not see the advanatge of that. Do we not just return a pointer to a large
     *        object then?
     *
     * \todo  All
     *
     * \param   index                 Index of the hyperedge to be returned.
     * \retval  hyEdge_geometry    Geometrical information about hyperedge.
     **********************************************************************************************/
    DataT& hyEdge_data(const hyEdge_index_t index)
    { return hyData_cont_.operator[](index); }
    /*!*********************************************************************************************
     * \brief   Returns the number of hyperedges making up the hypergraph.
     *
     * \retval  n_hyEdges          The total amount of hyperedges of a hypergraph.
     **********************************************************************************************/
    const hyEdge_index_t n_hyEdges() const  { return hyGraph_topology_->n_hyEdges(); }
    /*!*********************************************************************************************
     * \brief   Returns the number of hypernodes making up the hypergraph.
     *
     * \retval  n_hypernodes          The total amount of hypernodes of a hypergraph.
     **********************************************************************************************/
    const hyEdge_index_t n_hyNodes() const  { return hyNode_factory_.n_hyNodes(); }
    /*!*********************************************************************************************
     * \brief   Returns the total amount of degrees of freedom in the considered hypergraph.
     * 
     * \retval  n_global_dofs         The total amount of degreees of freedom in the considered
     *                                hypergraph.
     **********************************************************************************************/
    template < typename dof_index_t = unsigned int >
    const dof_index_t n_global_dofs() const  { return hyNode_factory_.n_global_dofs(); }
}; // end of class HDGHyperGraph
