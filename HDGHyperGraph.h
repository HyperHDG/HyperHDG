#ifndef HDGHYPERGRAPH_H
#define HDGHYPERGRAPH_H

#include "TypeDefs.h"
#include "HyperNodeFactory.h"
#include "Topo_Cubic.h"
#include "HyperGraph_Geometry.h"

/*!*************************************************************************************************
 * @brief   The class template uniting topology and geometry of a hypergraph with the topology of
 *          the skeleton space of the HDG method.
 *
 * \todo Is the name ok? It could be HDGHyperGraph and include the HDG loop given a local solver
 *       -> This has been implemented defining an iterator and using the std::for each in the
 *       DiffusionProblem implementation. I hope that this is even more general, since for_each
 *       is pretty general. What do you think?
 * 
 * The main class representing a hypergraph. It uses a class @c Topology to represent the collection
 * of nodes and edges as well as a class @c Geometry presenting the physical coordinates of the
 * edges. It behaves like a random access container of hyperedges and has additional access to its
 * nodes.
 *
 * In our abstraction, nodes only carry degrees of freedom. Thus, they can be obtained from one
 * object @c HyperNodeFactory for any graph. Their location, if such a notion is reasonable, must be
 * determined by that of the boundaries of an edge. The meaning of their degrees of freedom is
 * decided by the local solvers of the HDG method applied. The @c Geometry class may use degrees of
 * freedom of the nodes as well.
 *
 * @tparam  n_dofs_per_node The number of degrees of freedom of a single hypernode which is assumed
 *                          to be the same for all hypernodes.
 * @tparam  TopoT           Class that contains the topology of the hypergraph. This class is needs
 *                          to provide a getter function to the topological information of a
 *                          hyperedge of given index and can be arbitrarily implemented.
 * @tparam  GeomT           Class that contains the topology of the hypergraph. This class is needs
 *                          to provide a getter function to the topological information of a
 *                          hyperedge of given index and can be arbitrarily implemented.
 *
 * @authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * @authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
class HDGHyperGraph
{
  
  /*!***********************************************************************************************
   * @brief   The type for a hyperedge returned by @c operator[].
   *
   * This @c typedef @c struct is returned by the @c operator[] of an @c HDGHyperGraph. It contains
   * topological and geometrical information about a single hyperedge. It is therefore defined as
   * the @c value_type of class @c HDGHyperGraph.
   ************************************************************************************************/
  typedef struct HyperEdge
  {
    /*!*********************************************************************************************
     * @brief   Topological information of a hyperedge.
     *
     * A @c TopoT::value_type comprising the topological information of a hyperedge.
     **********************************************************************************************/
    typename TopoT::value_type topology;
    /*!*********************************************************************************************
     * @brief   Geometrical information of a hyperedge.
     *
     * A @c TopoT::value_type comprising the geometrical information of a hyperedge.
     **********************************************************************************************/
    typename GeomT::value_type geometry;
    /*!*********************************************************************************************
     * @brief   Construct the @c struct that contains geometrical and topological information on a
     *          hyperedge.
     *
     * Construct a @c struct @c HyperEdge which is the value type of a @c HDGHyperGraph and contains
     * topolopgical and geometrical information on a hyperedge.
     *
     * @param   topo    Topological information of a hyperedge.
     * @param   geom    Geometrical information of a hyperedge.
     **********************************************************************************************/
    HyperEdge(const typename TopoT::value_type& topo, const typename GeomT::value_type& geom)
      : topology(topo), geometry(geom) { };
  } value_type; // end of typedef struct HyperEdge
  
  /*!***********************************************************************************************
   * @brief   The type needed to construct an @c HDGHyperGraph.
   *
   * This @c typedef @c struct is needed to defaultly construct an @c HDGHyperGraph. It contains
   * topological and geometrical construction information about the hypergraph. It is therefore
   * defined as the @c constructor_value_type of class @c HDGHyperGraph.
   ************************************************************************************************/
  typedef struct HyperGraphConstructor
  {
    /*!*********************************************************************************************
     * @brief   Constructor arguments for topological information of a hypergraph.
     *
     * A @c TopoT::constructor_value_type comprising the information to construct topological
     * information of a hypergraph.
     **********************************************************************************************/
    typename TopoT::constructor_value_type topology;
    /*!*********************************************************************************************
     * @brief   Constructor arguments for geometrical information of a hypergraph.
     *
     * A @c TopoT::constructor_value_type comprising the information to construct geometrical
     * information of a hypergraph.
     **********************************************************************************************/
    typename GeomT::constructor_value_type geometry;
    /*!*********************************************************************************************
     * @brief   Construct the @c struct that contains geometrical and topological information to
     *          construct a hypergraph.
     *
     * Construct a @c struct @c HyperGraphConstructor which is needed to construct a
     * @c HDGHyperGraph and contains topolopgical and geometrical information to do so.
     *
     * @param   topo    Topological information to construct hypergraph.
     * @param   geom    Geometrical information to construct hypergraph.
     **********************************************************************************************/
    HyperGraphConstructor(const typename TopoT::constructor_value_type& topo,
                          const typename GeomT::constructor_value_type& geom)
      : topology(topo), geometry(geom) { };
  } constructor_value_type; // end of typedef struct HyperGraphConstructor
  
  /*!***********************************************************************************************
   * @brief   Iterator for @c struct @c HyperEdge returned by @c operator[].
   *
   * Iterator that allows to go through the hyperedges of a hypergraph forwards and backwards. This
   * iterator fulfills the preconditions to allow the use of @c std::for_each on the set of
   * hyperedges that are contained in the @c HDGHyperGraph.
   ************************************************************************************************/
  class iterator
  {
    private:
      /*!*******************************************************************************************
       * @brief   Reference to the @c HDGHyperGraph of the iterator.
       *
       * The @c HyperEdge is characterized via its respective @c HDGHypergraph (of which the
       * reference is saved) and its index who need to be members of the @c iterator.
       ********************************************************************************************/
      const HDGHyperGraph& hypergraph_;
      /*!*******************************************************************************************
       * @brief   Index of the @c HyperEdge of the iterator.
       *
       * The @c HyperEdge is characterized via its respective @c HDGHypergraph (of which the
       * reference is saved) and its index who need to be members of the @c iterator.
       ********************************************************************************************/
      hyperedge_index_type index_;
    public:
/*    
      using iterator_category = std::random_access_iterator_tag;
      using value_type = HDGHyperGraph::value_type;
      using difference_type = int;
      using pointer = HDGHyperGraph::value_type*;
      using reference = HDGHyperGraph::value_type&;
*/      
      /*!*******************************************************************************************
       * @brief   Construct an iterator from an @c HDGHyperGraph and an index.
       * 
       * Construct @c HDGHyperGraph::iterator by passing over an @c HDGHyperGraph object and the
       * index the iterator is supposed to dot at.
       *
       * @param   hypergraph    The @c HDGHyperGraph, the iterator refers to.
       * @param   index         Index of the object, the iterator dots at.
       ********************************************************************************************/
      iterator(const HDGHyperGraph& hypergraph, const hyperedge_index_type index)
        : hypergraph_(hypergraph), index_(index) {};
      /*!*******************************************************************************************
       * @brief   Copy--construct an iterator from another iterator.
       * 
       * Construct @c HDGHyperGraph::iterator as copy of another one.
       *
       * @param   other         Other @c iterator which is copied.
       ********************************************************************************************/
      iterator(const iterator& other)
        : hypergraph_(other.hypergraph_), index_(other.index_) {};
      /*!*******************************************************************************************
       * @brief   Copy--assign an iterator from another iterator.
       * 
       * Asign a given @c HDGHyperGraph::iterator to be a copy of another one
       *
       * @param   other         Other @c iterator which is copied.
       ********************************************************************************************/
      iterator& operator=(const iterator& other) = default;
      /*!*******************************************************************************************
       * @brief   Increment iterator and return incremented iterator.
       *
       * This function incements the iterator and returns the incremented iterator. Thus, no new
       * iterator needs to be constructed and only a reference needs to be returned. This makes the
       * function more performant compared to @c iterator @c operator++(int).
       * It is executed using @c ++iterator and @b not @c iterator++.
       *
       * @retval  incemented    The incremented iterator.
       ********************************************************************************************/
      iterator& operator++() { ++index_; return *this; };
      /*!*******************************************************************************************
       * @brief   Decrement iterator and return incremented iterator.
       *
       * This function decements the iterator and returns the decremented iterator. Thus, no new
       * iterator needs to be constructed and only a reference needs to be returned. This makes the
       * function more performant compared to @c iterator @c operator--(int).
       * It is executed using @c --iterator and @b not @c iterator--.
       *
       * @retval  decemented    The decremented iterator.
       ********************************************************************************************/
      iterator& operator--() { --index_; return *this; };
      /*!*******************************************************************************************
       * @brief   Increment iterator and return old iterator.
       *
       * This function incements the iterator and returns the old iterator. Thus, a new iterator
       * needs to be constructed and only a reference needs to be returned. This makes the function
       * less performant compared to @c iterator @c operator++().
       * It is executed using @c iterator++ and @b not @c ++iterator.
       *
       * @retval  incemented    The incremented iterator.
       ********************************************************************************************/
      iterator operator++(int) { return iterator(hypergraph_, index_++); };
      /*!*******************************************************************************************
       * @brief   Decrement iterator and return old iterator.
       *
       * This function decements the iterator and returns the old iterator. Thus, a new iterator
       * needs to be constructed and only a reference needs to be returned. This makes the function
       * less performant compared to @c iterator @c operator--().
       * It is executed using @c iterator-- and @b not @c --iterator.
       *
       * @retval  decemented    The decremented iterator.
       ********************************************************************************************/
      iterator operator--(int) { return iterator(hypergraph_, index_--); };
      /*!*******************************************************************************************
       * @brief   Dereference @c iterator to @c HyperEdge.
       *
       * This function dereferences the iterator and returns the @c HyperEdge this iterator dots at.
       *
       * @retval  hyperedge     The hyperedge described by the iterator.
       ********************************************************************************************/
      HDGHyperGraph::value_type operator*() { return hypergraph_[index_]; };
      /*!*******************************************************************************************
       * @brief   Check for equality with another iterator.
       *
       * This function checks whether the current iterator is equal to anoother @c iterator. In this
       * context equal means that they refer to the same @c HDGHyperGraph and have the same index.
       *
       * @param   other         @c iterator which is checked to be equal.
       * @retval  is_equal      @c boolean which is true if both iterators are equal and false
       *                        otherwise.
       ********************************************************************************************/
      bool operator==(const iterator& other)
      {
        return index_ == other.index_ 
               && std::addressof(hypergraph_) == std::addressof(other.hypergraph_);
      };
      /*!*******************************************************************************************
       * @brief   Check for unequality with another iterator.
       *
       * This function checks whether the current iterator is equal to anoother @c iterator. In this
       * context unequal means that they do not refer to the same @c HDGHyperGraph or do not have
       * the same index.
       *
       * @param   other         @c iterator which is checked to be unequal.
       * @retval  is_equal      @c boolean which is false if both iterators are equal and true
       *                        otherwise.
       ********************************************************************************************/
      bool operator!=(const iterator& other)
      { return index_ != other.index_ 
               || std::addressof(hypergraph_) != std::addressof(other.hypergraph_);
      };
  }; // end of class iterator
  
  private:
    /*!*********************************************************************************************
     * @brief   Topology of the hypergraph.
     *
     * This object contains the topology of the hypergraph, i.e., it encodes which hyperedges
     * connect which hypernodes.
     **********************************************************************************************/
    const TopoT hypergraph_topology_;
    /*!*********************************************************************************************
     * @brief   Geometry of the hypergraph.
     *
     * This object contains the geometry of the hypergraph, i.e., it encodes the geometry of the
     * various hyperedges.
     **********************************************************************************************/
    const GeomT hypergraph_geometry_;
    /*!*********************************************************************************************
     * @brief   Hypernode factory administrating the access to degrees of freedom.
     *
     * A @c HyperNodeFactory allowing to connect nodes of the hypergraph to degrees of freedom which
     * are located in some @c std::vector.
     **********************************************************************************************/
    const HyperNodeFactory<n_dofs_per_node> hypernode_factory_;
  public:
    /*!*********************************************************************************************
     * @brief   Construct @c HDGHyperGraph from @c constructor_value_type.
     *
     * This is one of two standard ways of constructing a hypergraph.
     * That is, a hypergraph is constructed by providing the necessary data in form of the
     * respective @c constructor_value_type.
     *
     * @param   constructor           Information needed to deduce topological and geometrical data
     *                                to construct a @c HDGHyperGraph.
     **********************************************************************************************/
    HDGHyperGraph(const constructor_value_type& construction_data);
    /*!*********************************************************************************************
     * @brief   Construct @c HDGHyperGraph from @c constructor_value_type.
     *
     * This is one of two standard ways of constructing a hypergraph.
     * That is, a hypergraph is constructed by providing the necessary data to construct its
     * topology and its geometry seperately in form of the respective @c constructor_value_type
     * (plural, two).
     *
     * @param   construct_topo        Information needed to deduce topological data.
     * @param   construct_geom        Information needed to deduce geometrical data.
     **********************************************************************************************/
    HDGHyperGraph(const typename TopoT::constructor_value_type& construct_topo,
                  const typename GeomT::constructor_value_type& construct_geom);
    /*!*********************************************************************************************
     * @brief   Subscript operator of a @c HDGHyperGraph.
     *
     * The subscript operator takes an index referring to an hyperedge and returns the respective
     * @c HyperEdge containing its topological and geometrical information. Thus, this operator can
     * be bypassed by using the functions @c hyperedge_topology (only returning the topological
     * data) and @c hyperedge_geometry (ony returning the geometrical data).
     *
     * @param   index                 Index of the @c HyperEdge to be returned.
     * @retval  hyperedge             The @c HyperEdge of the given index.
     **********************************************************************************************/
    const value_type operator[] (const hyperedge_index_type index) const;
    /*!*********************************************************************************************
     * @brief   Return iterator to first @c HyperEdge of @c HDGHyperGraph.
     *
     * This function returns an @c HDGHyperGraph::iterator that refers to the first @c HyperEdge of
     * the hypergraph (index = 0). Thus, it can be used to mark the starting point in @c for_each
     * loops.
     *
     * @retval  hyperedge             Iterator referring to first @c HyperEdge.
     **********************************************************************************************/
    typename HDGHyperGraph<n_dofs_per_node, TopoT, GeomT >::iterator begin() const;
    /*!*********************************************************************************************
     * @brief   Return iterator to the end of @c HyperEdge list.
     *
     * This function returns an @c HDGHyperGraph::iterator that refers to the position of an (non-
     * existing) @c HyperEdge of the hypergraph (index = n_hyperedges), i.e., the position  directly
     * after the last valid entry of the @c HDGHyperGraph. Thus, it can be used to mark the ending
     * point in @c for_each loops.
     *
     * @retval  hyperedge             Iterator referring to position behind last @c HyperEdge.
     **********************************************************************************************/
    typename HDGHyperGraph<n_dofs_per_node, TopoT, GeomT >::iterator end() const;
    /*!*********************************************************************************************
     * @brief   Return const reference to HyperNodeFactory.
     *
     * @todo    Why is this public? Why not get_hypernode()?
     *          -> Because a hypernode would have several functions and we decided not to introduce
     *          a hypernode class, but to only have a hypernode factory covering all those aspects.
     *          What we could do is to repeat all functions of the hypernode_factory in this class.
     * 
     * This function returns an @c HyperNodeFactory handling the access to the degrees of freedom
     * encoded in some @c std::vector.
     *
     * @retval  hypernode_factory     The @c HyperNodeFactory belonging the hypergraph.
     **********************************************************************************************/
    const HyperNodeFactory<n_dofs_per_node>& hypernode_factory() const;
    /*!*********************************************************************************************
     * @brief   Topological information of prescribed hyperedge.
     *
     * Return the topological information of a specific hyperedge identified via its index. This
     * function can be used to bypass the subscript operator which returns topological and geometric
     * information about a hyperedge of given index.
     *
     * @param   index                 Index of the hyperedge to be returned.
     * @retval  hyperedge_topology    Topological information about hyperedge.
     **********************************************************************************************/
    const typename TopoT::value_type hyperedge_topology(const hyperedge_index_type index) const;
    /*!*********************************************************************************************
     * @brief   Geometrical information of prescribed hyperedge.
     *
     * Return the geometrical information of a specific hyperedge identified via its index. This
     * function can be used to bypass the subscript operator which returns topological and geometric
     * information about a hyperedge of given index.
     *
     * @param   index                 Index of the hyperedge to be returned.
     * @retval  hyperedge_geometry    Geometrical information about hyperedge.
     **********************************************************************************************/
    const typename GeomT::value_type hyperedge_geometry(const hyperedge_index_type index) const;
    /*!*********************************************************************************************
     * @brief   Returns the number of hyperedges making up the hypergraph.
     *
     * @retval  n_hyperedges          The total amount of hyperedges of a hypergraph.
     **********************************************************************************************/
    const hyperedge_index_type n_hyperedges() const;
    /*!*********************************************************************************************
     * @brief   Returns the number of hypernodes making up the hypergraph.
     *
     * @retval  n_hypernodes          The total amount of hypernodes of a hypergraph.
     **********************************************************************************************/
    const hypernode_index_type n_hypernodes() const;
    /*!*********************************************************************************************
     * @brief   Returns the total amount of degrees of freedom in the considered hypergraph.
     * 
     * @retval  n_global_dofs         The total amount of degreees of freedom in the considered
     *                                hypergraph.
     **********************************************************************************************/
    const dof_index_type n_global_dofs() const;
    
    /*!*********************************************************************************************
     * @brief   Returns the template parameter representing the dimension of a hyperedge.
     *
     * @retval  hyperedge_dim         The dimension of a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyperedge_dimension() { return TopoT::hyperedge_dimension(); };
    /*!*********************************************************************************************
     * @brief   Returns the template parameter representing the dimension of the space.
     *
     * @retval  space_dim             The dimension of the space.
     **********************************************************************************************/
    static constexpr unsigned int space_dimension() { return TopoT::space_dimension(); };
    /*!*********************************************************************************************
     * @brief   Returns the template parameter representing the amount of dofs per node.
     *
     * @retval  n_dofs_per_node       The amount of degrees of freedom per node.
     **********************************************************************************************/
    static constexpr unsigned int n_dof_per_node() { return n_dofs_per_node; }
}; // end of class HDGHyperGraph

#endif // end of ifndef HYPERGRAPH_H
