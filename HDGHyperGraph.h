#ifndef HYPERGRAPH_H
#define HYPERGRAPH_H

#include "TypeDefs.h"
#include "HyperNodeFactory.h"
#include "HyperGraph_Topology.h"
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
 * @authors   Guido Kanschat, University of Heidelberg, 2019.
 * @authors   Andreas Rupp, University of Heidelberg, 2019.
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
  } value_type;
  
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
  };
  
  private:
    /*!*********************************************************************************************
     * @brief   Hypernode factory administrating the access to degrees of freedom.
     *
     * A @c HyperNodeFactory allowing to connect nodes of the hypergraph to degrees of freedom which
     * are located in some @c std::vector.
     **********************************************************************************************/
    const HyperNodeFactory<n_dofs_per_node> hypernode_factory_;
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
  public:
    HDGHyperGraph(const TopoT& hyperedge_getter);
    const value_type operator[] (const hyperedge_index_type index) const;
    
    typename HDGHyperGraph<n_dofs_per_node, TopoT, GeomT >::iterator begin() const;
    typename HDGHyperGraph<n_dofs_per_node, TopoT, GeomT >::iterator end() const;
    
    // Why is this public? Why not get_hypernode()?
    const HyperNodeFactory<n_dofs_per_node> hypernode_factory() const; // AR: No reference for performance?!
    const typename TopoT::value_type hyperedge_topology(const hyperedge_index_type index) const;
    const typename GeomT::value_type hyperedge_geometry(const hyperedge_index_type index) const;
    
    const hyperedge_index_type num_of_hyperedges() const;
    const hypernode_index_type num_of_hypernodes() const;
    const dof_index_type num_of_global_dofs() const;
    
    static constexpr unsigned int hyperedge_dimension() { return TopoT::hyperedge_dimension(); };
    static constexpr unsigned int space_dimension() { return TopoT::space_dimension(); };
    static constexpr unsigned int n_dof_per_node() { return n_dofs_per_node; }
};

#endif
