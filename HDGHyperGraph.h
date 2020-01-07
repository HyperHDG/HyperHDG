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
   * @brief   The type for a hyperedge returned by \c operator[].
   *
   * This \c typedef \c struct is returned by the \c operator[] of an \c HDGHyperGraph. It contains
   * topological and geometrical information about a single hyperedge. It is therefore defined as
   * the \c value_type of class \c HDGHyperGraph.
   ************************************************************************************************/
  typedef struct HyperEdge
  {
    /*!*********************************************************************************************
     * @brief   Topological information of a hyperedge.
     *
     * A \c TopoT::value_type comprising the topological information of a hyperedge.
     **********************************************************************************************/
    typename TopoT::value_type topology;
    /*!*********************************************************************************************
     * @brief   Geometrical information of a hyperedge.
     *
     * A \c TopoT::value_type comprising the geometrical information of a hyperedge.
     **********************************************************************************************/
    typename GeomT::value_type geometry;
    /*!*********************************************************************************************
     * @brief   Construct the \c struct that contains geometrical and topological information on a
     *          hyperedge.
     *
     * Construct a \c struct \c HyperEdge which is the value type of a \c HDGHyperGraph and contains
     * topolopgical and geometrical information on a hyperedge.
     *
     * @param   topo    Topological information of a hyperedge.
     * @param   geom    Geometrical information of a hyperedge.
     **********************************************************************************************/
    HyperEdge(const typename TopoT::value_type& topo, const typename GeomT::value_type& geom)
      : topology(topo), geometry(geom) { };
  } value_type;
  
  /*!***********************************************************************************************
   * @brief   Iterator for \c struct \c HyperEdge returned by \c operator[].
   *
   * Iterator that allows to go through the hyperedges of a hypergraph forwards and backwards. This
   * iterator fulfills the preconditions to allow the use of \c std::for_each on the set of
   * hyperedges that are contained in the \c HDGHyperGraph.
   ************************************************************************************************/
  class iterator
  {
    private:
      /*!*******************************************************************************************
       * @brief   Reference to the \c HDGHyperGraph of the iterator.
       *
       * The \c HyperEdge is characterized via its respective \c HDGHypergraph (of which the
       * reference is saved) and its index who need to be members of the \c iterator.
       ********************************************************************************************/
      const HDGHyperGraph& hypergraph_;
      /*!*******************************************************************************************
       * @brief   Index of the \c HyperEdge of the iterator.
       *
       * The \c HyperEdge is characterized via its respective \c HDGHypergraph (of which the
       * reference is saved) and its index who need to be members of the \c iterator.
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
       * @brief   Construct an iterator from an \c HDGHyperGraph and an index.
       *
       * \todo Finish this explanation and copy constructor and copy assignment.
       * 
       * Constructs a hypergraph from a \c std::array containing the elementens per spatial dimension
       * which is given as input data. The array has the correct length (as ensured by the involved
       * template parametzer \c space_dim.
       *
       * @param   num_elements    A @c std::array containing number of elements per spatial dimension.
       ********************************************************************************************/
      iterator(const HDGHyperGraph& hypergraph, const hyperedge_index_type index)
        : hypergraph_(hypergraph), index_(index) {};
      iterator(const iterator& other)
        : hypergraph_(other.hypergraph_), index_(other.index_) {};
      iterator& operator=(const iterator& other) = default;
      
      iterator& operator++() { ++index_; return *this; };
      iterator& operator--() { --index_; return *this; };
      iterator operator++(int) { return iterator(hypergraph_, index_++); };
      iterator operator--(int) { return iterator(hypergraph_, index_--); };
      HDGHyperGraph::value_type operator*() { return hypergraph_[index_]; };
      bool operator==(const iterator& other) { return index_ == other.index_ && std::addressof(hypergraph_) == std::addressof(other.hypergraph_); };
      bool operator!=(const iterator& other) { return index_ != other.index_ || std::addressof(hypergraph_) != std::addressof(other.hypergraph_); };
  };
  
  private:
    const HyperNodeFactory<n_dofs_per_node> hypernode_factory_;
    const TopoT hypergraph_topology_;
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
