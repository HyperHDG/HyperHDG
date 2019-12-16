/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef HYPERGRAPH_H
#define HYPERGRAPH_H

#include "TypeDefs.h"
#include "HyperNodeFactory.h"
#include "Topology.h"
#include "Geometry.h"


/**
 * \brief The class template uniting topology and geometry of a
 * hypergraph with the topology of the skeleton space of the HDG
 * method.
 *
 * \todo Is the name ok? It could be HDGHyperGraph and include the HDG loop given a local solver
 *
 * The main class representing a hypergraph. It uses a class
 * `Topology` to represent the collection of nodes and edges as well
 * as a class `Geometry` presenting the physical coordinates of the
 * edges. It behaves like a random access container of hyperedges and
 * has additional access to its nodes.
 *
 * In our abstraction, nodes only carry degrees of freedom. Thus, they
 * can be obtained from one object HyperNodeFactory for any
 * graph. Their location, if such a notion is reasonable, must be
 * determined by that of the boundaries of an edge. The meaning of
 * their degrees of freedom is decided by the local solvers of the HDG
 * method applied. The Geometry class may use degrees of freedom of
 * the nodes as well.
 */
template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
class HDGHyperGraph
{
  /**
   * The type for a hyperedge returned by operator[].
   */
  typedef struct HyperEdge
  {
    typename TopoT::value_type topology;
    typename GeomT::value_type geometry;
    HyperEdge(const typename TopoT::value_type& topo, const typename GeomT::value_type& geom)
      : topology(topo), geometry(geom) {};
  } value_type;
/*  
  class iterator
  {
    iterator
    
  };
*/  
  private:
    const HyperNodeFactory<amount_of_local_dofs> hypernode_factory_;
    const TopoT hypergraph_topology_;
    const GeomT hypergraph_geometry_;
  
  public:
    HDGHyperGraph(const TopoT& hyperedge_getter);
    const value_type operator[] (const hyperedge_index_type index) const;
    
    // Why is this public? Why not get_hypernode()?
    const HyperNodeFactory<amount_of_local_dofs> hypernode_factory() const; // AR: No reference for performance?!
    const typename TopoT::value_type hyperedge_topology(const hyperedge_index_type index) const;
    const typename GeomT::value_type hyperedge_geometry(const hyperedge_index_type index) const;
    
    const hyperedge_index_type num_of_hyperedges() const;
    const hypernode_index_type num_of_hypernodes() const;
    const dof_index_type num_of_global_dofs() const;
};

#endif
