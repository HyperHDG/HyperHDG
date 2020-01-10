/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "HDGHyperGraph.h"
#include <cassert>


#include "HDGHyperGraph.inst"


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
HDGHyperGraph(const TopoT& hypergraph_topology)
: hypergraph_topology_(hypergraph_topology), hypergraph_geometry_(hypergraph_topology),
  hypernode_factory_(hypergraph_topology.num_of_hypernodes())
{
  static_assert( TopoT::hyperedge_dimension() == GeomT::hyperedge_dimension() );
  assert( hypernode_factory_.num_of_hypernodes() == hypergraph_topology.num_of_hypernodes() );
  assert( hypernode_factory_.num_of_hypernodes() >= 2 );
  assert( hypergraph_topology.num_of_hyperedges() != 0 );
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
HDGHyperGraph(const constructor_value_type& construction_data)
: hypergraph_topology_(construction_data.topology), hypergraph_geometry_(construction_data.geometry),
  hypernode_factory_(hypergraph_topology_.num_of_hypernodes())
{
  static_assert( TopoT::hyperedge_dimension() == GeomT::hyperedge_dimension() );
  assert( hypernode_factory_.num_of_hypernodes() == hypergraph_topology.num_of_hypernodes() );
  assert( hypernode_factory_.num_of_hypernodes() >= 2 );
  assert( hypergraph_topology_.num_of_hyperedges() != 0 );
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
HDGHyperGraph(const typename TopoT::constructor_value_type& construct_topo,
              const typename GeomT::constructor_value_type& construct_geom)
: hypergraph_topology_(construct_topo), hypergraph_geometry_(construct_geom),
  hypernode_factory_(hypergraph_topology_.num_of_hypernodes())
{
  static_assert( TopoT::hyperedge_dimension() == GeomT::hyperedge_dimension() );
  assert( hypernode_factory_.num_of_hypernodes() == hypergraph_topology.num_of_hypernodes() );
  assert( hypernode_factory_.num_of_hypernodes() >= 2 );
  assert( hypergraph_topology_.num_of_hyperedges() != 0 );
}



template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
const typename HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::value_type
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
operator[](const hyperedge_index_type index) const
{
  return value_type(hyperedge_topology(index), hyperedge_geometry(index)); 
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
typename HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::iterator
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
begin() const
{
  return HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::iterator(*this, 0); 
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
typename HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::iterator
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
end() const
{
  return HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::iterator(*this, num_of_hyperedges()); 
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
const HyperNodeFactory<n_dofs_per_node>
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
hypernode_factory() const
{
  return hypernode_factory_;
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
const typename TopoT::value_type
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
hyperedge_topology(const hyperedge_index_type index) const
{
  return hypergraph_topology_.get_hyperedge(index);
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
const typename GeomT::value_type 
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
hyperedge_geometry(const hyperedge_index_type index) const
{
  return hypergraph_geometry_.get_hyperedge(index);
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
const hyperedge_index_type
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
num_of_hyperedges() const
{
  return hypergraph_topology_.num_of_hyperedges();
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
const hypernode_index_type
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
num_of_hypernodes() const
{
  return hypernode_factory_.num_of_hypernodes();
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
const dof_index_type 
HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
num_of_global_dofs() const
{
  return hypernode_factory_.num_of_global_dofs();
}

/*
template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
constexpr unsigned int HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
hyperedge_dimension()
{
  return TopoT::hyperedge_dimension();
}


template < unsigned int n_dofs_per_node, class TopoT, class GeomT >
constexpr unsigned int HDGHyperGraph< n_dofs_per_node, TopoT, GeomT >::
space_dimension()
{
  return TopoT::space_dimension();
}
*/
