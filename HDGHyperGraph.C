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


template class HDGHyperGraph < local_dof_amount_node(1 , 1), Topology::HyperGraph_Cubic< 1, 1 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 1 > >;
// template class HDGHyperGraph < local_dof_amount_node(1 , 2), Topology::HyperGraph_Cubic< 1, 1 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 1 > >;
// template class HDGHyperGraph < local_dof_amount_node(1 , 3), Topology::HyperGraph_Cubic< 1, 1 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 1 > >;
template class HDGHyperGraph < local_dof_amount_node(1 , 1), Topology::HyperGraph_Cubic< 1, 2 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 2 > >;
// template class HDGHyperGraph < local_dof_amount_node(1 , 2), Topology::HyperGraph_Cubic< 1, 2 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 2 > >;
// template class HDGHyperGraph < local_dof_amount_node(1 , 3), Topology::HyperGraph_Cubic< 1, 2 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 2 > >;
template class HDGHyperGraph < local_dof_amount_node(1 , 1), Topology::HyperGraph_Cubic< 1, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 3 > >;
// template class HDGHyperGraph < local_dof_amount_node(1 , 2), Topology::HyperGraph_Cubic< 1, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 3 > >;
// template class HDGHyperGraph < local_dof_amount_node(1 , 3), Topology::HyperGraph_Cubic< 1, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 1, 3 > >;
template class HDGHyperGraph < local_dof_amount_node(2 , 1), Topology::HyperGraph_Cubic< 2, 2 >, Geometry::HyperGraph_Cubic_UnitCube< 2, 2 > >;
template class HDGHyperGraph < local_dof_amount_node(2 , 2), Topology::HyperGraph_Cubic< 2, 2 >, Geometry::HyperGraph_Cubic_UnitCube< 2, 2 > >;
template class HDGHyperGraph < local_dof_amount_node(2 , 3), Topology::HyperGraph_Cubic< 2, 2 >, Geometry::HyperGraph_Cubic_UnitCube< 2, 2 > >;
template class HDGHyperGraph < local_dof_amount_node(2 , 1), Topology::HyperGraph_Cubic< 2, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 2, 3 > >;
template class HDGHyperGraph < local_dof_amount_node(2 , 2), Topology::HyperGraph_Cubic< 2, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 2, 3 > >;
template class HDGHyperGraph < local_dof_amount_node(2 , 3), Topology::HyperGraph_Cubic< 2, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 2, 3 > >;
template class HDGHyperGraph < local_dof_amount_node(3 , 1), Topology::HyperGraph_Cubic< 3, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 3, 3 > >;
template class HDGHyperGraph < local_dof_amount_node(3 , 2), Topology::HyperGraph_Cubic< 3, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 3, 3 > >;
template class HDGHyperGraph < local_dof_amount_node(3 , 3), Topology::HyperGraph_Cubic< 3, 3 >, Geometry::HyperGraph_Cubic_UnitCube< 3, 3 > >;


template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::
HDGHyperGraph(const TopoT& hypergraph_topology)
: hypernode_factory_(hypergraph_topology.num_of_hypernodes()), hypergraph_topology_(hypergraph_topology),
  hypergraph_geometry_(hypergraph_topology)
{
  assert( hypernode_factory_.num_of_hypernodes() == hypergraph_topology.num_of_hypernodes() );
  assert( hypernode_factory_.num_of_hypernodes() >= 2 );
  assert( hypergraph_topology.num_of_hyperedges() != 0 );
}


template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
const typename HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::value_type
HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::
operator[](const hyperedge_index_type index) const
{
  return value_type(hyperedge_topology(index), hyperedge_geometry(index)); 
}


template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
const HyperNodeFactory<amount_of_local_dofs>
HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::
hypernode_factory() const
{
  return hypernode_factory_;
}


template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
const typename TopoT::value_type
HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::
hyperedge_topology(const hyperedge_index_type index) const
{
  return hypergraph_topology_.get_hyperedge(index);
}


template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
const typename GeomT::value_type 
HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::
hyperedge_geometry(const hyperedge_index_type index) const
{
  return hypergraph_geometry_.get_hyperedge(index);
}


template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
const hyperedge_index_type
HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::
num_of_hyperedges() const
{
  return hypergraph_topology_.num_of_hyperedges();
}


template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
const hypernode_index_type
HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::
num_of_hypernodes() const
{
  return hypernode_factory_.num_of_hypernodes();
}


template < unsigned int amount_of_local_dofs, class TopoT, class GeomT >
const dof_index_type 
HDGHyperGraph< amount_of_local_dofs, TopoT, GeomT >::
num_of_global_dofs() const
{
  return hypernode_factory_.num_of_global_dofs();
}
