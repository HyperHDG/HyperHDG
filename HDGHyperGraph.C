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


template < unsigned int amount_of_local_dofs, class Topology, class Geometry >
HDGHyperGraph< amount_of_local_dofs, Topology, Geometry >::
HDGHyperGraph(const Topology& hyperedge_getter)
: hypernode_factory_(hyperedge_getter.num_of_hypernodes()), hyperedge_getter_(hyperedge_getter),
  hyperedge_geometry_(hyperedge_getter)
{
  assert( hypernode_factory_.num_of_hypernodes() == hyperedge_getter.num_of_hypernodes() );
  assert( hypernode_factory_.num_of_hypernodes() >= 2 );
  assert( hyperedge_getter.num_of_hyperedges() != 0 );
}


template < unsigned int amount_of_local_dofs, class Topology, class Geometry >
const HyperNodeFactory<amount_of_local_dofs>
HDGHyperGraph< amount_of_local_dofs, Topology, Geometry >::
hypernode_factory() const
{
  return hypernode_factory_;
}


template < unsigned int amount_of_local_dofs, class Topology, class Geometry >
const typename Topology::value_type
HDGHyperGraph< amount_of_local_dofs, Topology, Geometry >::
get_hyperedge(const hyperedge_index_type index) const
{
  return hyperedge_getter_.get_hyperedge(index);
}


template < unsigned int amount_of_local_dofs, class Topology, class Geometry >
const typename Geometry::value_type 
HDGHyperGraph< amount_of_local_dofs, Topology, Geometry >::
get_hyperedge_geometry(const hyperedge_index_type index) const
{
  return hyperedge_geometry_.get_hyperedge(index);
}


template < unsigned int amount_of_local_dofs, class Topology, class Geometry >
const hypernode_index_type
HDGHyperGraph< amount_of_local_dofs, Topology, Geometry >::
num_of_hypernodes() const
{
  return hypernode_factory_.num_of_hypernodes();
}


template < unsigned int amount_of_local_dofs, class Topology, class Geometry >
const hyperedge_index_type
HDGHyperGraph< amount_of_local_dofs, Topology, Geometry >::
num_of_hyperedges() const
{
  return hyperedge_getter_.num_of_hyperedges();
}


template < unsigned int amount_of_local_dofs, class Topology, class Geometry >
const dof_index_type 
HDGHyperGraph< amount_of_local_dofs, Topology, Geometry >::
num_of_global_dofs() const
{
  return hypernode_factory_.num_of_global_dofs();
}
