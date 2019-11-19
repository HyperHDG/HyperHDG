/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "HyperGraph.h"
#include <cassert>


template class HyperGraph < local_dof_amount_node(1 , 1), HyperGraph_Cubic< 1, 1 > >;
// template class HyperGraph < local_dof_amount_node(1 , 2), HyperGraph_Cubic< 1, 1 > >;
// template class HyperGraph < local_dof_amount_node(1 , 3), HyperGraph_Cubic< 1, 1 > >;
template class HyperGraph < local_dof_amount_node(1 , 1), HyperGraph_Cubic< 1, 2 > >;
// template class HyperGraph < local_dof_amount_node(1 , 2), HyperGraph_Cubic< 1, 2 > >;
// template class HyperGraph < local_dof_amount_node(1 , 3), HyperGraph_Cubic< 1, 2 > >;
template class HyperGraph < local_dof_amount_node(1 , 1), HyperGraph_Cubic< 1, 3 > >;
// template class HyperGraph < local_dof_amount_node(1 , 2), HyperGraph_Cubic< 1, 3 > >;
// template class HyperGraph < local_dof_amount_node(1 , 3), HyperGraph_Cubic< 1, 3 > >;
template class HyperGraph < local_dof_amount_node(2 , 1), HyperGraph_Cubic< 2, 2 > >;
template class HyperGraph < local_dof_amount_node(2 , 2), HyperGraph_Cubic< 2, 2 > >;
template class HyperGraph < local_dof_amount_node(2 , 3), HyperGraph_Cubic< 2, 2 > >;
template class HyperGraph < local_dof_amount_node(2 , 1), HyperGraph_Cubic< 2, 3 > >;
template class HyperGraph < local_dof_amount_node(2 , 2), HyperGraph_Cubic< 2, 3 > >;
template class HyperGraph < local_dof_amount_node(2 , 3), HyperGraph_Cubic< 2, 3 > >;
template class HyperGraph < local_dof_amount_node(3 , 1), HyperGraph_Cubic< 3, 3 > >;
template class HyperGraph < local_dof_amount_node(3 , 2), HyperGraph_Cubic< 3, 3 > >;
template class HyperGraph < local_dof_amount_node(3 , 3), HyperGraph_Cubic< 3, 3 > >;


template < unsigned int amount_of_local_dofs, class Topology >
HyperGraph< amount_of_local_dofs, Topology >::
HyperGraph(const Topology& hyperedge_getter)
: vertex_factory_(hyperedge_getter.num_of_vertices()), hyperedge_getter_(hyperedge_getter)
{
  assert( vertex_factory_.num_of_vertices() == hyperedge_getter.num_of_vertices() );
  assert( vertex_factory_.num_of_vertices() >= 2 );
  assert( hyperedge_getter.num_of_hyperedges() != 0 );
}


template < unsigned int amount_of_local_dofs, class Topology >
const VertexFactory<amount_of_local_dofs>
HyperGraph< amount_of_local_dofs, Topology >::
vertex_factory() const
{
  return vertex_factory_;
}


template < unsigned int amount_of_local_dofs, class Topology >
const typename Topology::value_type
HyperGraph< amount_of_local_dofs, Topology >::
get_hyperedge(const hyperedge_index_type index) const
{
  return hyperedge_getter_.get_hyperedge(index);
}


template < unsigned int amount_of_local_dofs, class Topology >
const joint_index_type
HyperGraph< amount_of_local_dofs, Topology >::
num_of_vertices() const
{
  return vertex_factory_.num_of_vertices();
}


template < unsigned int amount_of_local_dofs, class Topology >
const hyperedge_index_type
HyperGraph< amount_of_local_dofs, Topology >::
num_of_hyperedges() const
{
  return hyperedge_getter_.num_of_hyperedges();
}


template < unsigned int amount_of_local_dofs, class Topology >
const dof_index_type 
HyperGraph< amount_of_local_dofs, Topology >::
num_of_global_dofs() const
{
  return vertex_factory_.num_of_global_dofs();
}
