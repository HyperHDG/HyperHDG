/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "Geometry.h"
#include <cassert>

using namespace std;


template class HyperGraph_Cubic_UnitCube< 1, 1 >;
template class HyperGraph_Cubic_UnitCube< 1, 2 >;
template class HyperGraph_Cubic_UnitCube< 1, 3 >;
template class HyperGraph_Cubic_UnitCube< 2, 2 >;
template class HyperGraph_Cubic_UnitCube< 2, 3 >;
template class HyperGraph_Cubic_UnitCube< 3, 3 >;


template <unsigned int hyperedge_dim, unsigned int space_dim>
HyperGraph_Cubic_UnitCube<hyperedge_dim,space_dim>::
HyperGraph_Cubic_UnitCube(const HyperGraph_Cubic<hyperedge_dim,space_dim>& other)
: num_elements_(other.num_elements()) { }


template <unsigned int hyperedge_dim, unsigned int space_dim>
const HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim>
HyperGraph_Cubic_UnitCube<hyperedge_dim,space_dim>::
get_hyperedge(const hyperedge_index_type index) const
{
  assert ( index < num_of_hyperedges_ );
  return HyperEdge_Cubic_UnitCube<hyperedge_dim,space_dim>(index, num_elements_);
}
