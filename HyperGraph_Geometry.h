/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef HYPERGRAPH_GEOMETRY_H
#define HYPERGRAPH_GEOMETRY_H

#include "TypeDefs.h"
#include "HyperEdge_Geometry.h"
#include "HyperGraph_Topology.h"
#include <array>

namespace Geometry
{

template <unsigned int hyperedge_dim, unsigned int space_dim>
class HyperGraph_Cubic_UnitCube
{
  private:
    std::array<unsigned int, space_dim> num_elements_;
  public:
    typedef HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim> value_type;
    HyperGraph_Cubic_UnitCube(const Topology::HyperGraph_Cubic<hyperedge_dim,space_dim>& other);
    
    const HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim> get_hyperedge(const hyperedge_index_type index) const;
}; // end class HyperGraph_Cubic_UnitCube

} // end namespace Geometry

#endif // end ifndef HYPERGRAPH_GEOMETRY_H
