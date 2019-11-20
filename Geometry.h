/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef GEOMETRY_H
#define GEOMETRY_H

#include "TypeDefs.h"
#include "HyperEdge_Geometry.h"
#include "Topology.h"
#include <array>


template <unsigned int hyperedge_dim, unsigned int space_dim>
class HyperGraph_Cubic_UnitCube
{
  private:
    std::array<unsigned int, space_dim> num_elements_;
//    point_index_type num_of_vertices_;
//    hypernode_index_type num_of_hypernodes_;
  public:
    typedef HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim> value_type;
    HyperGraph_Cubic_UnitCube(const HyperGraph_Cubic<hyperedge_dim,space_dim>& other);
    
    const HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim> get_hyperedge(const hyperedge_index_type index) const;
    
//    const hyperedge_index_type num_of_hyperedges() const;
//    const point_index_type num_of_vertices() const;
};

#endif
