/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef TOPOLOGY_H
#define TOPOLOGY_H

#include "TypeDefs.h"
#include "HyperEdge.h"
#include <array>
#include <vector>


template <unsigned int hyperedge_dim, unsigned int space_dim>
class HyperGraph_Cubic
{
  private:
    std::array<unsigned int, space_dim> num_elements_;
    hyperedge_index_type num_of_hyperedges_;
    joint_index_type num_of_vertices_;
  public:
    typedef HyperEdge_Cubic<hyperedge_dim, space_dim> value_type;
    HyperGraph_Cubic(const std::vector<int>& num_elements);
    HyperGraph_Cubic(const std::array<unsigned int, space_dim>& num_elements);
    HyperGraph_Cubic(const HyperGraph_Cubic<hyperedge_dim,space_dim>& other);
    
    const HyperEdge_Cubic<hyperedge_dim, space_dim> get_hyperedge(const hyperedge_index_type index) const;
    
    const hyperedge_index_type num_of_hyperedges() const;
    const joint_index_type num_of_vertices() const;
};

#endif
