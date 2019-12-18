/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef HYPEREDGE_TOPOLOGY_H
#define HYPEREDGE_TOPOLOGY_H

#include "TypeDefs.h"
#include <array>

namespace Topology
{

template <unsigned int hyperedge_dim, unsigned int space_dim>
class HyperEdge_Cubic
{
  private:
    std::array<hypernode_index_type, 2*hyperedge_dim> hypernode_indices_;
    std::array<unsigned int, 2*hyperedge_dim> correct_hypernode_orientation_;
  public:
    HyperEdge_Cubic(const hyperedge_index_type index, const std::array<unsigned int, space_dim>& num_elements);
    const std::array<hypernode_index_type, 2*hyperedge_dim>& get_hypernode_indices() const;
};

}

#endif
