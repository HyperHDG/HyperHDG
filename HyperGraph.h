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

template < unsigned int amount_of_local_dofs, class Topology >
class HyperGraph
{
  private:
    const HyperNodeFactory<amount_of_local_dofs> hypernode_factory_;
    const Topology hyperedge_getter_;
  public:
    HyperGraph(const Topology& hyperedge_getter);
    
    const HyperNodeFactory<amount_of_local_dofs> hypernode_factory() const; // AR: No reference for performance?!
    const typename Topology::value_type get_hyperedge(const hyperedge_index_type index) const;
    
    const hypernode_index_type num_of_hypernodes() const;
    const hyperedge_index_type num_of_hyperedges() const;
    const dof_index_type num_of_global_dofs() const;
};

#endif
