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
#include "JointGetter.h"
#include "Topology.h"

template < class AbstractJointGetter, class Topology,
           class AbstractJoint >
class HyperGraph
{
  private:
    const AbstractJointGetter joint_getter_;
    const Topology hyperedge_getter_;
  public:
    HyperGraph(const AbstractJointGetter& joint_getter,
               const Topology& hyperedge_getter);
    HyperGraph();
    
    const AbstractJoint get_joint(const joint_index_type index) const;
    const typename Topology::value_type get_hyperedge(const hyperedge_index_type index) const;
    
    const joint_index_type num_of_joints() const;
    const hyperedge_index_type num_of_hyperedges() const;
    const dof_index_type num_of_global_dofs() const;
};

#endif
