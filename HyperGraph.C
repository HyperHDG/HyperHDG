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


template class HyperGraph
  < JointGetter_RegularQuad<1,1>, HyperGraph_Cubic< 1, 1 >, Joint_RegularQuad>;
template class HyperGraph
  < JointGetter_RegularQuad<1,2>, HyperGraph_Cubic< 1, 2 >, Joint_RegularQuad>;
template class HyperGraph 
  < JointGetter_RegularQuad<1,3>, HyperGraph_Cubic< 1, 3 >, Joint_RegularQuad>;
template class HyperGraph
  < JointGetter_RegularQuad<2,2>, HyperGraph_Cubic< 2, 2 >, Joint_RegularQuad>;
template class HyperGraph
  < JointGetter_RegularQuad<2,3>, HyperGraph_Cubic< 2, 3 >, Joint_RegularQuad>;
template class HyperGraph 
  < JointGetter_RegularQuad<3,3>, HyperGraph_Cubic< 3, 3 >, Joint_RegularQuad>;


template < class AbstractJointGetter, class Topology,
           class AbstractJoint >
HyperGraph< AbstractJointGetter, Topology,
                    AbstractJoint >::
HyperGraph(const AbstractJointGetter& joint_getter,
                   const Topology& hyperedge_getter)
: joint_getter_(joint_getter), hyperedge_getter_(hyperedge_getter) { }


template < class AbstractJointGetter, class Topology,
           class AbstractJoint >
const AbstractJoint
HyperGraph< AbstractJointGetter, Topology,
                    AbstractJoint >::
get_joint(const joint_index_type index) const
{
  return joint_getter_.get_joint(index);
}


template < class AbstractJointGetter, class Topology,
           class AbstractJoint >
const typename Topology::value_type
HyperGraph< AbstractJointGetter, Topology,
                    AbstractJoint >::
get_hyperedge(const hyperedge_index_type index) const
{
  return hyperedge_getter_.get_hyperedge(index);
}


template < class AbstractJointGetter, class Topology,
           class AbstractJoint >
const joint_index_type
HyperGraph< AbstractJointGetter, Topology,
                    AbstractJoint >::
num_of_joints() const
{
  return joint_getter_.num_of_joints();
}


template < class AbstractJointGetter, class Topology,
           class AbstractJoint >
const hyperedge_index_type
HyperGraph< AbstractJointGetter, Topology,
                    AbstractJoint >::
num_of_hyperedges() const
{
  return hyperedge_getter_.num_of_hyperedges();
}


template < class AbstractJointGetter, class Topology,
           class AbstractJoint >
const dof_index_type 
HyperGraph< AbstractJointGetter, Topology,
                    AbstractJoint >::
num_of_global_dofs() const
{
  return joint_getter_.num_of_global_dofs();
}
