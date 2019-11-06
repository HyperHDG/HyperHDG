/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "HyperGraphTopology.h"


template class HyperGraphTopology 
  < JointGetter_RegularQuad<1,1>, ConnectorGetter_RegularQuad<1,1>, Joint_RegularQuad, Connector_RegularQuad<1,1> >;
template class HyperGraphTopology 
  < JointGetter_RegularQuad<1,2>, ConnectorGetter_RegularQuad<1,2>, Joint_RegularQuad, Connector_RegularQuad<1,2> >;
template class HyperGraphTopology 
  < JointGetter_RegularQuad<1,3>, ConnectorGetter_RegularQuad<1,3>, Joint_RegularQuad, Connector_RegularQuad<1,3> >;
template class HyperGraphTopology 
  < JointGetter_RegularQuad<2,2>, ConnectorGetter_RegularQuad<2,2>, Joint_RegularQuad, Connector_RegularQuad<2,2> >;
template class HyperGraphTopology 
  < JointGetter_RegularQuad<2,3>, ConnectorGetter_RegularQuad<2,3>, Joint_RegularQuad, Connector_RegularQuad<2,3> >;
template class HyperGraphTopology 
  < JointGetter_RegularQuad<3,3>, ConnectorGetter_RegularQuad<3,3>, Joint_RegularQuad, Connector_RegularQuad<3,3> >;


template < class AbstractJointGetter, class AbstractConnectorGetter,
           class AbstractJoint,       class AbstractConnector >
HyperGraphTopology< AbstractJointGetter, AbstractConnectorGetter,
                    AbstractJoint,       AbstractConnector >::
HyperGraphTopology(const AbstractJointGetter& joint_getter,
                   const AbstractConnectorGetter& connector_getter)
: joint_getter_(joint_getter), connector_getter_(connector_getter) { }


template < class AbstractJointGetter, class AbstractConnectorGetter,
           class AbstractJoint,       class AbstractConnector >
const AbstractJoint
HyperGraphTopology< AbstractJointGetter, AbstractConnectorGetter,
                    AbstractJoint,       AbstractConnector >::
get_joint(const joint_index_type index) const
{
  return joint_getter_.get_joint(index);
}


template < class AbstractJointGetter, class AbstractConnectorGetter,
           class AbstractJoint,       class AbstractConnector >
const AbstractConnector 
HyperGraphTopology< AbstractJointGetter, AbstractConnectorGetter,
                    AbstractJoint,       AbstractConnector >::
get_connector(const connector_index_type index) const
{
  return connector_getter_.get_connector(index);
}


template < class AbstractJointGetter, class AbstractConnectorGetter,
           class AbstractJoint,       class AbstractConnector >
const joint_index_type
HyperGraphTopology< AbstractJointGetter, AbstractConnectorGetter,
                    AbstractJoint,       AbstractConnector >::
num_of_joints() const
{
  return joint_getter_.num_of_joints();
}


template < class AbstractJointGetter, class AbstractConnectorGetter,
           class AbstractJoint,       class AbstractConnector >
const connector_index_type
HyperGraphTopology< AbstractJointGetter, AbstractConnectorGetter,
                    AbstractJoint,       AbstractConnector >::
num_of_connectors() const
{
  return connector_getter_.num_of_connectors();
}


template < class AbstractJointGetter, class AbstractConnectorGetter,
           class AbstractJoint,       class AbstractConnector >
const dof_index_type 
HyperGraphTopology< AbstractJointGetter, AbstractConnectorGetter,
                    AbstractJoint,       AbstractConnector >::
num_of_global_dofs() const
{
  return joint_getter_.num_of_global_dofs();
}
