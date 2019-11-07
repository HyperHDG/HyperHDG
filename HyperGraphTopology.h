/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef DOMAIN_TOPOLOGY_H
#define DOMAIN_TOPOLOGY_H

#include "TypeDefs.h"
#include "JointGetter.h"
#include "ConnectorGetter.h"

template < class AbstractJointGetter, class AbstractConnectorGetter,
           class AbstractJoint >
class HyperGraphTopology
{
  private:
    const AbstractJointGetter joint_getter_;
    const AbstractConnectorGetter connector_getter_;
  public:
    HyperGraphTopology(const AbstractJointGetter& joint_getter,
                       const AbstractConnectorGetter& connector_getter);
    HyperGraphTopology();
    
    const AbstractJoint get_joint(const joint_index_type index) const;
    const typename AbstractConnectorGetter::value_type get_connector(const connector_index_type index) const;
    
    const joint_index_type num_of_joints() const;
    const connector_index_type num_of_connectors() const;
    const dof_index_type num_of_global_dofs() const;
};

#endif
