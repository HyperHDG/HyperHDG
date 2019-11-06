/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#ifndef JOINT_GETTER_H
#define JOINT_GETTER_H

#include "Joint.h"
#include "TypeDefs.h"

/*
template <class AbstractJoint>
class JointGetter
{
  private:
    vector<AbstractJoint> joints;
  public:
    AbstractJoint& get_joint(const unsigned int index);
    unsigned int num_of_joints();
}
*/
template <unsigned int connector_dim, unsigned int space_dim>
class JointGetter_RegularQuad
{
  private:
    joint_index_type num_of_joints_;
    const unsigned int amount_of_local_dofs_;
  public:
    JointGetter_RegularQuad(const unsigned int amount_of_local_dofs,
                            const unsigned int num_of_elements_in_x_dir,
                            const unsigned int num_of_elements_in_y_dir = 0,
                            const unsigned int num_of_elements_in_z_dir = 0);
    JointGetter_RegularQuad
      (const JointGetter_RegularQuad<connector_dim,space_dim>& other);
    const Joint_RegularQuad get_joint(const unsigned int index) const;
    const joint_index_type num_of_joints() const;
    const dof_index_type num_of_global_dofs() const;
};

#endif
