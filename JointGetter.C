/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "JointGetter.h"
#include <vector>
#include <cassert>

using namespace std;


template class JointGetter_RegularQuad<1,1>;
template class JointGetter_RegularQuad<1,2>;
template class JointGetter_RegularQuad<1,3>;
template class JointGetter_RegularQuad<2,2>;
template class JointGetter_RegularQuad<2,3>;
template class JointGetter_RegularQuad<3,3>;


template <unsigned int connector_dim, unsigned int space_dim>
JointGetter_RegularQuad<connector_dim,space_dim>::
JointGetter_RegularQuad(const unsigned int amount_of_local_dofs,
                        const unsigned int num_of_elements_in_x_dir,
                        const unsigned int num_of_elements_in_y_dir,
                        const unsigned int num_of_elements_in_z_dir)
: amount_of_local_dofs_(amount_of_local_dofs)
{
  static_assert( connector_dim >= 1, "Domains must have dimension larger than or equal to 1!" );
  static_assert( space_dim >= connector_dim, "A domain cannot live within a smaller space!" );
  static_assert( space_dim <= 3, "Only spaces up to dimension 3 are implemented!" );
  assert( amount_of_local_dofs_ >= 1 );
  assert( connector_dim > 1 || amount_of_local_dofs_ == 1 );
  
  vector<unsigned int> num_elements;
  num_elements.resize(space_dim);
  num_elements[0] = num_of_elements_in_x_dir;
  if (space_dim > 1)  num_elements[1] = num_of_elements_in_y_dir;
  if (space_dim > 2)  num_elements[2] = num_of_elements_in_z_dir;
  
  num_of_joints_ = 1;
  if (connector_dim == 1)
  {
    num_of_joints_ *= num_of_elements_in_x_dir + 1;
    if (space_dim > 1)  num_of_joints_ *= num_of_elements_in_y_dir + 1;
    if (space_dim > 2)  num_of_joints_ *= num_of_elements_in_z_dir + 1;
  }
  else if ( connector_dim == space_dim )
  {
    num_of_joints_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      num_of_joints_ += helper;
    }
  }
  else if ( connector_dim == space_dim - 1 )
  {
    num_of_joints_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      num_of_joints_ += helper;
    }
  }
  else  assert( 0 == 1 );
  
  assert( num_of_joints_ > 0 );
}


template <unsigned int connector_dim, unsigned int space_dim>
JointGetter_RegularQuad<connector_dim,space_dim>::
JointGetter_RegularQuad(const JointGetter_RegularQuad<connector_dim,space_dim>& other)
: num_of_joints_(other.num_of_joints_), amount_of_local_dofs_(other.amount_of_local_dofs_) { }


template <unsigned int connector_dim, unsigned int space_dim>
const Joint_RegularQuad JointGetter_RegularQuad<connector_dim,space_dim>::
get_joint(const unsigned int index) const
{
  assert( index < num_of_joints_ );
  return Joint_RegularQuad(index, amount_of_local_dofs_);
}


template <unsigned int connector_dim, unsigned int space_dim>
const joint_index_type JointGetter_RegularQuad<connector_dim,space_dim>::
num_of_joints() const
{
  return num_of_joints_;
}


template <unsigned int connector_dim, unsigned int space_dim>
const dof_index_type JointGetter_RegularQuad<connector_dim,space_dim>::
num_of_global_dofs() const
{
  return num_of_joints_ * amount_of_local_dofs_;
}
