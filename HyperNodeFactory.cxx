/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "HyperNodeFactory.hxx"
#include "HyAssert.hxx"

using namespace std;
#include "HyperNodeFactory.inst"


template <unsigned int n_dofs_per_node>
HyperNodeFactory<n_dofs_per_node>::
HyperNodeFactory(const hyNode_index_t n_hypernodes)
: n_hypernodes_(n_hypernodes) { }


template <unsigned int n_dofs_per_node>
HyperNodeFactory<n_dofs_per_node>::
HyperNodeFactory(const HyperNodeFactory<n_dofs_per_node>& other)
: n_hypernodes_(other.n_hypernodes_) { }


template <unsigned int n_dofs_per_node>
const hyNode_index_t HyperNodeFactory<n_dofs_per_node>::
n_hypernodes() const
{
  return n_hypernodes_;
}


template <unsigned int n_dofs_per_node>
const dof_index_type HyperNodeFactory<n_dofs_per_node>::
n_global_dofs() const
{
  return n_hypernodes_ * n_dofs_per_node;
}


template <unsigned int n_dofs_per_node>
array<dof_index_type, n_dofs_per_node> HyperNodeFactory<n_dofs_per_node>::
get_dof_indices(const hyNode_index_t hypernode_index) const
{
  dof_index_type initial_dof_index = hypernode_index * n_dofs_per_node;
  
  array<dof_index_type, n_dofs_per_node> dof_indices;
  for (unsigned int i = 0; i < n_dofs_per_node; ++i)
    dof_indices[i] = initial_dof_index + i;
  
  return dof_indices;
}


template <unsigned int n_dofs_per_node>
array<dof_value_t, n_dofs_per_node> HyperNodeFactory<n_dofs_per_node>::
get_dof_values(const hyNode_index_t hypernode_index,
               const vector<dof_value_t>& global_dof_vector) const
{
  dof_index_type initial_dof_index = hypernode_index * n_dofs_per_node;
  hy_assert( initial_dof_index >= 0
               && initial_dof_index + n_dofs_per_node <= global_dof_vector.size() ,
             "The initial dof index = " << initial_dof_index << "should be non-negative. Moreover, "
             << "the final index = " << initial_dof_index + n_dofs_per_node << " must not exceed "
             << "the size of the vector of global degrees of freedom." );
  
  array<dof_value_t, n_dofs_per_node> local_dof_values;
  for (unsigned int index = 0; index < n_dofs_per_node; ++index)
    local_dof_values[index] = global_dof_vector[initial_dof_index + index];
  
  return local_dof_values;
}


template <unsigned int n_dofs_per_node>
void HyperNodeFactory<n_dofs_per_node>::
add_to_dof_values(const hyNode_index_t hypernode_index,
                  vector<dof_value_t>& global_dof_vector,
                  const array<dof_value_t, n_dofs_per_node>& local_dof_vector) const
{
  dof_index_type initial_dof_index = hypernode_index * n_dofs_per_node;
  hy_assert( local_dof_vector.size() == n_dofs_per_node ,
             "The size of the local dof vector is " << local_dof_vector.size() << ", but should be "
             << "equal to the amount of local dofs, which is " << n_dofs_per_node << "." );
  hy_assert( initial_dof_index >= 0 &&
               initial_dof_index + n_dofs_per_node <= global_dof_vector.size() ,
             "The initial dof index = " << initial_dof_index << "should be non-negative. Moreover, "
             << "the final index = " << initial_dof_index + n_dofs_per_node << " must not exceed "
             << "the size of the vector of global degrees of freedom." );
  
  for(unsigned int index = 0; index < n_dofs_per_node; ++index)
    global_dof_vector[initial_dof_index + index] += local_dof_vector[index];
}


template <unsigned int n_dofs_per_node>
void HyperNodeFactory<n_dofs_per_node>::
set_dof_values(const hyNode_index_t hypernode_index,
               vector<dof_value_t>& global_dof_vector,
               const dof_value_t value) const
{
  dof_index_type initial_dof_index = hypernode_index * n_dofs_per_node;
  hy_assert( initial_dof_index >= 0 &&
               initial_dof_index + n_dofs_per_node <= global_dof_vector.size() ,
             "The initial dof index = " << initial_dof_index << "should be non-negative. Moreover, "
             << "the final index = " << initial_dof_index + n_dofs_per_node << " must not exceed "
             << "the size of the vector of global degrees of freedom." );
  
  for(unsigned int index = 0; index < n_dofs_per_node; ++index)
    global_dof_vector[initial_dof_index + index] = value;
}
