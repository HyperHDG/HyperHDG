/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "HyperNodeFactory.h"
#include <cassert>

using namespace std;


template class HyperNodeFactory<compute_n_dofs_per_node(1, 1, 1)>;
template class HyperNodeFactory<compute_n_dofs_per_node(2, 1, 1)>;
template class HyperNodeFactory<compute_n_dofs_per_node(2, 2, 1)>;
template class HyperNodeFactory<compute_n_dofs_per_node(2, 3, 1)>;
// template class HyperNodeFactory<compute_n_dofs_per_node(3, 1, 1)>;
template class HyperNodeFactory<compute_n_dofs_per_node(3, 2, 1)>;
template class HyperNodeFactory<compute_n_dofs_per_node(3, 3, 1)>;
// template class HyperNodeFactory<compute_n_dofs_per_node(1, 1, 2)>;
// template class HyperNodeFactory<compute_n_dofs_per_node(2, 1, 2)>;
template class HyperNodeFactory<compute_n_dofs_per_node(2, 2, 2)>;
template class HyperNodeFactory<compute_n_dofs_per_node(2, 3, 2)>;
// template class HyperNodeFactory<compute_n_dofs_per_node(3, 1, 2)>;
template class HyperNodeFactory<compute_n_dofs_per_node(3, 2, 2)>;
template class HyperNodeFactory<compute_n_dofs_per_node(3, 3, 2)>;
// template class HyperNodeFactory<compute_n_dofs_per_node(1, 1, 3)>;
// template class HyperNodeFactory<compute_n_dofs_per_node(2, 1, 3)>;
// template class HyperNodeFactory<compute_n_dofs_per_node(2, 2, 3)>;
template class HyperNodeFactory<compute_n_dofs_per_node(2, 3, 3)>;
// template class HyperNodeFactory<compute_n_dofs_per_node(3, 1, 3)>;
template class HyperNodeFactory<compute_n_dofs_per_node(3, 2, 3)>;
template class HyperNodeFactory<compute_n_dofs_per_node(3, 3, 3)>;


template <unsigned int n_dofs_per_node>
HyperNodeFactory<n_dofs_per_node>::
HyperNodeFactory(const hypernode_index_type num_of_hypernodes)
: num_of_hypernodes_(num_of_hypernodes) { }


template <unsigned int n_dofs_per_node>
HyperNodeFactory<n_dofs_per_node>::
HyperNodeFactory(const HyperNodeFactory<n_dofs_per_node>& other)
: num_of_hypernodes_(other.num_of_hypernodes_) { }


template <unsigned int n_dofs_per_node>
const hypernode_index_type HyperNodeFactory<n_dofs_per_node>::
num_of_hypernodes() const
{
  return num_of_hypernodes_;
}


template <unsigned int n_dofs_per_node>
const dof_index_type HyperNodeFactory<n_dofs_per_node>::
num_of_global_dofs() const
{
  return num_of_hypernodes_ * n_dofs_per_node;
}


template <unsigned int n_dofs_per_node>
array<dof_index_type, n_dofs_per_node> HyperNodeFactory<n_dofs_per_node>::
get_dof_indices(const hypernode_index_type hypernode_index) const
{
  dof_index_type initial_dof_index = hypernode_index * n_dofs_per_node;
  
  array<dof_index_type, n_dofs_per_node> dof_indices;
  for (unsigned int i = 0; i < n_dofs_per_node; ++i)
    dof_indices[i] = initial_dof_index + i;
  
  return dof_indices;
}


template <unsigned int n_dofs_per_node>
array<dof_value_type, n_dofs_per_node> HyperNodeFactory<n_dofs_per_node>::
get_dof_values(const hypernode_index_type hypernode_index, const vector<dof_value_type>& global_dof_vector) const
{
  dof_index_type initial_dof_index = hypernode_index * n_dofs_per_node;
  assert( initial_dof_index + n_dofs_per_node <= global_dof_vector.size() );
  
  array<dof_value_type, n_dofs_per_node> local_dof_values;
  for (unsigned int index = 0; index < n_dofs_per_node; ++index)
    local_dof_values[index] = global_dof_vector[initial_dof_index + index];
  
  return local_dof_values;
}


template <unsigned int n_dofs_per_node>
void HyperNodeFactory<n_dofs_per_node>::
add_to_dof_values(const hypernode_index_type hypernode_index, vector<dof_value_type>& global_dof_vector,
                  const array<dof_value_type, n_dofs_per_node>& local_dof_vector) const
{
  dof_index_type initial_dof_index = hypernode_index * n_dofs_per_node;
  assert( local_dof_vector.size() == n_dofs_per_node );
  assert( initial_dof_index + n_dofs_per_node <= global_dof_vector.size() );
  
  for(unsigned int index = 0; index < n_dofs_per_node; ++index)
    global_dof_vector[initial_dof_index + index] += local_dof_vector[index];
}


template <unsigned int n_dofs_per_node>
void HyperNodeFactory<n_dofs_per_node>::
set_dof_values(const hypernode_index_type hypernode_index, vector<dof_value_type>& global_dof_vector,
               const dof_value_type value) const
{
  dof_index_type initial_dof_index = hypernode_index * n_dofs_per_node;
  assert( initial_dof_index + n_dofs_per_node <= global_dof_vector.size() );
  for(unsigned int index = 0; index < n_dofs_per_node; ++index)
    global_dof_vector[initial_dof_index + index] = value;
}
