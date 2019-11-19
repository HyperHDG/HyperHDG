/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "VertexFactory.h"
#include <cassert>

using namespace std;


template class VertexFactory<local_dof_amount_node(1, 1)>;
template class VertexFactory<local_dof_amount_node(2, 1)>;
template class VertexFactory<local_dof_amount_node(2, 2)>;
template class VertexFactory<local_dof_amount_node(2, 3)>;
// template class VertexFactory<local_dof_amount_node(3, 1)>;
template class VertexFactory<local_dof_amount_node(3, 2)>;
template class VertexFactory<local_dof_amount_node(3, 3)>;


template <unsigned int amount_of_local_dofs>
VertexFactory<amount_of_local_dofs>::
VertexFactory(const joint_index_type num_of_vertices)
: num_of_vertices_(num_of_vertices) { }


template <unsigned int amount_of_local_dofs>
VertexFactory<amount_of_local_dofs>::
VertexFactory(const VertexFactory<amount_of_local_dofs>& other)
: num_of_vertices_(other.num_of_vertices_) { }


template <unsigned int amount_of_local_dofs>
const joint_index_type VertexFactory<amount_of_local_dofs>::
num_of_vertices() const
{
  return num_of_vertices_;
}


template <unsigned int amount_of_local_dofs>
const dof_index_type VertexFactory<amount_of_local_dofs>::
num_of_global_dofs() const
{
  return num_of_vertices_ * amount_of_local_dofs;
}


template <unsigned int amount_of_local_dofs>
array<dof_index_type, amount_of_local_dofs> VertexFactory<amount_of_local_dofs>::
get_dof_indices(const joint_index_type joint_index) const
{
  dof_index_type initial_dof_index = joint_index * amount_of_local_dofs;
  
  array<dof_index_type, amount_of_local_dofs> dof_indices;
  for (unsigned int i = 0; i < amount_of_local_dofs; ++i)
    dof_indices[i] = initial_dof_index + i;
  
  return dof_indices;
}


template <unsigned int amount_of_local_dofs>
array<dof_value_type, amount_of_local_dofs> VertexFactory<amount_of_local_dofs>::
get_dof_values(const joint_index_type joint_index, const vector<dof_value_type>& global_dof_vector) const
{
  dof_index_type initial_dof_index = joint_index * amount_of_local_dofs;
  assert( initial_dof_index + amount_of_local_dofs <= global_dof_vector.size() );
  
  array<dof_value_type, amount_of_local_dofs> local_dof_values;
  for (unsigned int i = 0; i < amount_of_local_dofs; ++i)
    local_dof_values[i] = global_dof_vector[initial_dof_index + i];
  
  return local_dof_values;
}


template <unsigned int amount_of_local_dofs>
void VertexFactory<amount_of_local_dofs>::
add_to_dof_values(const joint_index_type joint_index, vector<dof_value_type>& global_dof_vector,
                  const array<dof_value_type, amount_of_local_dofs>& local_dof_vector) const
{
  dof_index_type initial_dof_index = joint_index * amount_of_local_dofs;
  assert( local_dof_vector.size() == amount_of_local_dofs );
  assert( initial_dof_index + amount_of_local_dofs <= global_dof_vector.size() );
  
  for(dof_index_type index = 0; index < amount_of_local_dofs; ++index)
    global_dof_vector[index + initial_dof_index] += local_dof_vector[index];
}


template <unsigned int amount_of_local_dofs>
void VertexFactory<amount_of_local_dofs>::
set_dof_values(const joint_index_type joint_index, vector<dof_value_type>& global_dof_vector,
               const double value) const
{
  dof_index_type initial_dof_index = joint_index * amount_of_local_dofs;
  assert( initial_dof_index + amount_of_local_dofs <= global_dof_vector.size() );
  for(dof_index_type index = 0; index < amount_of_local_dofs; ++index)
    global_dof_vector[index + initial_dof_index] = value;
}
