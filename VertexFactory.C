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


template class VertexFactory<local_dof_amount(1, 1)>;
template class VertexFactory<local_dof_amount(2, 1)>;
template class VertexFactory<local_dof_amount(2, 2)>;
template class VertexFactory<local_dof_amount(2, 3)>;
// template class VertexFactory<local_dof_amount(3, 1)>;
template class VertexFactory<local_dof_amount(3, 2)>;
template class VertexFactory<local_dof_amount(3, 3)>;


/*
constexpr const unsigned int local_dof_amount(const unsigned int hyperedge_dim, const unsigned int polynomial_degree)
{
  return pow(polynomial_degree + 1 , hyperedge_dim - 1);
}
*/

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
vector<dof_index_type> VertexFactory<amount_of_local_dofs>::
get_dof_indices(const joint_index_type joint_index) const
{
  dof_index_type initial_dof_index = joint_index * amount_of_local_dofs;
  vector<dof_index_type> dof_indices(amount_of_local_dofs ,initial_dof_index);
  unsigned int local_increment = 0;
  
  for_each(dof_indices.begin(), dof_indices.end(), [&local_increment](dof_index_type& glob_index)
  { glob_index += local_increment++; });
  
  return dof_indices;
}


template <unsigned int amount_of_local_dofs>
vector<dof_value_type> VertexFactory<amount_of_local_dofs>::
get_dof_values(const joint_index_type joint_index, const vector<dof_value_type>& global_dof_vector) const
{
  dof_index_type initial_dof_index = joint_index * amount_of_local_dofs;
  assert( initial_dof_index + amount_of_local_dofs_ <= global_dof_vector.size() );
  vector<dof_value_type>::const_iterator first = global_dof_vector.begin() + initial_dof_index;
  vector<dof_value_type>::const_iterator last = first + amount_of_local_dofs;
  vector<dof_value_type> local_dof_values(first, last);
  return local_dof_values;
}


template <unsigned int amount_of_local_dofs>
void VertexFactory<amount_of_local_dofs>::
add_to_dof_values(const joint_index_type joint_index, vector<dof_value_type>& global_dof_vector,
                  const vector<dof_value_type>& local_dof_vector) const
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
