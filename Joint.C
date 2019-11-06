/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "Joint.h"
#include <algorithm>
#include <cassert>

using namespace std;


Joint_RegularQuad::Joint_RegularQuad(const unsigned int index, const unsigned int amount_of_local_dofs)
: initial_dof_index_(index * amount_of_local_dofs), amount_of_local_dofs_(amount_of_local_dofs) { }


vector<dof_index_type> Joint_RegularQuad::get_dof_indices() const
{
  vector<dof_index_type> dof_indices(amount_of_local_dofs_,initial_dof_index_);
  unsigned int local_increment = 0;
  
  for_each(dof_indices.begin(), dof_indices.end(), [&local_increment](dof_index_type& glob_index)
  { glob_index += local_increment++; });
  
  return dof_indices;
}


vector<dof_value_type> Joint_RegularQuad::
get_dof_values(const vector<dof_value_type>& global_dof_vector) const
{
  assert( initial_dof_index_ + amount_of_local_dofs_ <= global_dof_vector.size() );
  vector<dof_value_type>::const_iterator first = global_dof_vector.begin() + initial_dof_index_;
  vector<dof_value_type>::const_iterator last = first + amount_of_local_dofs_;
  vector<dof_value_type> local_dof_values(first, last);
  return local_dof_values;
}


void Joint_RegularQuad::add_to_dof_values(vector<dof_value_type>& global_dof_vector,
                                          const vector<dof_value_type>& local_dof_vector) const
{
  assert( local_dof_vector.size() == amount_of_local_dofs_ );
  assert( initial_dof_index_ + amount_of_local_dofs_ <= global_dof_vector.size() );
  for(dof_index_type index = 0; index < amount_of_local_dofs_; ++index)
    global_dof_vector[index + initial_dof_index_] += local_dof_vector[index];
}


void Joint_RegularQuad::set_dof_values(vector<dof_value_type>& global_dof_vector,
                                       const double value) const
{
  assert( initial_dof_index_ + amount_of_local_dofs_ <= global_dof_vector.size() );
  for(dof_index_type index = 0; index < amount_of_local_dofs_; ++index)
    global_dof_vector[index + initial_dof_index_] = value;
}


vector<double> Joint_RegularQuad::abs_det_of_jacobian_at_quad
  (const vector<double>& local_quadrature) const
{
  vector<double> det_at_quad(local_quadrature.size(), 1.);
  return det_at_quad;
}
