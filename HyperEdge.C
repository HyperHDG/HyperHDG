/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "HyperEdge.h"
#include <algorithm>
#include <cassert>

#include <iostream>

using namespace std;


template class HyperEdge_Cubic<1,1>;
template class HyperEdge_Cubic<1,2>;
template class HyperEdge_Cubic<1,3>;
template class HyperEdge_Cubic<2,2>;
template class HyperEdge_Cubic<2,3>;
template class HyperEdge_Cubic<3,3>;


template<unsigned int space_dim>
array<joint_index_type, 2> line_to_point_index(const array<unsigned int, space_dim>& num_lines, const hyperedge_index_type index)
{
  assert( num_lines.size() == space_dim );
  unsigned int orientation;
  hyperedge_index_type num_elements_in_direction;
  hyperedge_index_type number_with_lower_orientation = 0;
  hyperedge_index_type index_helper = index;
  
  array<joint_index_type, 2> point_indices;
  point_indices.fill(0);
  
  array<hyperedge_index_type, space_dim> num_lines_with_orientation;
  num_lines_with_orientation.fill(1);
  
  for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
      if (dim_m == dim_n)  num_lines_with_orientation[dim_m] *= num_lines[dim_n];
      else                 num_lines_with_orientation[dim_m] *= num_lines[dim_n] + 1;
  
  for ( orientation = 0;
        number_with_lower_orientation + num_lines_with_orientation[orientation] <= index ;
        ++orientation)
  {
    number_with_lower_orientation += num_lines_with_orientation[orientation];
    assert( orientation <= space_dim );
  }
  
  array<hyperedge_index_type, space_dim> local_indices;
  local_indices.fill(0);
  
  index_helper -= number_with_lower_orientation;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    num_elements_in_direction = num_lines[(dim + orientation) % space_dim] + (dim != 0);
    local_indices[(dim + orientation) % space_dim] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[(dim + orientation) % space_dim];
    index_helper /= num_elements_in_direction;
  }
  assert( index_helper == 0 );
  
  for(unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
  {
    unsigned int helper = 1;
    for (unsigned int dim_n = 0; dim_n < dim_m; ++dim_n)
      helper *= num_lines[dim_n] + 1;
    point_indices[0] += local_indices[dim_m] * helper;
    if (dim_m == orientation)  point_indices[1] += (local_indices[dim_m] + 1) * helper;
    else                       point_indices[1] += local_indices[dim_m] * helper;
  }
  
  return point_indices;
}


template<unsigned int space_dim>
array<joint_index_type, 4> square_to_line_index(const array<unsigned int, space_dim>& num_squares, const hyperedge_index_type index)
{
  assert( num_squares.size() == space_dim );
  unsigned int orientation;
  hyperedge_index_type num_elements_in_direction;
  hyperedge_index_type number_with_lower_orientation = 0;
  hyperedge_index_type index_helper = index;
  
  array<joint_index_type, 4> line_indices;
  line_indices.fill(0);
  
  array<hyperedge_index_type, space_dim> num_squares_with_orientation;
  num_squares_with_orientation.fill(1);
  
  for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
      if ( dim_m != dim_n )     num_squares_with_orientation[dim_m] *= num_squares[dim_n];
      else if (space_dim == 2)  num_squares_with_orientation[dim_m] *= num_squares[dim_n];
      else if (space_dim == 3)  num_squares_with_orientation[dim_m] *= num_squares[dim_n] + 1;
  
  for ( orientation = 0;
        number_with_lower_orientation + num_squares_with_orientation[orientation] <= index ;
        ++orientation)
  {
    number_with_lower_orientation += num_squares_with_orientation[orientation];
    assert( orientation <= space_dim );
  }
  
  array<hyperedge_index_type, space_dim> local_indices;
  local_indices.fill(0);
  
  index_helper -= number_with_lower_orientation;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    if (space_dim == 3)  num_elements_in_direction = num_squares[(dim + orientation) % space_dim] + (dim == 0);
    else                 num_elements_in_direction = num_squares[(dim + orientation) % space_dim];
    local_indices[(dim + orientation) % space_dim] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[(dim + orientation) % space_dim];
    index_helper /= num_elements_in_direction;
  }
  assert( index_helper == 0 );
  
  array<hyperedge_index_type, space_dim> num_lines_with_orientation;
  for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
      if (dim_m == dim_n)  num_lines_with_orientation[dim_m] *= num_squares[dim_n];
      else                 num_lines_with_orientation[dim_m] *= num_squares[dim_n] + 1;
  
  unsigned int local_line_index = 0;
  for (unsigned int line_orientation = 0; line_orientation < space_dim; ++line_orientation)
  {
    if (space_dim == 3 && line_orientation == orientation)  continue;
    number_with_lower_orientation = 0;
    for (unsigned int dim = 0; dim < line_orientation; ++dim)
      number_with_lower_orientation += num_lines_with_orientation[dim];
    for (int dim = space_dim - 1; dim >= 0; --dim)
    {
      num_elements_in_direction = num_squares[(dim + line_orientation) % space_dim] + (dim != 0);
      line_indices[local_line_index] *= num_elements_in_direction;
      line_indices[local_line_index] += local_indices[(dim + line_orientation) % space_dim];
      line_indices[local_line_index + 1] *= num_elements_in_direction;
      if (space_dim == 2) line_indices[local_line_index + 1] += local_indices[(dim + line_orientation) % space_dim] + (dim != 0);
      else  line_indices[local_line_index + 1] += local_indices[(dim + line_orientation) % space_dim] 
                                                  + ((dim + line_orientation) % space_dim != orientation && dim != 0);
    }
    line_indices[local_line_index] += number_with_lower_orientation;
    line_indices[local_line_index + 1] += number_with_lower_orientation;
    local_line_index += 2;
  }
  
  return line_indices;
}


template<unsigned int space_dim>
array<joint_index_type, 6> cube_to_square_index(const array<unsigned int, space_dim>& num_cubes, const hyperedge_index_type index)
{
  assert( num_cubes.size() == space_dim );
  hyperedge_index_type num_elements_in_direction;
  hyperedge_index_type number_with_lower_orientation = 0;
  hyperedge_index_type index_helper = index;
  
  array<joint_index_type, 6> square_indices;
  square_indices.fill(0);
  
  array<hyperedge_index_type, space_dim> local_indices;
  local_indices.fill(0);
  
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    num_elements_in_direction = num_cubes[dim];
    local_indices[dim] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[dim];
    index_helper /= num_elements_in_direction;
  }
  assert( index_helper == 0 );
  
  array<hyperedge_index_type, space_dim> num_squares_with_orientation;
  num_squares_with_orientation.fill(1);
  
  for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
      if ( dim_m != dim_n )     num_squares_with_orientation[dim_m] *= num_cubes[dim_n];
      else if (space_dim == 2)  num_squares_with_orientation[dim_m] *= num_cubes[dim_n];
      else if (space_dim == 3)  num_squares_with_orientation[dim_m] *= num_cubes[dim_n] + 1;
  
  unsigned int local_square_index = 0;
  for (unsigned int square_orientation = 0; square_orientation < space_dim; ++square_orientation)
  {
    number_with_lower_orientation = 0;
    for (unsigned int dim = 0; dim < square_orientation; ++dim)
      number_with_lower_orientation += num_squares_with_orientation[dim];
    for (int dim = space_dim - 1; dim >= 0; --dim)
    {
      num_elements_in_direction = num_cubes[(dim + square_orientation) % space_dim] + (dim == 0);
      square_indices[local_square_index] *= num_elements_in_direction;
      square_indices[local_square_index] += local_indices[(dim + square_orientation) % space_dim];
      square_indices[local_square_index + 1] *= num_elements_in_direction;
      square_indices[local_square_index + 1] += local_indices[(dim + square_orientation) % space_dim] + (dim == 0);
    }
    square_indices[local_square_index] += number_with_lower_orientation;
    square_indices[local_square_index + 1] += number_with_lower_orientation;
    local_square_index += 2;
  }
  
  return square_indices;
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
HyperEdge_Cubic<hyperedge_dim,space_dim>::
HyperEdge_Cubic(const hyperedge_index_type index, const array<unsigned int, space_dim>& num_elements,
                      const hyperedge_index_type num_of_hyperedges)
{
//  joint_indices_.resize(2*hyperedge_dim);
//  correct_joint_orientation_.resize(2*hyperedge_dim);
  for (unsigned int local_joint = 0; local_joint < 2 * hyperedge_dim; ++local_joint)
    correct_joint_orientation_[local_joint] = true;
  if constexpr ( hyperedge_dim == 1 )  joint_indices_ = line_to_point_index<space_dim>(num_elements, index);
  else if constexpr ( hyperedge_dim == 2 )  joint_indices_ = square_to_line_index<space_dim>(num_elements, index);
  else if constexpr ( hyperedge_dim == 3 )  joint_indices_ = cube_to_square_index<space_dim>(num_elements, index);    
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
const array<joint_index_type, 2*hyperedge_dim>&
HyperEdge_Cubic<hyperedge_dim,space_dim>::get_joint_indices() const
{
  return joint_indices_;
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
std::vector<double> HyperEdge_Cubic<hyperedge_dim,space_dim>::
abs_det_of_jacobian_at_quad(const vector<double>& local_quadrature) const
{
  vector<double> det_at_quad(local_quadrature.size(), 1.);
  return det_at_quad;
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
vector< vector<double> > HyperEdge_Cubic<hyperedge_dim,space_dim>::
inv_of_transposed_jacobian_at_quad(const vector<double>& local_quadrature) const
{
  vector< vector<double> > jac_at_quad(local_quadrature.size());
  for_each(jac_at_quad.begin(), jac_at_quad.end(), [](vector<double>& local_jac)
  {
    local_jac.resize(1);
    local_jac[0] = 1.;
  });
  return jac_at_quad;
}
