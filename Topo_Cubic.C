/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */

#include "Topo_Cubic.h"
#include "HyAssert.h"
#include <algorithm>

using namespace std;
using namespace Topology;
#include "Topo_Cubic.inst"


/*
 * HyperEdge functions!
 */


template<unsigned int space_dim>
array<hypernode_index_type, 2> line_to_point_index
(const array<unsigned int, space_dim>& num_lines, const hyperedge_index_type index)
{
  hy_assert( num_lines.size() == space_dim , "The size of the handed over parmeter does not fit!" );
  unsigned int orientation;
  hyperedge_index_type num_elements_in_direction;
  hyperedge_index_type number_with_lower_orientation = 0;
  hyperedge_index_type index_helper = index;
  
  array<hypernode_index_type, 2> point_indices;
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
    hy_assert( orientation <= space_dim , "Orientation is a space_dim and connot exceed it." );
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
  hy_assert( index_helper == 0 , "No lines should be left any more!" );
  
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
array<hypernode_index_type, 4> square_to_line_index
(const array<unsigned int, space_dim>& num_squares, const hyperedge_index_type index)
{
  hy_assert( num_squares.size() == space_dim , "The size of the handed over parmeter does not fit!" );
  unsigned int orientation;
  hyperedge_index_type num_elements_in_direction;
  hyperedge_index_type number_with_lower_orientation = 0;
  hyperedge_index_type index_helper = index;
  
  array<hypernode_index_type, 4> line_indices;
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
    hy_assert( orientation <= space_dim , "Orientation is a space_dim and connot exceed it." );
  }
  
  array<hyperedge_index_type, space_dim> local_indices;
  local_indices.fill(0);
  
  index_helper -= number_with_lower_orientation;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    if (space_dim == 3)  num_elements_in_direction = num_squares[(dim + orientation) % space_dim]
                                                       + (dim == 0);
    else                 num_elements_in_direction = num_squares[(dim + orientation) % space_dim];
    local_indices[(dim + orientation) % space_dim] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[(dim + orientation) % space_dim];
    index_helper /= num_elements_in_direction;
  }
  hy_assert( index_helper == 0 , "No squares should be left any more!" );
  
  array<hyperedge_index_type, space_dim> num_lines_with_orientation;
  num_lines_with_orientation.fill(1);
  
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
array<hypernode_index_type, 6> cube_to_square_index
(const array<unsigned int, space_dim>& num_cubes, const hyperedge_index_type index)
{
  hy_assert( num_cubes.size() == space_dim , "The size of the handed over parmeter does not fit!" );
  hyperedge_index_type num_elements_in_direction;
  hyperedge_index_type number_with_lower_orientation = 0;
  hyperedge_index_type index_helper = index;
  
  array<hypernode_index_type, 6> square_indices;
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
  hy_assert( index_helper == 0 , "No cubes should be left any more!" );
  
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
Cubic<hyperedge_dim,space_dim>::hyperedge::
hyperedge(const hyperedge_index_type index, const array<unsigned int, space_dim>& num_elements)
{
  for (unsigned int local_hypernode = 0; local_hypernode < 2 * hyperedge_dim; ++local_hypernode)
    correct_hypernode_orientation_[local_hypernode] = 1;
  if constexpr ( hyperedge_dim == 1 )       hypernode_indices_ = line_to_point_index<space_dim>(num_elements, index);
  else if constexpr ( hyperedge_dim == 2 )  hypernode_indices_ = square_to_line_index<space_dim>(num_elements, index);
  else if constexpr ( hyperedge_dim == 3 )  hypernode_indices_ = cube_to_square_index<space_dim>(num_elements, index);    
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
const array<hypernode_index_type, 2*hyperedge_dim>&
Cubic<hyperedge_dim,space_dim>::hyperedge::get_hypernode_indices() const
{
  return hypernode_indices_;
}


/*
 * HyperGraph functions!
 */


template <unsigned int hyperedge_dim, unsigned int space_dim>
Cubic<hyperedge_dim,space_dim>::
Cubic(const array<unsigned int, space_dim>& num_elements)
: num_elements_(num_elements)
{
  static_assert( hyperedge_dim >= 1, "Domains must have dimension larger than or equal to 1!" );
  static_assert( space_dim >= hyperedge_dim, "A domain cannot live within a smaller space!" );
  static_assert( space_dim <= 3, "Only spaces up to dimension 3 are implemented!" );
    
  // Set n_hyperedges_
  n_hyperedges_ = 1;
  if ( hyperedge_dim == space_dim )
    for (unsigned int dim = 0; dim < space_dim; ++dim)  n_hyperedges_ *= num_elements[dim];
  else if ( hyperedge_dim == space_dim - 1 )
  {
    n_hyperedges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      n_hyperedges_ += helper;
    }
  }
  else if ( hyperedge_dim == space_dim - 2 )
  {
    n_hyperedges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      n_hyperedges_ += helper;
    }
  }
  else  hy_assert( 0 == 1 , "Internal error when trying to construct a hypergraph topology.");
  hy_assert( n_hyperedges_ > 0 , "An empty hypergraph is being constructed." );
  
  // Set n_hypernodes
  n_hypernodes_ = 1;
  if (hyperedge_dim == 1)
  {
    n_hypernodes_ *= num_elements[0] + 1;
    if (space_dim > 1)  n_hypernodes_ *= num_elements[1] + 1;
    if (space_dim > 2)  n_hypernodes_ *= num_elements[2] + 1;
  }
  else if ( hyperedge_dim == space_dim )
  {
    n_hypernodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      n_hypernodes_ += helper;
    }
  }
  else if ( hyperedge_dim == space_dim - 1 )
  {
    n_hypernodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      n_hypernodes_ += helper;
    }
  }
  else  hy_assert( 0 == 1 , "Internal error when trying to construct a hypergraph topology." );
  hy_assert( n_hypernodes_ > 0 , "An empty hypergraph is being constructed." );
}

template <unsigned int hyperedge_dim, unsigned int space_dim>
Cubic<hyperedge_dim,space_dim>::
Cubic(const constructor_value_type& num_elements)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim) num_elements_[dim] = num_elements[dim];
  
  static_assert( hyperedge_dim >= 1, "Domains must have dimension larger than or equal to 1!" );
  static_assert( space_dim >= hyperedge_dim, "A domain cannot live within a smaller space!" );
  static_assert( space_dim <= 3, "Only spaces up to dimension 3 are implemented!" );
    
  // Set n_hyperedges_
  n_hyperedges_ = 1;
  if ( hyperedge_dim == space_dim )
    for (unsigned int dim = 0; dim < space_dim; ++dim)  n_hyperedges_ *= num_elements[dim];
  else if ( hyperedge_dim == space_dim - 1 )
  {
    n_hyperedges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      n_hyperedges_ += helper;
    }
  }
  else if ( hyperedge_dim == space_dim - 2 )
  {
    n_hyperedges_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      n_hyperedges_ += helper;
    }
  }
  else  hy_assert( 0 == 1 , "Internal error when trying to construct a hypergraph topology." );
  hy_assert( n_hyperedges_ > 0 , "An empty hypergraph is being constructed." );
  
  // Set n_hypernodes
  n_hypernodes_ = 1;
  if (hyperedge_dim == 1)
  {
    n_hypernodes_ *= num_elements[0] + 1;
    if (space_dim > 1)  n_hypernodes_ *= num_elements[1] + 1;
    if (space_dim > 2)  n_hypernodes_ *= num_elements[2] + 1;
  }
  else if ( hyperedge_dim == space_dim )
  {
    n_hypernodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n] + 1;
        else                 helper *= num_elements[dim_n];
      n_hypernodes_ += helper;
    }
  }
  else if ( hyperedge_dim == space_dim - 1 )
  {
    n_hypernodes_ = 0;
    for (unsigned int dim_m = 0; dim_m < space_dim; ++dim_m)
    {
      int helper = 1;
      for (unsigned int dim_n = 0; dim_n < space_dim; ++dim_n)
        if (dim_m == dim_n)  helper *= num_elements[dim_n];
        else                 helper *= num_elements[dim_n] + 1;
      n_hypernodes_ += helper;
    }
  }
  else  hy_assert( 0 == 1 , "Internal error when trying to construct a hypergraph topology." );
  hy_assert( n_hypernodes_ > 0 , "An empty hypergraph is being constructed." );
}

template <unsigned int hyperedge_dim, unsigned int space_dim>
Cubic<hyperedge_dim,space_dim>::
Cubic(const Cubic<hyperedge_dim,space_dim>& other)
: num_elements_(other.num_elements_), n_hyperedges_(other.n_hyperedges_),
  n_hypernodes_(other.n_hypernodes_) { }


template <unsigned int hyperedge_dim, unsigned int space_dim>
const typename Cubic<hyperedge_dim,space_dim>::hyperedge
Cubic<hyperedge_dim,space_dim>::
get_hyperedge(const hyperedge_index_type index) const
{
  hy_assert( index >= 0 && index < n_hyperedges_ ,
             "The index of an hyperedge must be non-negative and smaller than the total amount of "
             << "hyperedges, which is " << n_hyperedges_ << ". Nonetheless, the " << index <<
             "-th hyperedge is tried to be accessed." );
  return hyperedge(index, num_elements_);
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
const array<unsigned int, space_dim>& 
Cubic<hyperedge_dim,space_dim>::
num_elements() const
{
  return num_elements_;
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
const hyperedge_index_type
Cubic<hyperedge_dim,space_dim>::
n_hyperedges() const
{
  return n_hyperedges_;
}

template <unsigned int hyperedge_dim, unsigned int space_dim>
const hypernode_index_type
Cubic<hyperedge_dim,space_dim>::
n_hypernodes() const
{
  return n_hypernodes_;
}

/*
template <unsigned int hyperedge_dim, unsigned int space_dim>
constexpr unsigned int Cubic<hyperedge_dim,space_dim>::hyperedge_dimension()
{
  return hyperedge_dim;
}

template <unsigned int hyperedge_dim, unsigned int space_dim>
constexpr unsigned int Cubic<hyperedge_dim,space_dim>::space_dimension()
{
  return space_dim;
}
*/
