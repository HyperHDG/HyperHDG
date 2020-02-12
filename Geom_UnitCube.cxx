/*!*************************************************************************************************
 * \file    Geom_UnitCube.C
 * \brief   Implementation of Geom_UnitCube.h (for more details refer to that file).
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/

#include "Geom_UnitCube.hxx"
#include "HyAssert.hxx"

using namespace std;
using namespace Geometry;

typedef unsigned int hyNode_index_t;

/*
 * Relevant instances of Unit cubes!
 */


template class UnitCube< 1, 1 >;
template class UnitCube< 1, 2 >;
template class UnitCube< 1, 3 >;
template class UnitCube< 2, 2 >;
template class UnitCube< 2, 3 >;
template class UnitCube< 3, 3 >;


/*
 * HyperEdge functions!
 */


template<unsigned int space_dim>
array<Point<space_dim>, 2> line_to_points(const array<unsigned int, space_dim>& num_lines, const hyEdge_index_t index)
{
  hy_assert( num_lines.size() == space_dim , "The size of the handed over parmeter does not fit!" );
  unsigned int orientation;
  hyEdge_index_t num_elements_in_direction;
  hyEdge_index_t number_with_lower_orientation = 0;
  hyEdge_index_t index_helper = index;
  
  array<Point<space_dim>, 2> point_indices;
  point_indices.fill(Point<space_dim>());
  
  array<hyEdge_index_t, space_dim> num_lines_with_orientation;
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
  
  array<hyEdge_index_t, space_dim> local_indices;
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
    double helper = num_lines[dim_m];
    point_indices[0][dim_m] = local_indices[dim_m] / helper;
    if (dim_m == orientation)  point_indices[1][dim_m] = (local_indices[dim_m] + 1) / helper;
    else                       point_indices[1][dim_m] = local_indices[dim_m] / helper;
  }
  
  return point_indices;
}


template<unsigned int space_dim>
array<hyNode_index_t, 4> square_to_line_index(const array<unsigned int, space_dim>& num_squares, const hyEdge_index_t index)
{
  hy_assert( num_squares.size() == space_dim , "The size of the handed over parmeter does not fit!" );
  unsigned int orientation;
  hyEdge_index_t num_elements_in_direction;
  hyEdge_index_t number_with_lower_orientation = 0;
  hyEdge_index_t index_helper = index;
  
  array<hyNode_index_t, 4> line_indices;
  line_indices.fill(0);
  
  array<hyEdge_index_t, space_dim> num_squares_with_orientation;
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
  
  array<hyEdge_index_t, space_dim> local_indices;
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
  hy_assert( index_helper == 0 , "No lines should be left any more!" );
  
  array<hyEdge_index_t, space_dim> num_lines_with_orientation;
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
array<hyNode_index_t, 6> cube_to_square_index(const array<unsigned int, space_dim>& num_cubes, const hyEdge_index_t index)
{
  hy_assert( num_cubes.size() == space_dim , "The size of the handed over parmeter does not fit!" );
  hyEdge_index_t num_elements_in_direction;
  hyEdge_index_t number_with_lower_orientation = 0;
  hyEdge_index_t index_helper = index;
  
  array<hyNode_index_t, 6> square_indices;
  square_indices.fill(0);
  
  array<hyEdge_index_t, space_dim> local_indices;
  local_indices.fill(0);
  
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    num_elements_in_direction = num_cubes[dim];
    local_indices[dim] = index_helper % num_elements_in_direction;
    index_helper -= local_indices[dim];
    index_helper /= num_elements_in_direction;
  }
  hy_assert( index_helper == 0 , "No lines should be left any more!" );
  
  array<hyEdge_index_t, space_dim> num_squares_with_orientation;
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


template <unsigned int hyEdge_dim, unsigned int space_dim>
UnitCube<hyEdge_dim,space_dim>::hyEdge::
hyEdge(const hyEdge_index_t index, const array<unsigned int, space_dim>& num_elements)
{
  if constexpr ( hyEdge_dim == 1 )
  {
    points_ = line_to_points<space_dim>(num_elements, index);
  }
  else if constexpr ( hyEdge_dim == 2 )
  {
//    hy_assert( 0 == 1 , "This should not be executed since elasticity is not yet defined for real HYPERgraphs!" );
    array<hyNode_index_t, 4> line_indices = square_to_line_index<space_dim>(num_elements, index);
    array<Point<space_dim>, 2> points1 = line_to_points<space_dim>(num_elements, line_indices[0]);
    array<Point<space_dim>, 2> points2 = line_to_points<space_dim>(num_elements, line_indices[1]);
    points_[0] = points1[0];  points_[1] = points1[1];  points_[2] = points2[0];  points_[3] = points2[1];
  }
  else if constexpr ( hyEdge_dim == 3 )
  {
//    hy_assert( 0 == 1 , "This should not be executed since elasticity is not yet defined for real HYPERgraphs!" );
    array<hyNode_index_t, 6> square_indices = cube_to_square_index<space_dim>(num_elements, index);
    array<hyNode_index_t, 4> line_indices1 = square_to_line_index<space_dim>(num_elements, square_indices[0]);
    array<hyNode_index_t, 4> line_indices2 = square_to_line_index<space_dim>(num_elements, square_indices[1]);
    array<Point<space_dim>, 2> points1 = line_to_points<space_dim>(num_elements, line_indices1[0]);
    array<Point<space_dim>, 2> points2 = line_to_points<space_dim>(num_elements, line_indices1[1]);
    array<Point<space_dim>, 2> points3 = line_to_points<space_dim>(num_elements, line_indices2[0]);
    array<Point<space_dim>, 2> points4 = line_to_points<space_dim>(num_elements, line_indices2[1]);
    points_[0] = points1[0];  points_[1] = points1[1];  points_[2] = points2[0];  points_[3] = points2[1];
    points_[4] = points3[0];  points_[5] = points3[1];  points_[6] = points4[0];  points_[7] = points4[1];
  }
//  sort(points_.begin(), points_.end());
}


template <unsigned int hyEdge_dim, unsigned int space_dim>
Point<space_dim> UnitCube<hyEdge_dim,space_dim>::hyEdge::
normal(unsigned int index) const
{
  Point<space_dim> normal;
  if (index == 0)  normal = points_[0] - points_[1];
  else             normal = points_[1] - points_[0];
  normal /= norm_2(normal);
  return normal;
}

/*
const array<hyNode_index_t, 2*hyEdge_dim>&
UnitCube<hyEdge_dim,space_dim>::get_hyNode_indices() const
{
  return hyNode_indices_;
}


template <unsigned int hyEdge_dim, unsigned int space_dim>
std::vector<double> UnitCube<hyEdge_dim,space_dim>::
abs_det_of_jacobian_at_quad(const vector<double>& local_quadrature) const
{
  vector<double> det_at_quad(local_quadrature.size(), 1.);
  return det_at_quad;
}


template <unsigned int hyEdge_dim, unsigned int space_dim>
vector< vector<double> > UnitCube<hyEdge_dim,space_dim>::
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
*/
