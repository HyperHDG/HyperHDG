/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "Point.h"
#include "HyAssert.h"
#include <cmath>

using namespace std;
#include "Point.inst"


/*
template<unsigned int space_dim>
double distance(const Point<space_dim>& left, const Point<space_dim>& right)
{
  point_coord_type distance = 0.;
  for (unsigned int n = 0; n < space_dim; ++n)
    distance += (left[n] - right[n]) * (left[n] - right[n]);
  return sqrt(distance);
}
*/

template<unsigned int space_dim>
Point<space_dim>::Point()
{ 
  coordinates_.fill(0.);
}


template<unsigned int space_dim>
Point<space_dim>::Point (const array<point_coord_type, space_dim>& coordinates)
: coordinates_(coordinates)
{
  static_assert( 0 < space_dim && space_dim < 4 ,
                 "Dimension of a point is supposed to be between 1 and 3 (included), but is not!" );
  hy_assert( coordinates.size() == space_dim ,
             "Size of coordinates array is " << coordinates.size() << " and should be equal to the "
             << "template parameter dimension, which is " << space_dim << "." );
}


template<unsigned int space_dim>
point_coord_type& Point<space_dim>::operator[] (const unsigned int coord_entry)
{
  hy_assert( 0 <= coord_entry && coord_entry < space_dim ,
             "You can only access entries of a point's coordinates that have a non-negaitive index "
             << "that is smaller than the space dimension (which is " << space_dim << "). However, "
             << "you tried to access the " << coord_entry << "-th entry." );
  return coordinates_[coord_entry];
}


template<unsigned int space_dim>
bool Point<space_dim>::operator== (const Point<space_dim>& other_point) const
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    if (coordinates_[dim] != other_point.coordinates_[dim])  return false;
  return true;
}


template<unsigned int space_dim>
bool Point<space_dim>::operator!= (const Point<space_dim>& other_point) const
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    if (coordinates_[dim] != other_point.coordinates_[dim])  return true;
  return false;
}


template<unsigned int space_dim>
bool Point<space_dim>::operator<(const Point<space_dim>& other_point) const
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    if (coordinates_[dim] < other_point.coordinates_[dim])      return true;
    else if (coordinates_[dim] > other_point.coordinates_[dim]) return false;
  return false;
}
