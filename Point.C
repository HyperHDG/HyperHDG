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
#include <cmath>
#include <cassert>

using namespace std;


template class Point<1>;
template class Point<2>;
template class Point<3>;

/*
template double distance<1>(const Point<1>& left, const Point<1>& right);
template double distance<2>(const Point<2>& left, const Point<2>& right);
template double distance<3>(const Point<3>& left, const Point<3>& right);


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
Point<space_dim>::Point (const array<point_coord_type, space_dim>& coordinates) : coordinates_(coordinates)
{ 
  assert(coordinates.size() == dim);
  assert(0 < dim && dim < 4);
}


template<unsigned int space_dim>
point_coord_type& Point<space_dim>::operator[] (const unsigned int coord_entry)
{
  assert(0 <= coord_entry && coord_entry <= space_dim);
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
