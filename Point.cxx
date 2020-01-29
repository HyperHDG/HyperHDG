/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "Point.hxx"
#include "HyAssert.hxx"
#include <cmath>

using namespace std;
#include "Point.inst"


template<unsigned int space_dim>
pt_coord_t norm_2(const Point<space_dim>& point)
{
  pt_coord_t norm = 0.;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    norm += point.coordinate(dim) * point.coordinate(dim);
  return sqrt(norm);
}


template<unsigned int space_dim>
pt_coord_t distance_2(const Point<space_dim>& left, const Point<space_dim>& right)
{
  pt_coord_t distance = 0.;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    distance += (left.coordinate(dim) - right.coordinate(dim)) 
                * (left.coordinate(dim) - right.coordinate(dim));
  return sqrt(distance);
}


template<unsigned int space_dim>
Point<space_dim>::Point()
{ 
  coordinates_.fill(0.);
}


template<unsigned int space_dim>
Point<space_dim>::Point (const array<pt_coord_t, space_dim>& coordinates)
: coordinates_(coordinates)
{
  static_assert( 0 < space_dim && space_dim < 4 ,
                 "Dimension of a point is supposed to be between 1 and 3 (included), but is not!" );
  hy_assert( coordinates.size() == space_dim ,
             "Size of coordinates array is " << coordinates.size() << " and should be equal to the "
             << "template parameter dimension, which is " << space_dim << "." );
}


template<unsigned int space_dim>
Point<space_dim>::Point (const Point<space_dim>& other) // copy constructor
{
  static_assert( 0 < space_dim && space_dim < 4 ,
                 "Dimension of a point is supposed to be between 1 and 3 (included), but is not!" );
  hy_assert( coordinates_.size() == space_dim ,
             "Size of coordinates array is " << coordinates_.size() << " and should be equal to the"
             << " template parameter dimension, which is " << space_dim << "." );
  for(unsigned int dim = 0; dim < space_dim; ++dim)
    coordinates_[dim] = other.coordinate(dim);
}


template<unsigned int space_dim>
Point<space_dim>::Point (Point<space_dim>&& other) noexcept // move constructor
: coordinates_(std::move(other.coordinates_))
{
  static_assert( 0 < space_dim && space_dim < 4 ,
                 "Dimension of a point is supposed to be between 1 and 3 (included), but is not!" );
  hy_assert( coordinates_.size() == space_dim ,
             "Size of coordinates array is " << coordinates_.size() << " and should be equal to the"
             << " template parameter dimension, which is " << space_dim << "." );
}
  

template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator= (const Point<space_dim>& other) // copy assignement
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    coordinates_[dim] = other.coordinate(dim);
  return *this;
}


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator= (Point<space_dim>&& other) noexcept // move assignment
{
  std::swap(coordinates_, other.coordinates_);
  return *this;
}


template<unsigned int space_dim>
pt_coord_t& Point<space_dim>::operator[] (const unsigned int coord_entry)
{
  hy_assert( 0 <= coord_entry && coord_entry < space_dim ,
             "You can only access entries of a point's coordinates that have a non-negaitive index "
             << "that is smaller than the space dimension (which is " << space_dim << "). However, "
             << "you tried to access the " << coord_entry << "-th entry." );
  return coordinates_[coord_entry];
}


template<unsigned int space_dim>
pt_coord_t Point<space_dim>::coordinate(const unsigned int coord_entry) const
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


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator*=(const pt_coord_t scale_fac)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    coordinates_[dim] *= scale_fac;
  return *this;
}


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator/=(const pt_coord_t scale_denom)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    coordinates_[dim] /= scale_denom;
  return *this;
}


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator+=(const pt_coord_t additum)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    coordinates_[dim] += additum;
  return *this;
}


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator-=(const pt_coord_t subtractum)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    coordinates_[dim] -= subtractum;
  return *this;
}


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator+=(const Point<space_dim>& other)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    coordinates_[dim] += other.coordinate(dim);
  return *this;
}


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator-=(const Point<space_dim>& other)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    coordinates_[dim] -= other.coordinate(dim);
  return *this;
}


template<unsigned int space_dim>
Point<space_dim> Point<space_dim>::operator+(const Point<space_dim>& other) const
{
  Point difference(coordinates_);
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    difference[dim] += other.coordinate(dim);
  return difference;
}


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator+(Point<space_dim>&& other) noexcept
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    other[dim] += coordinates_[dim];
  return other;
}


template<unsigned int space_dim>
Point<space_dim> Point<space_dim>::operator-(const Point<space_dim>& other) const
{
  Point difference(coordinates_);
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    difference[dim] -= other.coordinate(dim);
  return difference;
}


template<unsigned int space_dim>
Point<space_dim>& Point<space_dim>::operator-(Point<space_dim>&& other) noexcept
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    other[dim] = coordinates_[dim] - other[dim];
  return other;
}
