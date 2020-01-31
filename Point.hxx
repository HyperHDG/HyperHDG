/*!*************************************************************************************************
 * \file    Point.hxx
 * \brief   This class implements a point in a d-dimensional space.
 * 
 * This file implements a point in a \f$d\f$-dimensional space, where the \f$d\f$ is given by the
 * template parameter \c space_dim.
 * 
 * \tparam  space_dim           The dimension of the space, the object is located in.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/

#ifndef POINT_HXX
#define POINT_HXX

#include "TypeDefs.hxx"
#include "HyAssert.hxx"
#include <cmath>
#include <array>
#include <ostream>

/*!*************************************************************************************************
 * \brief   This class implements a point in a d-dimensional space.
 * 
 * This class implements a point in a \f$d\f$-dimensional space, where the \f$d\f$ is given by the
 * template parameter \c space_dim.
 * 
 * \tparam  space_dim           The dimension of the space, the object is located in.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dim>
class Point
{
  private:
    /*!*********************************************************************************************
     * \brief   Array containint the coordinates of the point.
     * 
     * A \c std::array conatining the i-th coordinate of the point as its i-th entry.
     **********************************************************************************************/
    std::array<pt_coord_t, space_dim> coordinates_;
  public:
    /*!*********************************************************************************************
     * \brief   Empty constructor for a point.
     * 
     * Fills coordinates of the point with zeros.
     **********************************************************************************************/
    Point();
    /*!*********************************************************************************************
     * \brief   Construct point from array of coordinates.
     * 
     * Fills the point's array of coordinates with the input parameter. 
     * 
     * \param   coordinates   A \c std::array containing the coordinates of the point.
     **********************************************************************************************/
    Point(const std::array<pt_coord_t, space_dim>& coordinates);
    
    Point(const Point<space_dim>& other); // copy constructor
    Point(Point<space_dim>&& other) noexcept; // move constructor
    Point<space_dim>& operator= (const Point<space_dim>& other); // copy assignement
    Point<space_dim>& operator= (Point<space_dim>&& other) noexcept; // move assignment
    
    /*!*********************************************************************************************
     * \brief   Return single coordinate of a constant point.
     * 
     * \param   coord_entry   An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    \c pt_coord_t describing the coord_entry'th
     *                        coordinate.
     **********************************************************************************************/
    pt_coord_t operator[](const unsigned int coord_entry) const;
    /*!*********************************************************************************************
     * \brief   Return reference to single coordinate of a point.
     * 
     * \param   coord_entry   An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    A reference to a \c pt_coord_t describing the coord_entry'th
     *                        coordinate.
     **********************************************************************************************/
    pt_coord_t& operator[](const unsigned int coord_entry);
    /*!*********************************************************************************************
     * \brief   Return reference to const single coordinate of a point.
     * 
     * \param   coord_entry   An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    A reference to a \c pt_coord_t describing the coord_entry'th
     *                        coordinate.
     **********************************************************************************************/
    pt_coord_t coordinate(const unsigned int coord_entry) const;
    /*!*********************************************************************************************
     * \brief   Find out whether two points have (exactly) the same coordinates.
     * 
     * This function compares the point to another point and returns true if and only if both points
     * have exactly (that is not only with respect to some rounding errors) the same coordinates.
     * 
     * \param   other_point   Another \c Point<space_dim> that is to be dicriminate from.
     * \retval  isEqual       A \c boolean which is true if both points have the same coordinates.
     **********************************************************************************************/
    bool operator==(const Point<space_dim>& other_point) const;
    /*!*********************************************************************************************
     * \brief   Find out whether two points have (exactly) the same coordinates.
     * 
     * This function compares the point to another point and returns false if and only if both
     * points have exactly (that is not only with respect to some rounding errors) the same
     * coordinates.
     * 
     * \param   other_point   Another \c Point<space_dim> that is to be dicriminate from.
     * \retval  isEqual       A \c boolean which is false if both points have the same coordinates.
     **********************************************************************************************/
    bool operator!=(const Point<space_dim>& other_point) const;
    /*!*********************************************************************************************
     * \brief   Find out whether the point is "smaller than" another point.
     * 
     * This function compares the point to another point and returns true if and only if the lowest
     * ranked coordinate (according to the coordinate index) where the both points are not equal of
     * the given point is smaller than that of the other point. It is false, if both points are
     * equal.
     * 
     * \param   other_point   Another \c Point<space_dim> that is to be dicriminate from.
     * \retval  isEqual       A \c boolean which is true if the left point is strictly smaller than
     *                        the right one.
     **********************************************************************************************/
    bool operator<(const Point<space_dim>& other_point) const;
    
    Point<space_dim>& operator*=(const pt_coord_t scale_fac);
    Point<space_dim>& operator/=(const pt_coord_t scale_denom);
    Point<space_dim>& operator+=(const pt_coord_t additum);
    Point<space_dim>& operator-=(const pt_coord_t subtractum);
    
    Point<space_dim>& operator+=(const Point<space_dim>& other);
    Point<space_dim>& operator-=(const Point<space_dim>& other);
    
    Point<space_dim> operator+(const Point<space_dim>& other) const;
    Point<space_dim>& operator+(Point<space_dim>&& other) noexcept;
    Point<space_dim> operator-(const Point<space_dim>& other) const;
    Point<space_dim>& operator-(Point<space_dim>&& other) noexcept;
    
}; // end of class Point

template<unsigned int space_dim>
pt_coord_t norm_2(const Point<space_dim>& point)
{
  pt_coord_t norm = 0.;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    norm += point.coordinate(dim) * point.coordinate(dim);
  return std::sqrt(norm);
}


template<unsigned int space_dim>
pt_coord_t distance_2(const Point<space_dim>& left, const Point<space_dim>& right)
{
  pt_coord_t distance = 0.;
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    distance += (left.coordinate(dim) - right.coordinate(dim)) 
                * (left.coordinate(dim) - right.coordinate(dim));
  return std::sqrt(distance);
}


template<unsigned int space_dim>
std::ostream& operator<< (std::ostream& stream, const Point<space_dim>& pt)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    stream << " " << pt.coordinate(dim) << " ";
  return stream;
}


template<unsigned int space_dim>
Point<space_dim>::Point()
{ 
  coordinates_.fill(0.);
}


template<unsigned int space_dim>
Point<space_dim>::Point (const std::array<pt_coord_t, space_dim>& coordinates)
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
pt_coord_t Point<space_dim>::operator[] (const unsigned int coord_entry) const
{
  hy_assert( 0 <= coord_entry && coord_entry < space_dim ,
             "You can only access entries of a point's coordinates that have a non-negaitive index "
             << "that is smaller than the space dimension (which is " << space_dim << "). However, "
             << "you tried to access the " << coord_entry << "-th entry." );
  return coordinates_[coord_entry];
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

// End of implementation of member functions

template<unsigned int space_dim>
pt_coord_t norm_2(const Point<space_dim>& point); 


template<unsigned int space_dim>
pt_coord_t distance_2(const Point<space_dim>& left, const Point<space_dim>& right);


template<unsigned int space_dim>
std::ostream& operator<< (std::ostream& stream, const Point<space_dim>& pt);

#endif // end of ifndef POINT_HXX
