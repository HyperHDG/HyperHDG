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

#include <TypeDefs.hxx>
#include <HyAssert.hxx>

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
    Point()
    {
      coordinates_.fill(0.);
    }
    /*!*********************************************************************************************
     * \brief   Construct point from array of coordinates.
     * 
     * Fills the point's array of coordinates with the input parameter. 
     * 
     * \param   coordinates   A \c std::array containing the coordinates of the point.
     **********************************************************************************************/
    Point(const std::array<pt_coord_t, space_dim>& coordinates)
    : coordinates_(coordinates) { }
    /*!*********************************************************************************************
     * \brief   Copy constructor.
     **********************************************************************************************/
    Point(const Point<space_dim>& other)
    : coordinates_(other.coordinates_) { }
    /*!*********************************************************************************************
     * \brief   Move constructor.
     **********************************************************************************************/
    Point(Point<space_dim>&& other) noexcept
    : coordinates_(std::move(other.coordinates_)) { }
    /*!*********************************************************************************************
     * \brief   Copy assignment.
     **********************************************************************************************/
    Point<space_dim>& operator= (const Point<space_dim>& other)
    {
      coordinates_ = other.coordinates_;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Move assignment.
     **********************************************************************************************/
    Point<space_dim>& operator= (Point<space_dim>&& other) noexcept
    {
      std::swap(coordinates_, other.coordinates_);
      return *this;
    }
    
    /*!*********************************************************************************************
     * \brief   Return reference to single coordinate of a point.
     * 
     * \param   coord_entry   An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    A reference to a \c pt_coord_t describing the coord_entry'th
     *                        coordinate.
     **********************************************************************************************/
    pt_coord_t& operator[](const unsigned int coord_entry)
    {
      hy_assert( 0 <= coord_entry && coord_entry < space_dim ,
                 "You can only access entries of a point's coordinates that have a non-negaitive "
                 << "index that is smaller than the space dimension (which is " << space_dim << ")."
                 << " However, you tried to access the " << coord_entry << "-th entry." );
      return coordinates_[coord_entry];
    }
    /*!*********************************************************************************************
     * \brief   Return single coordinate of a constant point.
     * 
     * \param   coord_entry   An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    \c pt_coord_t describing the coord_entry'th
     *                        coordinate.
     **********************************************************************************************/
    const pt_coord_t& operator[](const unsigned int coord_entry) const
    {
      hy_assert( 0 <= coord_entry && coord_entry < space_dim ,
                 "You can only access entries of a point's coordinates that have a non-negaitive "
                 << "index that is smaller than the space dimension (which is " << space_dim << ")."
                 << " However, you tried to access the " << coord_entry << "-th entry." );
      return coordinates_[coord_entry];
    }
    /*!*********************************************************************************************
     * \brief   Find out whether two points have (exactly) the same coordinates.
     * 
     * This function compares the point to another point and returns true if and only if both points
     * have exactly (that is not only with respect to some rounding errors) the same coordinates.
     * 
     * \param   other_point   Another \c Point<space_dim> that is to be dicriminate from.
     * \retval  isEqual       A \c boolean which is true if both points have the same coordinates.
     **********************************************************************************************/
    bool operator==(const Point<space_dim>& other_point) const
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        if (coordinates_[dim] != other_point[dim])  return false;
      return true;
    }
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
    bool operator!=(const Point<space_dim>& other_point) const
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        if (coordinates_[dim] != other_point[dim])  return true;
      return false;
    }
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
    bool operator<(const Point<space_dim>& other_point) const
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        if (coordinates_[dim] < other_point[dim])      return true;
        else if (coordinates_[dim] > other_point[dim]) return false;
      return false;
    }
    
    
    Point<space_dim>& operator*=(const pt_coord_t scale_fac)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        coordinates_[dim] *= scale_fac;
      return *this;
    }
    
    Point<space_dim>& operator/=(const pt_coord_t scale_denom)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        coordinates_[dim] /= scale_denom;
      return *this;
    }
    
    Point<space_dim>& operator+=(const pt_coord_t additum)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        coordinates_[dim] += additum;
      return *this;
    }

    Point<space_dim>& operator-=(const pt_coord_t subtractum)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        coordinates_[dim] -= subtractum;
      return *this;
    }
    
    Point<space_dim>& operator+=(const Point<space_dim>& other)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        coordinates_[dim] += other[dim];
      return *this;
    }
    
    Point<space_dim>& operator-=(const Point<space_dim>& other)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        coordinates_[dim] -= other[dim];
      return *this;
    }
    
    pt_coord_t operator*(const Point<space_dim>& other) const
    {
      pt_coord_t scalar_product = 0.;
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        scalar_product += coordinates_[dim] * other[dim];
      return scalar_product;
    }
    
}; // end of class Point

/*!*************************************************************************************************
 * \brief   Computes Euclidean norm of a \c Point.
 * 
 * \todo    Fill information, when details clarified with Guido.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim>
pt_coord_t norm_2(const Point<space_dim>& pt)
{
  return std::sqrt( pt * pt );
}

/*!*************************************************************************************************
 * \brief   Computes Euclidean distance between two \c Point.
 * 
 * \todo    Fill information, when details clarified with Guido.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim>
pt_coord_t distance_2(const Point<space_dim>& left, const Point<space_dim>& right)
{
  return norm_2( left - right );
}

/*!*************************************************************************************************
 * \brief   Fill \c stream with \c Point.
 * 
 * \todo    Fill information, when details clarified with Guido.
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim>
std::ostream& operator<< (std::ostream& stream, const Point<space_dim>& pt)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)
    stream << " " << pt[dim] << " ";
  return stream;
}

/*!*************************************************************************************************
 * \brief   Add two \c Point.
 * 
 * \todo    Discuss with Guido, why this way of implementation is to be preferred and why the assert
 *          in all related functions are never thrown (even if Point a + b + c is artificially added
 *          in executed code)! Compare the advice in
 *          https://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim>
Point<space_dim> operator+(const Point<space_dim>& left, const Point<space_dim>& right)
{
  Point<space_dim> sum(left);
  for (unsigned int dim = 0; dim < space_dim; ++dim)  sum[dim] += right[dim];
  return sum;
}

template<unsigned int space_dim>
Point<space_dim> operator+(Point<space_dim>&& left, const Point<space_dim>& right)
{
  hy_assert( 0 == 1 ,
             "This function is never called, although this way of implementation is recommended on "
             << "https://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter" );
  for (unsigned int dim = 0; dim < space_dim; ++dim)  left[dim] += right[dim];
  return left;
}

template<unsigned int space_dim>
Point<space_dim> operator+(const Point<space_dim>& left, Point<space_dim>&& right)
{
  hy_assert( 0 == 1 ,
             "This function is never called, although this way of implementation is recommended on "
             << "https://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter" );
  for (unsigned int dim = 0; dim < space_dim; ++dim)  right[dim] += left[dim];
  return right;
}

template<unsigned int space_dim>
Point<space_dim> operator+(Point<space_dim>&&left , Point<space_dim>&& right)
{
  hy_assert( 0 == 1 ,
             "This function is never called, although this way of implementation is recommended on "
             << "https://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter" );
  for (unsigned int dim = 0; dim < space_dim; ++dim)  left[dim] += right[dim];
  return left;
}

/*!*************************************************************************************************
 * \brief   Subtract two \c Point.
 * 
 * \todo    Discuss with Guido, why this way of implementation is to be preferred and why the assert
 *          in all related functions are never thrown (even if Point a - b - c is artificially added
 *          in executed code)! Compare the advice in
 *          https://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter
 * 
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim>
Point<space_dim> operator-(const Point<space_dim>& left, const Point<space_dim>& right)
{
  Point<space_dim> difference(left);
  for (unsigned int dim = 0; dim < space_dim; ++dim)  difference[dim] -= right[dim];
  return difference;
}

template<unsigned int space_dim>
Point<space_dim> operator-(Point<space_dim>&& left, const Point<space_dim>& right)
{
  hy_assert( 0 == 1 ,
             "This function is never called, although this way of implementation is recommended on "
             << "https://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter" );
  for (unsigned int dim = 0; dim < space_dim; ++dim)  left[dim] -= right[dim];
  return left;
}

template<unsigned int space_dim>
Point<space_dim> operator-(const Point<space_dim>& left, Point<space_dim>&& right)
{
  hy_assert( 0 == 1 ,
             "This function is never called, although this way of implementation is recommended on "
             << "https://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter" );
  for (unsigned int dim = 0; dim < space_dim; ++dim)  right[dim] = left[dim] - right[dim];
  return right;
}

template<unsigned int space_dim>
Point<space_dim> operator-(Point<space_dim>&&left , Point<space_dim>&& right)
{
  hy_assert( 0 == 1 ,
             "This function is never called, although this way of implementation is recommended on "
             << "https://stackoverflow.com/questions/11726171/numeric-vector-operator-overload-rvalue-reference-parameter" );
  for (unsigned int dim = 0; dim < space_dim; ++dim)  left[dim] -= right[dim];
  return left;
}

#endif // end of ifndef POINT_HXX
