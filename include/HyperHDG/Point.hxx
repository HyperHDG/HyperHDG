#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/HyAssert.hxx>

#include <array>
#include <cmath>
#include <ostream>

// #include <exception>

/*!*************************************************************************************************
 * \brief   Exception to be thrown if division by zero appears.
 *
 * \todo    Is this the way to do it intended by Guido? If so, should we include <exception>? It
 *          works perfecly without the include. I suspect that array (or another by that point
 *          included package includes exception?!
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019.
 * \authors   Andreas Rupp, Heidelberg University, 2019.
 **************************************************************************************************/
struct PointDivByZeroException : public std::exception
{
  const char * what () const throw () { return "Attempted division by zero."; }
};
/*!*************************************************************************************************
 * \brief   This class implements a point in a d-dimensional space.
 * 
 * This class implements a point in a \f$d\f$-dimensional space, where the \f$d\f$ is given by the
 * template parameter \c space_dim.
 * 
 * \tparam  pt_coord_t        Floating point type specification. Default is float.
 * \tparam  space_dim         The dimension of the space, the object is located in.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int space_dim, typename pt_coord_t = float >
class Point
{
  private:
    /*!*********************************************************************************************
     * \brief   Array containing the coordinates of the point.
     * 
     * A \c std::array conatining the i-th coordinate of the point as its i-th entry.
     **********************************************************************************************/
    std::array<pt_coord_t, space_dim> coordinates_;
  public:

    // Constructors and assignment operators:

    /*!*********************************************************************************************
     * \brief   Empty constructor for a point.
     * 
     * Fills coordinates of the point with zeros.
     **********************************************************************************************/
    Point() { coordinates_.fill(0.); }
    /*!*********************************************************************************************
     * \brief   Construct point from array of coordinates.
     * 
     * Fills the point's array of coordinates with the input parameter. 
     * 
     * \param   coordinates   A \c std::array containing the coordinates of the point.
     **********************************************************************************************/
    Point(const std::array<pt_coord_t, space_dim>& coordinates) : coordinates_(coordinates) { }
    /*!*********************************************************************************************
     * \brief   Copy constructor.
     **********************************************************************************************/
    Point(const Point<space_dim,pt_coord_t>& other) : coordinates_(other.coordinates_) { }
    /*!*********************************************************************************************
     * \brief   Move constructor.
     **********************************************************************************************/
    Point(Point<space_dim,pt_coord_t>&& other) noexcept
    : coordinates_(std::move(other.coordinates_)) { }
    /*!*********************************************************************************************
     * \brief   Copy assignment.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator= (const Point<space_dim,pt_coord_t>& other)
    { coordinates_ = other.coordinates_; return *this; }
    /*!*********************************************************************************************
     * \brief   Move assignment.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator= (Point<space_dim,pt_coord_t>&& other) noexcept
    { std::swap(coordinates_, other.coordinates_); return *this; }
    
    // Random access operators:

    /*!*********************************************************************************************
     * \brief   Return single coordinate of a constant point.
     * 
     * \param   coord_entry   An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    \c pt_coord_t describing the coord_entry'th coordinate.
     **********************************************************************************************/
    pt_coord_t operator[](const unsigned int coord_entry) const
    {
      hy_assert( 0 <= coord_entry && coord_entry < space_dim ,
                 "You can only access entries of a point's coordinates that have a non-negaitive "
                 << "index that is smaller than the space dimension (which is " << space_dim << ")."
                 << " However, you tried to access the " << coord_entry << "-th entry." );
      return coordinates_[coord_entry];
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

    // Comparison operators:

    /*!*********************************************************************************************
     * \brief   Find out whether two points have (exactly) the same coordinates.
     * 
     * This function compares the point to another point and returns true if and only if both points
     * have exactly (that is not only with respect to some rounding errors) the same coordinates.
     * 
     * \param   other_point   Another \c Point<space_dim> that is to be dicriminate from.
     * \retval  isEqual       A \c boolean which is true if both points have the same coordinates.
     **********************************************************************************************/
    bool operator==(const Point<space_dim,pt_coord_t>& other_point) const
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        if (coordinates_[dim] != other_point[dim])  return false;
      return true;
    }
    /*!*********************************************************************************************
     * \brief   Find out whether two points do not have (exactly) the same coordinates.
     * 
     * This function compares the point to another point and returns false if and only if both
     * points have exactly (that is not only with respect to some rounding errors) the same
     * coordinates.
     * 
     * \param   other_point   Another \c Point<space_dim> that is to be dicriminate from.
     * \retval  isEqual       A \c boolean which is false if both points have the same coordinates.
     **********************************************************************************************/
    bool operator!=(const Point<space_dim,pt_coord_t>& other_point) const
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
    bool operator<(const Point<space_dim,pt_coord_t>& other_point) const
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        if (coordinates_[dim] < other_point[dim])      return true;
        else if (coordinates_[dim] > other_point[dim]) return false;
      return false;
    }
    
    // Operators updating a Point by a scalar:

    /*!*********************************************************************************************
     * \brief   Add scalar to a given point.
     * 
     * \param   scalar        Floating point that is added to all of the point's coordinates.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator+=(const pt_coord_t scalar)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)  coordinates_[dim] += scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Subtract scalar from a given point.
     * 
     * \param   scalar        Floating point that is subtracted from all of the point's coordinates.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator-=(const pt_coord_t scalar)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)  coordinates_[dim] -= scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Multiply scalar a given point.
     * 
     * \param   scalar        Floating point that is multiplied with all of the point's coordinates.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator*=(const pt_coord_t scalar)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)  coordinates_[dim] *= scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Divide given point by a scalar.
     * 
     * \param   scalar        Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator/=(const pt_coord_t scalar)
    {
      if (scalar == 0.)  throw PointDivByZeroException();
//      hy_assert( scalar != 0. , "Division by zeros is not well-defined!" );
      for (unsigned int dim = 0; dim < space_dim; ++dim)  coordinates_[dim] /= scalar;
      return *this;
    }
    
    // Operators updating a Point by another Point:

    /*!*********************************************************************************************
     * \brief   Add point to given point.
     * 
     * \param   scalar        Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator+=(const Point<space_dim,pt_coord_t>& other)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)  coordinates_[dim] += other[dim];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Subtract other point from point.
     * 
     * \param   scalar        Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator-=(const Point<space_dim,pt_coord_t>& other)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)  coordinates_[dim] -= other[dim];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Hadamard product with other point.
     * 
     * \param   scalar        Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator*=(const Point<space_dim,pt_coord_t>& other)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)  coordinates_[dim] *= other[dim];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Hadamard division by other point.
     * 
     * \param   scalar        Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    Point<space_dim,pt_coord_t>& operator/=(const Point<space_dim,pt_coord_t>& other)
    {
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        if (other[dim] == 0.)  throw PointDivByZeroException();
//        hy_assert( other[dim] != 0. ,
//                   "Division by a point with " << dim << "-th component being 0 is prohibited!" );
        coordinates_[dim] /= other[dim];
      }
      return *this;
    }
    
    // Fundamental functions returning scalar from two Points:

    /*!*********************************************************************************************
     * \brief   Euclidean scalar product with other point.
     * 
     * \param   scalar        Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_point    The updated point.
     **********************************************************************************************/
    pt_coord_t operator*(const Point<space_dim,pt_coord_t>& other) const
    {
      pt_coord_t scalar_product = 0.;
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        scalar_product += coordinates_[dim] * other[dim];
      return scalar_product;
    }
}; // end of class Point

// Fundamental functions returning Point from two Points:

/*!*************************************************************************************************
 * \brief   Add two \c Point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator+
(const Point<space_dim,pt_coord_t>& left, const Point<space_dim,pt_coord_t>& right)
{
  Point<space_dim,pt_coord_t> sum(left);
  return sum += right;
}
/*!*************************************************************************************************
 * \brief   Subtract two \c Point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator-
(const Point<space_dim,pt_coord_t>& left, const Point<space_dim,pt_coord_t>& right)
{
  Point<space_dim,pt_coord_t> difference(left);
  return difference -= right;
}
/*!*************************************************************************************************
 * \brief   Hadamard product of two \c Point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> hada_prod
(const Point<space_dim,pt_coord_t>& left, const Point<space_dim,pt_coord_t>& right)
{
  Point<space_dim,pt_coord_t> product(left);
  return product *= right;
}
/*!*************************************************************************************************
 * \brief   Hadamard division two \c Point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template<unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> hada_divi
(const Point<space_dim,pt_coord_t>& left, const Point<space_dim,pt_coord_t>& right)
{
  Point<space_dim,pt_coord_t> quotient(left);
  return quotient /= right;;
}

// Fundamental functions returning Point from a scalar and a Point:

/*!*************************************************************************************************
 * \brief   Add point to scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator+
(const pt_coord_t& scalar, const Point<space_dim,pt_coord_t>& pt)
{
  Point<space_dim,pt_coord_t> sum(pt);
  return sum += scalar;
}
/*!*************************************************************************************************
 * \brief   Add scalar to point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator+
(const Point<space_dim,pt_coord_t>& pt, const pt_coord_t& scalar)
{
  Point<space_dim,pt_coord_t> sum(pt);
  return sum += scalar;
}
/*!*************************************************************************************************
 * \brief   Subtract point from scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator-
(const pt_coord_t& scalar, const Point<space_dim,pt_coord_t>& pt)
{
  Point<space_dim,pt_coord_t> difference(pt);
  for (unsigned int dim = 0; dim < space_dim; ++dim)  difference[dim] = scalar - pt[dim];
  return difference;
}
/*!*************************************************************************************************
 * \brief   Subtract scalar from point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator-
(const Point<space_dim,pt_coord_t>& pt, const pt_coord_t& scalar)
{
  Point<space_dim,pt_coord_t> difference(pt);
  return difference -= scalar;
}
/*!*************************************************************************************************
 * \brief   Multiply scalar with point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator*
(const pt_coord_t& scalar, const Point<space_dim,pt_coord_t>& pt)
{
  Point<space_dim,pt_coord_t> product(pt);
  return product *= scalar;
}
/*!*************************************************************************************************
 * \brief   Multiply point with scalar.
 *  
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator*
(const Point<space_dim,pt_coord_t>& pt, const pt_coord_t& scalar)
{
  Point<space_dim,pt_coord_t> product(pt);
  return product *= scalar;
}
/*!*************************************************************************************************
 * \brief   Divide scalar by \c Point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator/
(const pt_coord_t& scalar, const Point<space_dim,pt_coord_t>& pt)
{
  Point<space_dim,pt_coord_t> quotient(pt);
  for (unsigned int dim = 0; dim < space_dim; ++dim)
  {
    if (pt[dim] == 0.)  throw PointDivByZeroException();
//    hy_assert( pt[dim] != 0. ,
//               "Divison by a point whith " << dim << "-th component bein zero is not allowed!" );
    quotient[dim] = scalar / pt[dim];
  }
  return quotient;
}
/*!*************************************************************************************************
 * \brief   Divide \c Point by scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
Point<space_dim,pt_coord_t> operator/
(const Point<space_dim,pt_coord_t>& pt, const pt_coord_t& scalar)
{
  Point<space_dim,pt_coord_t> quotient(pt);
  return quotient /= scalar;
}

// Norms of Points:

/*!*************************************************************************************************
 * \brief   Absolute sum norm of point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
pt_coord_t norm_1(const Point<space_dim,pt_coord_t>& pt)
{
  pt_coord_t norm = 0.;
  for (unsigned int dim = 0; dim < space_dim; ++dim)  norm += std::abs( pt[dim] );
  return norm;
}
/*!*************************************************************************************************
 * \brief   Computes Euclidean norm of a point.
 * 
 * \todo    Fill information, when details clarified with Guido.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
pt_coord_t norm_2(const Point<space_dim,pt_coord_t>& pt)
{
  return std::sqrt( pt * pt );
}
/*!*************************************************************************************************
 * \brief   Maximum norm of point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
pt_coord_t norm_infty(const Point<space_dim,pt_coord_t>& pt)
{
  pt_coord_t norm = std::abs( pt[0] );
  for (unsigned int dim = 1; dim < space_dim; ++dim)  norm = std::max( norm, std::abs(pt[dim]) );
  return norm;
}
/*!*************************************************************************************************
 * \brief   p-norm of point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
pt_coord_t norm_p(const Point<space_dim,pt_coord_t>& pt, const pt_coord_t power)
{
  pt_coord_t norm = 0.;
  for (unsigned int dim = 0; dim < space_dim; ++dim)  norm += std::pow( std::abs(pt[dim]) , power );
  return std::pow( norm , 1. / power );
}

// Output of Point:

/*!*************************************************************************************************
 * \brief   Fill \c stream with \c Point.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int space_dim, typename pt_coord_t >
std::ostream& operator<< (std::ostream& stream, const Point<space_dim,pt_coord_t>& pt)
{
  for (unsigned int dim = 0; dim < space_dim; ++dim)  stream << " " << pt[dim] << " ";
  return stream;
}
