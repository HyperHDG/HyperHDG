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
 *          works perfecly without the include. I suspect that array (or another by that SmallVec
 *          included package includes exception?!
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019.
 * \authors   Andreas Rupp, Heidelberg University, 2019.
 **************************************************************************************************/
struct SmallVecDivByZeroException : public std::exception
{
  const char * what () const throw () { return "Attempted division by zero."; }
};

/*!*************************************************************************************************
 * \brief   This class implements a SmallVec in a d-dimensional space.
 * 
 * This class implements a SmallVec in a \f$d\f$-dimensional space, where the \f$d\f$ is given by
 * the template parameter \c vec_dim.
 * 
 * \tparam  vec_entry_t         Floating point type specification. Default is double.
 * \tparam  vec_dim             The dimension of the space, the object is located in.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int vec_dim, typename vec_entry_t = double >
class SmallVec
{
  private:
    /*!*********************************************************************************************
     * \brief   Array containing the coordinates of the SmallVec.
     * 
     * A \c std::array conatining the i-th coordinate of the SmallVec as its i-th entry.
     **********************************************************************************************/
    std::array<vec_entry_t, vec_dim> coordinates_;
  public:

    // Constructors and assignment operators:

    /*!*********************************************************************************************
     * \brief   Empty constructor for a SmallVec.
     * 
     * Fills coordinates of the SmallVec with zeros.
     **********************************************************************************************/
    SmallVec() { coordinates_.fill(0.); }
    /*!*********************************************************************************************
     * \brief   Construct SmallVec from array of coordinates.
     * 
     * Fills the SmallVec's array of coordinates with the input parameter. 
     * 
     * \param   coordinates   A \c std::array containing the coordinates of the SmallVec.
     **********************************************************************************************/
    SmallVec(const std::array<vec_entry_t, vec_dim>& coordinates) : coordinates_(coordinates) { }
    /*!*********************************************************************************************
     * \brief   Copy constructor.
     **********************************************************************************************/
    SmallVec(const SmallVec<vec_dim,vec_entry_t>& other) : coordinates_(other.coordinates_) { }
    /*!*********************************************************************************************
     * \brief   Move constructor.
     **********************************************************************************************/
    SmallVec(SmallVec<vec_dim,vec_entry_t>&& other) noexcept
    : coordinates_(std::move(other.coordinates_)) { }
    /*!*********************************************************************************************
     * \brief   Copy assignment.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator= (const SmallVec<vec_dim,vec_entry_t>& other)
    { coordinates_ = other.coordinates_; return *this; }
    /*!*********************************************************************************************
     * \brief   Move assignment.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator= (SmallVec<vec_dim,vec_entry_t>&& other) noexcept
    { std::swap(coordinates_, other.coordinates_); return *this; }
    
    // Random access operators:

    /*!*********************************************************************************************
     * \brief   Return single coordinate of a constant SmallVec.
     * 
     * \param   coord_entry     An \c unsigned \c int referring to the coordinate that is to be
     *                          returned.
     * \retval  coordinate      \c vec_entry_t describing the coord_entry'th coordinate.
     **********************************************************************************************/
    vec_entry_t operator[](const unsigned int coord_entry) const
    {
      hy_assert( 0 <= coord_entry && coord_entry < vec_dim ,
                 "You can only access entries of a SmallVec's coordinates that have non-negaitive "
                 << "index that is smaller than the space dimension (which is " << vec_dim << ")."
                 << " However, you tried to access the " << coord_entry << "-th entry." );
      return coordinates_[coord_entry];
    }
    /*!*********************************************************************************************
     * \brief   Return reference to single coordinate of a SmallVec.
     * 
     * \param   coord_entry     An \c unsigned \c int referring to the coordinate that is to be
     *                          returned.
     * \retval  coordinate      A reference to a \c vec_entry_t describing the coord_entry'th
     *                          coordinate.
     **********************************************************************************************/
    vec_entry_t& operator[](const unsigned int coord_entry)
    {
      hy_assert( 0 <= coord_entry && coord_entry < vec_dim ,
                 "You can only access entries of a SmallVec's coordinates that have non-negaitive "
                 << "index that is smaller than the space dimension (which is " << vec_dim << ")."
                 << " However, you tried to access the " << coord_entry << "-th entry." );
      return coordinates_[coord_entry];
    }

    // Comparison operators:

    /*!*********************************************************************************************
     * \brief   Find out whether two SmallVecs have (exactly) the same coordinates.
     * 
     * This function compares the SmallVec to another SmallVec and returns true if and only if both
     * SmallVecs have exactly (that is not only with respect to some rounding errors) the same
     * coordinates.
     * 
     * \param   other_SmallVec  Another \c SmallVec<vec_dim> that is to be dicriminate from.
     * \retval  isEqual         A \c boolean which is true if both SmallVecs have the same coords.
     **********************************************************************************************/
    bool operator==(const SmallVec<vec_dim,vec_entry_t>& other_SmallVec) const
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)
        if (coordinates_[dim] != other_SmallVec[dim])  return false;
      return true;
    }
    /*!*********************************************************************************************
     * \brief   Find out whether two SmallVecs do not have (exactly) the same coordinates.
     * 
     * This function compares the SmallVec to another SmallVec and returns false if and only if both
     * SmallVecs have exactly (that is not only with respect to some rounding errors) the same
     * coordinates.
     * 
     * \param   other_SmallVec  Another \c SmallVec<vec_dim> that is to be dicriminate from.
     * \retval  isEqual         A \c boolean which is false if both SmallVecs have the same coords.
     **********************************************************************************************/
    bool operator!=(const SmallVec<vec_dim,vec_entry_t>& other_SmallVec) const
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)
        if (coordinates_[dim] != other_SmallVec[dim])  return true;
      return false;
    }
    /*!*********************************************************************************************
     * \brief   Find out whether the SmallVec is "smaller than" another SmallVec.
     * 
     * This function compares the SmallVec to another SmallVec and returns true if and only if the
     * lowest ranked coordinate (according to the coordinate index) where the both SmallVecs are not
     * equal of the given SmallVec is smaller than that of the other SmallVec. It is false, if both
     * SmallVecs are equal.
     * 
     * \param   other_SmallVec  Another \c SmallVec<vec_dim> that is to be dicriminate from.
     * \retval  isEqual         A \c boolean which is true if the left SmallVec is strictly smaller
     *                          than the right one.
     **********************************************************************************************/
    bool operator<(const SmallVec<vec_dim,vec_entry_t>& other_SmallVec) const
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)
        if (coordinates_[dim] < other_SmallVec[dim])      return true;
        else if (coordinates_[dim] > other_SmallVec[dim]) return false;
      return false;
    }
    
    // Operators updating a SmallVec by a scalar:

    /*!*********************************************************************************************
     * \brief   Add scalar to a given SmallVec.
     * 
     * \param   scalar          Floating point that is added to all of the SmallVec's coordinates.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator+=(const vec_entry_t scalar)
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)  coordinates_[dim] += scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Subtract scalar from a given SmallVec.
     * 
     * \param   scalar          Floating point that is subtracted from all of the SmallVec's
     *                          coordinates.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator-=(const vec_entry_t scalar)
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)  coordinates_[dim] -= scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Multiply scalar a given SmallVec.
     * 
     * \param   scalar          Floating point that is multiplied with all of the SmallVec's
     *                          coordinates.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator*=(const vec_entry_t scalar)
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)  coordinates_[dim] *= scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Divide given SmallVec by a scalar.
     * 
     * \param   scalar          Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator/=(const vec_entry_t scalar)
    {
      if (scalar == 0.)  throw SmallVecDivByZeroException();
      for (unsigned int dim = 0; dim < vec_dim; ++dim)  coordinates_[dim] /= scalar;
      return *this;
    }
    
    // Operators updating a SmallVec by another SmallVec:

    /*!*********************************************************************************************
     * \brief   Add SmallVec to given SmallVec.
     * 
     * \param   scalar          Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator+=(const SmallVec<vec_dim,vec_entry_t>& other)
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)  coordinates_[dim] += other[dim];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Subtract other SmallVec from SmallVec.
     * 
     * \param   scalar          Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator-=(const SmallVec<vec_dim,vec_entry_t>& other)
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)  coordinates_[dim] -= other[dim];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Hadamard product with other SmallVec.
     * 
     * \param   scalar          Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator*=(const SmallVec<vec_dim,vec_entry_t>& other)
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)  coordinates_[dim] *= other[dim];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Hadamard division by other SmallVec.
     * 
     * \param   scalar          Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    SmallVec<vec_dim,vec_entry_t>& operator/=(const SmallVec<vec_dim,vec_entry_t>& other)
    {
      for (unsigned int dim = 0; dim < vec_dim; ++dim)
      {
        if (other[dim] == 0.)  throw SmallVecDivByZeroException();
        coordinates_[dim] /= other[dim];
      }
      return *this;
    }
    
    // Fundamental functions returning scalar from two SmallVecs:

    /*!*********************************************************************************************
     * \brief   Euclidean scalar product with other SmallVec.
     * 
     * \param   scalar          Floating point (\f$\neq 0\f$) all coordinates are divided by.
     * \retval  this_SmallVec   The updated SmallVec.
     **********************************************************************************************/
    vec_entry_t operator*(const SmallVec<vec_dim,vec_entry_t>& other) const
    {
      vec_entry_t scalar_product = 0.;
      for (unsigned int dim = 0; dim < vec_dim; ++dim)
        scalar_product += coordinates_[dim] * other[dim];
      return scalar_product;
    }
}; // end of class SmallVec

// Fundamental functions returning SmallVec from two SmallVecs:

/*!*************************************************************************************************
 * \brief   Add two \c SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template<unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator+
(const SmallVec<vec_dim,vec_entry_t>& left, const SmallVec<vec_dim,vec_entry_t>& right)
{
  SmallVec<vec_dim,vec_entry_t> sum(left);
  return sum += right;
}
/*!*************************************************************************************************
 * \brief   Subtract two \c SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template<unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator-
(const SmallVec<vec_dim,vec_entry_t>& left, const SmallVec<vec_dim,vec_entry_t>& right)
{
  SmallVec<vec_dim,vec_entry_t> difference(left);
  return difference -= right;
}
/*!*************************************************************************************************
 * \brief   Hadamard product of two \c SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template<unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> hada_prod
(const SmallVec<vec_dim,vec_entry_t>& left, const SmallVec<vec_dim,vec_entry_t>& right)
{
  SmallVec<vec_dim,vec_entry_t> product(left);
  return product *= right;
}
/*!*************************************************************************************************
 * \brief   Hadamard division two \c SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template<unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> hada_divi
(const SmallVec<vec_dim,vec_entry_t>& left, const SmallVec<vec_dim,vec_entry_t>& right)
{
  SmallVec<vec_dim,vec_entry_t> quotient(left);
  return quotient /= right;;
}

// Fundamental functions returning SmallVec from a scalar and a SmallVec:

/*!*************************************************************************************************
 * \brief   Add SmallVec to scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator+
(const vec_entry_t& scalar, const SmallVec<vec_dim,vec_entry_t>& pt)
{
  SmallVec<vec_dim,vec_entry_t> sum(pt);
  return sum += scalar;
}
/*!*************************************************************************************************
 * \brief   Add scalar to SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator+
(const SmallVec<vec_dim,vec_entry_t>& pt, const vec_entry_t& scalar)
{
  SmallVec<vec_dim,vec_entry_t> sum(pt);
  return sum += scalar;
}
/*!*************************************************************************************************
 * \brief   Subtract SmallVec from scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator-
(const vec_entry_t& scalar, const SmallVec<vec_dim,vec_entry_t>& pt)
{
  SmallVec<vec_dim,vec_entry_t> difference(pt);
  for (unsigned int dim = 0; dim < vec_dim; ++dim)  difference[dim] = scalar - pt[dim];
  return difference;
}
/*!*************************************************************************************************
 * \brief   Subtract scalar from SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator-
(const SmallVec<vec_dim,vec_entry_t>& pt, const vec_entry_t& scalar)
{
  SmallVec<vec_dim,vec_entry_t> difference(pt);
  return difference -= scalar;
}
/*!*************************************************************************************************
 * \brief   Multiply scalar with SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator*
(const vec_entry_t& scalar, const SmallVec<vec_dim,vec_entry_t>& pt)
{
  SmallVec<vec_dim,vec_entry_t> product(pt);
  return product *= scalar;
}
/*!*************************************************************************************************
 * \brief   Multiply SmallVec with scalar.
 *  
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator*
(const SmallVec<vec_dim,vec_entry_t>& pt, const vec_entry_t& scalar)
{
  SmallVec<vec_dim,vec_entry_t> product(pt);
  return product *= scalar;
}
/*!*************************************************************************************************
 * \brief   Divide scalar by \c SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator/
(const vec_entry_t& scalar, const SmallVec<vec_dim,vec_entry_t>& pt)
{
  SmallVec<vec_dim,vec_entry_t> quotient(pt);
  for (unsigned int dim = 0; dim < vec_dim; ++dim)
  {
    if (pt[dim] == 0.)  throw SmallVecDivByZeroException();
    quotient[dim] = scalar / pt[dim];
  }
  return quotient;
}
/*!*************************************************************************************************
 * \brief   Divide \c SmallVec by scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
SmallVec<vec_dim,vec_entry_t> operator/
(const SmallVec<vec_dim,vec_entry_t>& pt, const vec_entry_t& scalar)
{
  SmallVec<vec_dim,vec_entry_t> quotient(pt);
  return quotient /= scalar;
}

// Norms of SmallVecs:

/*!*************************************************************************************************
 * \brief   Absolute sum norm of SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
vec_entry_t norm_1(const SmallVec<vec_dim,vec_entry_t>& pt)
{
  vec_entry_t norm = 0.;
  for (unsigned int dim = 0; dim < vec_dim; ++dim)  norm += std::abs( pt[dim] );
  return norm;
}
/*!*************************************************************************************************
 * \brief   Computes Euclidean norm of a SmallVec.
 * 
 * \todo    Fill information, when details clarified with Guido.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
vec_entry_t norm_2(const SmallVec<vec_dim,vec_entry_t>& pt)
{
  return std::sqrt( pt * pt );
}
/*!*************************************************************************************************
 * \brief   Maximum norm of SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
vec_entry_t norm_infty(const SmallVec<vec_dim,vec_entry_t>& pt)
{
  vec_entry_t norm = std::abs( pt[0] );
  for (unsigned int dim = 1; dim < vec_dim; ++dim)  norm = std::max( norm, std::abs(pt[dim]) );
  return norm;
}
/*!*************************************************************************************************
 * \brief   p-norm of SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
vec_entry_t norm_p(const SmallVec<vec_dim,vec_entry_t>& pt, const vec_entry_t power)
{
  vec_entry_t norm = 0.;
  for (unsigned int dim = 0; dim < vec_dim; ++dim)  norm += std::pow( std::abs(pt[dim]) , power );
  return std::pow( norm , 1. / power );
}

// Output of SmallVec:

/*!*************************************************************************************************
 * \brief   Fill \c stream with \c SmallVec.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int vec_dim, typename vec_entry_t >
std::ostream& operator<< (std::ostream& stream, const SmallVec<vec_dim,vec_entry_t>& pt)
{
  for (unsigned int dim = 0; dim < vec_dim; ++dim)  stream << " " << pt[dim] << " ";
  return stream;
}


// Derived classes:

/*!*************************************************************************************************
 * \brief   A Point is a SmallVec, where the standard vec_entry_t is float.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int space_dim, typename pt_coord_t = float >
using Point = SmallVec<space_dim, pt_coord_t>;