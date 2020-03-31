#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/HyAssert.hxx>
#include <HyperHDG/LapackWrapper.hxx>

#include <array>
#include <cmath>
#include <ostream>

// #include <exception>

/*!*************************************************************************************************
 * \brief   Exception to be thrown if division by zero appears.
 *
 * \todo    Is this the way to do it intended by Guido? If so, should we include <exception>? It
 *          works perfecly without the include. I suspect that array (or another by that SmallMat
 *          included package includes exception?!
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019.
 * \authors   Andreas Rupp, Heidelberg University, 2019.
 **************************************************************************************************/
struct SmallMatDivByZeroException : public std::exception
{
  const char * what () const throw () { return "Attempted division by zero."; }
};

/*!*************************************************************************************************
 * \brief   This class implements a small matrix..
 * 
 * This class implements a SmallMat in a \f$d\f$-dimensional space, where the \f$d\f$ is given by
 * the template parameter \c n_rows,n_cols.
 * 
 * \tparam  mat_entry_t        Floating point type specification. Default is double.
 * \tparam  n_rows,n_cols         The dimension of the space, the object is located in.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rows, unsigned int n_cols = n_rows, typename mat_entry_t = double >
class SmallMat
{
  public:
    /*!*********************************************************************************************
     * \brief   Return size a SmallMat.
     * 
     * \retval  size            Amount of entries that can be stored within the matrix
     **********************************************************************************************/
    static constexpr unsigned int size()  { return n_rows * n_cols; }
    /*!*********************************************************************************************
     * \brief   Return dimensions a SmallMat.
     * 
     * \retval  dimensions      Number of rows and number of columns of the matrix.
     **********************************************************************************************/
    static constexpr std::array<unsigned int, 2> dimensions()
    { return std::array<unsigned int, 2> {n_rows, n_cols}; }
  private:
    /*!*********************************************************************************************
     * \brief   Array containing the entries of the SmallMat.
     * 
     * A \c std::array conatining the i-th coordinate of the SmallMat as its i-th entry.
     **********************************************************************************************/
    std::array<mat_entry_t, size()> entries_;
    /*!*********************************************************************************************
     * \brief   Translate row and column indices to local index of entry in matrix.
     * 
     * Local \f$ n \times n \f$ matrices are encoded as arrays of size \f$n^2\f$. This function
     * translates a row and a column index into the index of the long array, where the corresponding
     * entry is located. Note that this is done column-wise (not row-wise as usually), to have the
     * correct format for LAPACK.
     *
     * The function is static inline, since it is used in the constructor's initializer list.
     *
     * \param   row           Row index of local mtatrix entry.
     * \param   column        Column index of local matrix entry.
     * \retval  index         Overall index of local matrix entry.
     **********************************************************************************************/
    static inline unsigned int loc_matrix_index(const unsigned int row, const unsigned int column)
    {
      hy_assert( 0 <= row && row < n_rows ,
                 "Row index must be >= 0 and smaller than number of rows." );
      hy_assert( 0 <= column && column < n_cols ,
                 "Column index must be >= 0 and smaller than number of columns." );

      return column * n_rows + row;  // Encoded like this for easy use of LAPACK!
    }
  public:

    // Constructors and assignment operators:

    /*!*********************************************************************************************
     * \brief   Empty constructor for a SmallMat.
     * 
     * Fills entries of the SmallMat with zeros.
     **********************************************************************************************/
    SmallMat() { entries_.fill(0.); }
    /*!*********************************************************************************************
     * \brief   Construct SmallMat from array of entries.
     * 
     * Fills the SmallMat's array of entries with the input parameter. 
     * 
     * \param   entries   A \c std::array containing the entries of the SmallMat.
     **********************************************************************************************/
    SmallMat(const std::array<mat_entry_t, size()>& entries) : entries_(entries) { }
    /*!*********************************************************************************************
     * \brief   Move constructor from array.
     **********************************************************************************************/
    SmallMat(std::array<mat_entry_t, size()>&& entries) noexcept
    : entries_(std::move(entries)) { }
    /*!*********************************************************************************************
     * \brief   Copy constructor.
     **********************************************************************************************/
    SmallMat(const SmallMat<n_rows,n_cols,mat_entry_t>& other) : entries_(other.entries_) { }
    /*!*********************************************************************************************
     * \brief   Move constructor.
     **********************************************************************************************/
    SmallMat(SmallMat<n_rows,n_cols,mat_entry_t>&& other) noexcept
    : entries_(std::move(other.entries_)) { }
    /*!*********************************************************************************************
     * \brief   Copy assignment.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator=
    ( const SmallMat<n_rows,n_cols,mat_entry_t>& other )
    { entries_ = other.entries_; return *this; }
    /*!*********************************************************************************************
     * \brief   Move assignment.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator=
    ( SmallMat<n_rows,n_cols,mat_entry_t>&& other ) noexcept
    { std::swap(entries_, other.entries_); return *this; }
    
    // Return array with data:

    /*!*********************************************************************************************
     * \brief   Return data array.
     **********************************************************************************************/
    std::array<mat_entry_t, size()>& data()  { return entries_; }
    /*!*********************************************************************************************
     * \brief   Return data array.
     **********************************************************************************************/
    const std::array<mat_entry_t, size()>& data() const  { return entries_; }

    // Random access operators:

    /*!*********************************************************************************************
     * \brief   Return single entry of a constant SmallMat.
     * 
     * \param   index         An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    \c mat_entry_t describing the coord_entry'th coordinate.
     **********************************************************************************************/
    mat_entry_t operator()(const unsigned int row, const unsigned int column) const
    { return operator[](loc_matrix_index(row, column)); }
    /*!*********************************************************************************************
     * \brief   Return reference to single entry of a SmallMat.
     * 
     * \param   index         An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    A reference to a \c mat_entry_t describing the coord_entry'th
     *                        coordinate.
     **********************************************************************************************/
    mat_entry_t& operator()(const unsigned int row, const unsigned int column)
    { return operator[](loc_matrix_index(row, column)); }

    /*!*********************************************************************************************
     * \brief   Return single entry of a constant SmallMat.
     * 
     * \param   index         An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    \c mat_entry_t describing the coord_entry'th coordinate.
     **********************************************************************************************/
    mat_entry_t operator[](const unsigned int index) const
    {
      hy_assert( 0 <= index && index < size() ,
                 "You can only access entries of a SmallMat's entries that have non-negaitive "
                 << "index that is smaller than its size (which is " << size() << ")."
                 << " However, you tried to access the " << index << "-th entry." );
      return entries_[index];
    }
    /*!*********************************************************************************************
     * \brief   Return reference to single entry of a SmallMat.
     * 
     * \param   index         An \c unsigned \c int referring to the coordinate that is to be
     *                        returned.
     * \retval  coordinate    A reference to a \c mat_entry_t describing the coord_entry'th
     *                        coordinate.
     **********************************************************************************************/
    mat_entry_t& operator[](const unsigned int index)
    {
      hy_assert( 0 <= index && index < size() ,
                 "You can only access entries of a SmallMat's entries that have non-negaitive "
                 << "index that is smaller than its size (which is " << size() << ")."
                 << " However, you tried to access the " << index << "-th entry." );
      return entries_[index];
    }

    // Comparison operators:

    /*!*********************************************************************************************
     * \brief   Find out whether two SmallMats have (exactly) the same entries.
     * 
     * This function compares the SmallMat to another SmallMat and returns true if and only if both
     * SmallMats have exactly (that is not only with respect to some rounding errors) the same
     * entries.
     * 
     * \param   other_SmallMat  Another \c SmallMat<n_rows,n_cols> that is to be dicriminate from.
     * \retval  isEqual         A \c boolean which is true if both SmallMats have the same coords.
     **********************************************************************************************/
    bool operator==(const SmallMat<n_rows,n_cols,mat_entry_t>& other_SmallMat) const
    {
      for (unsigned int index = 0; index < size(); ++index)
        if (entries_[index] != other_SmallMat[index])  return false;
      return true;
    }
    /*!*********************************************************************************************
     * \brief   Find out whether two SmallMats do not have (exactly) the same entries.
     * 
     * This function compares the SmallMat to another SmallMat and returns false if and only if both
     * SmallMats have exactly (that is not only with respect to some rounding errors) the same
     * entries.
     * 
     * \param   other_SmallMat  Another \c SmallMat<n_rows,n_cols> that is to be dicriminate from.
     * \retval  isEqual         A \c boolean which is false if both SmallMats have the same coords.
     **********************************************************************************************/
    bool operator!=(const SmallMat<n_rows,n_cols,mat_entry_t>& other_SmallMat) const
    {
      for (unsigned int index = 0; index < size(); ++index)
        if (entries_[index] != other_SmallMat[index])  return true;
      return false;
    }
    /*!*********************************************************************************************
     * \brief   Find out whether the SmallMat is "smaller than" another SmallMat.
     * 
     * This function compares the SmallMat to another SmallMat and returns true if and only if the
     * lowest ranked coordinate (according to the coordinate index) where the both SmallMats are not
     * equal of the given SmallMat is smaller than that of the other SmallMat. It is false, if both
     * SmallMats are equal.
     * 
     * \param   other_SmallMat  Another \c SmallMat<n_rows,n_cols> that is to be dicriminate from.
     * \retval  isEqual         A \c boolean which is true if the left SmallMat is strictly smaller
     *                          than the right one.
     **********************************************************************************************/
    bool operator<(const SmallMat<n_rows,n_cols,mat_entry_t>& other_SmallMat) const
    {
      for (unsigned int index = 0; index < size(); ++index)
        if (entries_[index] < other_SmallMat[index])      return true;
        else if (entries_[index] > other_SmallMat[index]) return false;
      return false;
    }
    
    // Operators updating a SmallMat by a scalar:

    /*!*********************************************************************************************
     * \brief   Add scalar to a given SmallMat.
     * 
     * \param   scalar        Floating SmallMat that is added to all of the SmallMat's entries.
     * \retval  this_SmallMat    The updated SmallMat.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator+=(const mat_entry_t scalar)
    {
      for (unsigned int index = 0; index < size(); ++index)  entries_[index] += scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Subtract scalar from a given SmallMat.
     * 
     * \param   scalar           Floating point that is subtracted from all of the SmallMat's
     *                           entries.
     * \retval  this_SmallMat    The updated SmallMat.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator-=(const mat_entry_t scalar)
    {
      for (unsigned int index = 0; index < size(); ++index)  entries_[index] -= scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Multiply scalar a given SmallMat.
     * 
     * \param   scalar           Floating point that is multiplied with all of the SmallMat's
     *                           entries.
     * \retval  this_SmallMat    The updated SmallMat.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator*=(const mat_entry_t scalar)
    {
      for (unsigned int index = 0; index < size(); ++index)  entries_[index] *= scalar;
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Divide given SmallMat by a scalar.
     * 
     * \param   scalar        Floating SmallMat (\f$\neq 0\f$) all entries are divided by.
     * \retval  this_SmallMat    The updated SmallMat.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator/=(const mat_entry_t scalar)
    {
      if (scalar == 0.)  throw SmallMatDivByZeroException();
      for (unsigned int index = 0; index < size(); ++index)  entries_[index] /= scalar;
      return *this;
    }
    
    // Operators updating a SmallMat by another SmallMat:

    /*!*********************************************************************************************
     * \brief   Add SmallMat to given SmallMat.
     * 
     * \param   scalar        Floating SmallMat (\f$\neq 0\f$) all entries are divided by.
     * \retval  this_SmallMat    The updated SmallMat.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator+=
    ( const SmallMat<n_rows,n_cols,mat_entry_t>& other )
    {
      for (unsigned int index = 0; index < size(); ++index)  entries_[index] += other[index];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Subtract other SmallMat from SmallMat.
     * 
     * \param   scalar        Floating SmallMat (\f$\neq 0\f$) all entries are divided by.
     * \retval  this_SmallMat    The updated SmallMat.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator-=
    ( const SmallMat<n_rows,n_cols,mat_entry_t>& other )
    {
      for (unsigned int index = 0; index < size(); ++index)  entries_[index] -= other[index];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Hadamard product with other SmallMat.
     * 
     * \param   scalar        Floating SmallMat (\f$\neq 0\f$) all entries are divided by.
     * \retval  this_SmallMat    The updated SmallMat.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator*=
    ( const SmallMat<n_rows,n_cols,mat_entry_t>& other )
    {
      for (unsigned int index = 0; index < size(); ++index)  entries_[index] *= other[index];
      return *this;
    }
    /*!*********************************************************************************************
     * \brief   Hadamard division by other SmallMat.
     * 
     * \param   scalar        Floating SmallMat (\f$\neq 0\f$) all entries are divided by.
     * \retval  this_SmallMat    The updated SmallMat.
     **********************************************************************************************/
    SmallMat<n_rows,n_cols,mat_entry_t>& operator/=
    ( const SmallMat<n_rows,n_cols,mat_entry_t>& other )
    {
      for (unsigned int index = 0; index < size(); ++index)
      {
        if (other[index] == 0.)  throw SmallMatDivByZeroException();
        entries_[index] /= other[index];
      }
      return *this;
    }
}; // end of class SmallMat

 // Fundamental functions returning scalar from two SmallVecs:

/*!*************************************************************************************************
 * \brief   Euclidean scalar product with other SmallVec.
 * 
 * \param   scalar          Floating point (\f$\neq 0\f$) all entries are divided by.
 * \retval  this_SmallVec   The updated SmallVec.
 **************************************************************************************************/
template < unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
mat_entry_t scalar_product
(const SmallMat<n_rows,n_cols,mat_entry_t>& left, const SmallMat<n_rows,n_cols,mat_entry_t>& right)
{
  mat_entry_t scalar_product = 0.;
  for (unsigned int index = 0; index < left.size(); ++index)
    scalar_product += left[index] * right[index];
  return scalar_product;
}

// Fundamental functions returning SmallMat from two SmallMats:

/*!*************************************************************************************************
 * \brief   Add two \c SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator+
(const SmallMat<n_rows,n_cols,mat_entry_t>& left, const SmallMat<n_rows,n_cols,mat_entry_t>& right)
{
  SmallMat<n_rows,n_cols,mat_entry_t> sum(left);
  return sum += right;
}
/*!*************************************************************************************************
 * \brief   Subtract two \c SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator-
(const SmallMat<n_rows,n_cols,mat_entry_t>& left, const SmallMat<n_rows,n_cols,mat_entry_t>& right)
{
  SmallMat<n_rows,n_cols,mat_entry_t> difference(left);
  return difference -= right;
}
/*!*************************************************************************************************
 * \brief   Hadamard product of two \c SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> hada_prod
(const SmallMat<n_rows,n_cols,mat_entry_t>& left, const SmallMat<n_rows,n_cols,mat_entry_t>& right)
{
  SmallMat<n_rows,n_cols,mat_entry_t> product(left);
  return product *= right;
}
/*!*************************************************************************************************
 * \brief   Hadamard division two \c SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> hada_divi
(const SmallMat<n_rows,n_cols,mat_entry_t>& left, const SmallMat<n_rows,n_cols,mat_entry_t>& right)
{
  SmallMat<n_rows,n_cols,mat_entry_t> quotient(left);
  return quotient /= right;
}

// Standard matrix matrix multiplication and sove linear systems of equations:

/*!*************************************************************************************************
 * \brief   Standard matrix vector multiplication.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rowsA, unsigned int n_colsA, unsigned int n_colsB, typename mat_entry_t >
SmallMat<n_rowsA,n_colsB,mat_entry_t> operator*
( const SmallMat<n_rowsA,n_colsA,mat_entry_t>& A, const SmallMat<n_colsA,n_colsB,mat_entry_t>& B )
{
  SmallMat<n_rowsA,n_colsB,mat_entry_t> result;
  for (unsigned int colB = 0; colB < n_colsB; ++colB)
    for (unsigned int colA = 0; colA < n_colsA; ++colA)
      for (unsigned int rowA = 0; rowA < n_rowsA; ++rowA)
        result(rowA,colB) += A(rowA,colA) * B(colA,colB);
  return result;
}
/*!*************************************************************************************************
 * \brief   Standard matrix vector multiplication.
 * 
 * \todo    Implement this.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rowsA, unsigned int n_colsB, typename mat_entry_t >
SmallMat<n_rowsA,n_colsB,mat_entry_t> operator/
( SmallMat<n_rowsA,n_colsB,mat_entry_t>& b, SmallMat<n_rowsA,n_rowsA,mat_entry_t>& A )
{ return lapack_solve<n_rowsA,n_colsB>(A.data(), b.data()); }
/*!*************************************************************************************************
 * \brief   Standard matrix vector multiplication.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rowsA, unsigned int n_colsB, typename mat_entry_t >
SmallMat<n_rowsA,n_colsB,mat_entry_t> operator/
( const SmallMat<n_rowsA,n_colsB,mat_entry_t>& b, const SmallMat<n_rowsA,n_rowsA,mat_entry_t>& A )
{
  SmallMat<n_rowsA,n_colsB,mat_entry_t> helperb(b);
  SmallMat<n_rowsA,n_rowsA,mat_entry_t> helperA(A);
  return helperb / helperA;
}

// Fundamental functions returning SmallMat from a scalar and a SmallMat:

/*!*************************************************************************************************
 * \brief   Add SmallMat to scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator+
(const mat_entry_t& scalar, const SmallMat<n_rows,n_cols,mat_entry_t>& mat)
{
  SmallMat<n_rows,n_cols,mat_entry_t> sum(mat);
  return sum += scalar;
}
/*!*************************************************************************************************
 * \brief   Add scalar to SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator+
(const SmallMat<n_rows,n_cols,mat_entry_t>& mat, const mat_entry_t& scalar)
{
  SmallMat<n_rows,n_cols,mat_entry_t> sum(mat);
  return sum += scalar;
}
/*!*************************************************************************************************
 * \brief   Subtract SmallMat from scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator-
(const mat_entry_t& scalar, const SmallMat<n_rows,n_cols,mat_entry_t>& mat)
{
  SmallMat<n_rows,n_cols,mat_entry_t> difference(mat);
  for (unsigned int index = 0; index < mat.size(); ++index)
    difference[index] = scalar - mat[index];
  return difference;
}
/*!*************************************************************************************************
 * \brief   Subtract scalar from SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator-
(const SmallMat<n_rows,n_cols,mat_entry_t>& mat, const mat_entry_t& scalar)
{
  SmallMat<n_rows,n_cols,mat_entry_t> difference(mat);
  return difference -= scalar;
}
/*!*************************************************************************************************
 * \brief   Multiply scalar with SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator*
(const mat_entry_t& scalar, const SmallMat<n_rows,n_cols,mat_entry_t>& mat)
{
  SmallMat<n_rows,n_cols,mat_entry_t> product(mat);
  return product *= scalar;
}
/*!*************************************************************************************************
 * \brief   Multiply SmallMat with scalar.
 *  
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator*
(const SmallMat<n_rows,n_cols,mat_entry_t>& mat, const mat_entry_t& scalar)
{
  SmallMat<n_rows,n_cols,mat_entry_t> product(mat);
  return product *= scalar;
}
/*!*************************************************************************************************
 * \brief   Divide scalar by \c SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator/
(const mat_entry_t& scalar, const SmallMat<n_rows,n_cols,mat_entry_t>& mat)
{
  SmallMat<n_rows,n_cols,mat_entry_t> quotient(mat);
  for (unsigned int index = 0; index < mat.size(); ++index)
  {
    if (mat[index] == 0.)  throw SmallMatDivByZeroException();
    quotient[index] = scalar / mat[index];
  }
  return quotient;
}
/*!*************************************************************************************************
 * \brief   Divide \c SmallMat by scalar.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
SmallMat<n_rows,n_cols,mat_entry_t> operator/
(const SmallMat<n_rows,n_cols,mat_entry_t>& mat, const mat_entry_t& scalar)
{
  SmallMat<n_rows,n_cols,mat_entry_t> quotient(mat);
  return quotient /= scalar;
}

// Norms of SmallMats:

/*!*************************************************************************************************
 * \brief   Absolute sum norm of SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
mat_entry_t norm_1(const SmallMat<n_rows,n_cols,mat_entry_t>& mat)
{
  static_assert( n_cols == 1 , "This is only implemented for vectors." );
  mat_entry_t norm = 0.;
  for (unsigned int index = 0; index < mat.size(); ++index)  norm += std::abs( mat[index] );
  return norm;
}
/*!*************************************************************************************************
 * \brief   Computes Euclidean norm of a SmallMat.
 * 
 * \todo    Fill information, when details clarified with Guido.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
mat_entry_t norm_2(const SmallMat<n_rows,n_cols,mat_entry_t>& mat)
{
  static_assert( n_cols == 1 , "This is only implemented for vectors." );
  return std::sqrt( scalar_product(mat, mat) );
}
/*!*************************************************************************************************
 * \brief   Maximum norm of SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
mat_entry_t norm_infty(const SmallMat<n_rows,n_cols,mat_entry_t>& mat)
{
  static_assert( n_cols == 1 , "This is only implemented for vectors." );
  mat_entry_t norm = std::abs( mat[0] );
  for (unsigned int index = 1; index < mat.size(); ++index)
    norm = std::max( norm, std::abs(mat[index]) );
  return norm;
}
/*!*************************************************************************************************
 * \brief   p-norm of SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
mat_entry_t norm_p(const SmallMat<n_rows,n_cols,mat_entry_t>& mat, const mat_entry_t power)
{
  static_assert( n_cols == 1 , "This is only implemented for vectors." );
  mat_entry_t norm = 0.;
  for (unsigned int index = 0; index < mat.size(); ++index)
    norm += std::pow( std::abs(mat[index]) , power );
  return std::pow( norm , 1. / power );
}

// Output of SmallMat:

/*!*************************************************************************************************
 * \brief   Fill \c stream with \c SmallMat.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template< unsigned int n_rows, unsigned int n_cols, typename mat_entry_t >
std::ostream& operator<< (std::ostream& stream, const SmallMat<n_rows,n_cols,mat_entry_t>& mat)
{
  for (unsigned int row = 0; row < n_rows; ++row)
  {
    for (unsigned int col = 0; col < n_cols; ++col)
        stream << " " << mat(row,col) << " ";
    stream << std::endl;
  }
  return stream;
}


// Derived classes:

/*!*************************************************************************************************
 * \brief   A SmallSquareMat is a SmallMat, which is square.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rows, typename mat_entry_t = double >
using SmallSquareMat = SmallMat<n_rows, n_rows, mat_entry_t>;

/*!*************************************************************************************************
 * \brief   A SmallVec is a SmallMat, where the standard mat_entry_t is float.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rows, typename mat_entry_t = double >
using SmallVec = SmallMat<n_rows, 1, mat_entry_t>;

/*!*************************************************************************************************
 * \brief   A Point is a SmallMat, where the standard mat_entry_t is float.
 * 
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int n_rows, typename mat_entry_t = float >
using Point = SmallVec<n_rows, mat_entry_t>;