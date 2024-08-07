#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hy_assert.hxx>

namespace Mapping
{
/*!*************************************************************************************************
 * \brief   Mapping of a unit hypercube to a parallelotope --- can also be used for simplices.
 *
 * The affine-linear mapping of a \c hyEdge_dimT dimensional unit square to an parallelotope which
 * lives in \c space_dimT dimensions.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT, unsigned int space_dimT, typename map_float_t>
class Linear
{
 public:
  /*!***********************************************************************************************
   * \brief   Returns the template parameter representing the dimension of a hyperedge.
   *
   * \retval  hyEdge_dimT   The dimension of a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Returns the template parameter representing the dimension of the space.
   *
   * \retval  space_dimT       The dimension of the space.
   ************************************************************************************************/
  static constexpr unsigned int space_dim() { return space_dimT; }

 private:
  /*!***********************************************************************************************
   * \brief   The translation of the affine mapping from reference to physical element.
   ************************************************************************************************/
  SmallVec<space_dimT, map_float_t> translation_;
  /*!***********************************************************************************************
   * \brief   The matrix of the affine mapping from reference to physical element.
   ************************************************************************************************/
  SmallMat<space_dimT, hyEdge_dimT, map_float_t> matrix_;
  /*!***********************************************************************************************
   * \brief   Matrix Q of QR decomposition of the matrix of the affine transformation.
   *
   * This matrix is normalized to have det(Q) = +1.
   ************************************************************************************************/
  SmallSquareMat<space_dimT, map_float_t> matrix_q_;
  /*!***********************************************************************************************
   * \brief   Matrix R of QR decomposition of the matrix of the affine transformation.
   *
   * This matrix actually is a rectangular matrix. Due to our assumption that the space dimension
   * is larger than or equal to the dimension of a hyperedge, the matrix \c matrix and its R of
   * the QR decomposition have at least as many rows, as they have columns. Nonetheless, R is an
   * upper triangular matrix. This allows to remove the zero rows (actually below our saved R)
   * without loss of generality, but saving some space.
   *
   * This matrix is normalized to have non-negative diagonal entries, except for the very first.
   ************************************************************************************************/
  SmallSquareMat<hyEdge_dimT, map_float_t> matrix_r_;

 public:
  /*!***********************************************************************************************
   * \brief   Construct affine-linear mapping.
   ************************************************************************************************/
  Linear() = default;
  /*!***********************************************************************************************
   * \brief   Construct affine-linear mapping.
   *
   * \param   translation   The translation of the afine-linear mapping.
   * \param   matrix        The matrix of the affine-linear mapping.
   ************************************************************************************************/
  Linear(const Point<space_dimT, map_float_t>& translation,
         const SmallMat<space_dimT, hyEdge_dimT, map_float_t>& matrix)
  : translation_(translation), matrix_(matrix)
  {
    qr_decomp(matrix, matrix_q_, matrix_r_);
  }
  /*!***********************************************************************************************
   * \brief   The functional determinant of the affine-linear mappring.
   *
   * The determinant, in general is only defined for square matrices. In our case, the determinant
   * is related to the ratio of Lebesque measure of the reference cell and the Haussdorff measure
   * of the physical cell. This induces a natural generalization to rectangular matrices (with at
   * least as many rows as columns), where the generalized determinant can be defined as the
   * product of the diagonal entries of matrix R of the QR decomposition if normalized as it is
   * described above.
   *
   * This determinant corresponds (up to its sign) to the general Haussdorff transformation
   * formula of measures by a factor of g = sqrt( (D Phi)^T (D Phi) ).
   * The difference in the sign will become important for the transformation of gradients.
   ************************************************************************************************/
  map_float_t functional_determinant_hyEdge() const
  {
    map_float_t determinant = 1.;
    for (unsigned int i = 0; i < hyEdge_dimT; ++i)
      determinant *= matrix_r_(i, i);
    return determinant;
  }
  /*!***********************************************************************************************
   * \brief   The functional determinant of the affine-linear mapping of a hypernode.
   *
   * For details consider the description of \c functional_determinant_hyEdge(), and that a
   * hypernode is spanned by all, but one columns of the matrix for the affine-linear mapping.
   *
   * \param   index         Index of the vector in the matrix which is \b not related to the node.
   ************************************************************************************************/
  map_float_t functional_determinant_hyNode(const unsigned int index) const
  {
    if constexpr (hyEdge_dimT == 1)
      return 1.;
    else
    {
      SmallMat<space_dimT, hyEdge_dimT - 1, map_float_t> mat_face;
      for (unsigned int i = 0; i < hyEdge_dimT; ++i)
        if (i != index)
          mat_face.set_column(i - (i > index), matrix_.get_column(i));
      return determinant(mat_face);
    }
  }
  /*!***********************************************************************************************
   * \brief   Return vector representing matrix column of specified index.
   *
   * \param   col         The index of the column that is to be returned.
   ************************************************************************************************/
  SmallVec<space_dimT, map_float_t> matrix_column(const unsigned int col) const
  {
    return matrix_.get_column(col);
  }

  /*!***********************************************************************************************
   * \brief   Map n_vec points from reference to physical element.
   *
   * \param   points        Matrix whose columns consist of the points to be mapped.
   * \retval  phy_points    Matrix whose columns consist of the mapped points.
   ************************************************************************************************/
  template <unsigned int n_vec>
  SmallMat<space_dimT, n_vec, map_float_t> map_reference_to_physical(
    const SmallMat<hyEdge_dimT, n_vec, map_float_t>& points) const
  {
    return matrix_ * points + translation_;
  }
  /*!***********************************************************************************************
   * \brief   Map n_vec points from physical to element.
   *
   * \param   phy_points    Matrix whose columns consist of the mapped points.
   * \retval  points        Matrix whose columns consist of the points to be mapped.
   ************************************************************************************************/
  template <unsigned int n_vec>
  SmallMat<hyEdge_dimT, n_vec, map_float_t> map_physical_to_reference(
    const SmallMat<space_dimT, n_vec, map_float_t>& phy_points) const
  {
    return (phy_points - translation_) / matrix_;
  }
  /*!***********************************************************************************************
   * \brief   Return matrix Q of the QR decomposition.
   ************************************************************************************************/
  const SmallSquareMat<space_dimT, map_float_t>& mat_q() const { return matrix_q_; }
  /*!***********************************************************************************************
   * \brief   Return matrix R of the QR decomposition.
   ************************************************************************************************/
  const SmallSquareMat<hyEdge_dimT, map_float_t>& mat_r() const { return matrix_r_; }
  /*!***********************************************************************************************
   * \brief   Return local normal of given index.
   *
   * Return outer unit normal with respect to the hypernode which is spanned by all columns of R,
   * but the column of the given index. This is an element of the same dimension as the reference
   * square.
   *
   * \param   index         Index of the vector in the matrix which is \b not related to the node.
   ************************************************************************************************/
  SmallVec<hyEdge_dimT, map_float_t> local_normal(const unsigned int index) const
  {
    hy_assert(index < hyEdge_dimT,
              "The index of the searched normal must not be bigger than their amount.");
    if constexpr (hyEdge_dimT == 1)
    {
      if (matrix_r_(0, 0) > 0)
        return SmallVec<hyEdge_dimT, map_float_t>(-1.);
      else
        return SmallVec<hyEdge_dimT, map_float_t>(+1.);
    }
    else
    {
      SmallMat<hyEdge_dimT, hyEdge_dimT - 1, map_float_t> other_vectors;
      for (unsigned int i = 0; i < hyEdge_dimT; ++i)
        if (i != index)
          other_vectors.set_column(i - (i > index), matrix_r_.get_column(i));

      SmallVec<hyEdge_dimT, map_float_t> normal =
        qr_decomp_q(other_vectors).get_column(hyEdge_dimT - 1);
      map_float_t scalar_pdct = scalar_product(normal, matrix_r_.get_column(index));
      hy_assert(scalar_pdct != 0., "Scalar product must not be zero!");
      if (scalar_pdct > 0.)
        normal *= -1.;
      return normal;
    }
  }
  /*!***********************************************************************************************
   * \brief   Return inner normal of given index.
   *
   * Return outer unit normal with respect to the hypernode which is spanned by all columns of the
   * transformation matrix but the column of the given index. The vector has to be in the span of
   * the columns of the transformation matrix. This is an element of the same dimension as the
   * full space.
   *
   * \param   index         Index of the vector in the matrix which is \b not related to the node.
   ************************************************************************************************/
  SmallVec<space_dimT, map_float_t> inner_normal(const unsigned int index) const
  {
    // hy_assert(index < hyEdge_dimT,
    //           "The index of the inner normal must not be bigger than their amount.");
    // if constexpr (space_dimT == 1)
    //   return SmallVec<space_dimT, map_float_t>(1.);

    // SmallMat<space_dimT, space_dimT - 1, map_float_t> other_vectors;
    // for (unsigned int i = 0; i < space_dimT - hyEdge_dimT; ++i)
    //   other_vectors.set_column(i, outer_normal(i));
    // for (unsigned int i = 0; i < hyEdge_dimT; ++i)
    //   if (i != index)
    //     other_vectors.set_column(i + space_dimT - hyEdge_dimT - (i > index),
    //     matrix_.get_column(i));

    // SmallVec<space_dimT, map_float_t> normal =
    //   qr_decomp_q(other_vectors).get_column(space_dimT - 1);
    // map_float_t scalar_pdct = scalar_product(normal, matrix_.get_column(index));
    // hy_assert(scalar_pdct != 0., "Scalar product must not be zero.");
    // if (scalar_pdct > 0)
    //   normal *= -1.;
    // return normal;
    hy_assert(index < hyEdge_dimT,
              "The index of the outer normal must not be bigger than their amount.");

    return matrix_q_.get_column(index);
  }
  /*!***********************************************************************************************
   * \brief   Return outer normal of given index.
   *
   * Return unit normal with respect to the hyperedge within the full space.
   *
   * \param   index         Index of the normal.
   ************************************************************************************************/
  SmallVec<space_dimT, map_float_t> outer_normal(const unsigned int index) const
  {
    hy_assert(index < space_dimT - hyEdge_dimT,
              "The index of the outer normal must not be bigger than their amount.");

    return matrix_q_.get_column(index + hyEdge_dimT);
  }
  /*!***********************************************************************************************
   * \brief   Return matrix associated to affine-linear transformation.
   ************************************************************************************************/
  const SmallMat<space_dimT, hyEdge_dimT, map_float_t>& matrix() const { return matrix_; }
  /*!***********************************************************************************************
   * \brief   Return shifting vector associated to affine-linear transformation.
   ************************************************************************************************/
  const SmallVec<space_dimT, map_float_t>& translation() const { return translation_; }
};  // end class File

}  // end of namespace Mapping
