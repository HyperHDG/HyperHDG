#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/DenseLA.hxx>
#include <array>

/*!*************************************************************************************************
 * \brief   Mapping of a unit square to a parallelopipedon.
 *
 * \todo    Some functions of this class only work for one-dimensional hyperedges. These need to be
 *          generalized to arbitrary domensions of hyperedges.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template 
< unsigned int hyEdge_dimT, unsigned int space_dimT, typename map_float_t >
class Parallelopipedon
{
  public:
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of a hyperedge.
     *
     * \retval  hyEdge_dimT   The dimension of a hyperedge.
     **********************************************************************************************/
    static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
    /*!*********************************************************************************************
     * \brief   Returns the template parameter representing the dimension of the space.
     *
     * \retval  space_dimT       The dimension of the space.
     **********************************************************************************************/
    static constexpr unsigned int space_dim() { return space_dimT; }

  private:
    const Point<space_dimT,map_float_t> translation_;
    
    const SmallMat<space_dimT,hyEdge_dimT,map_float_t> matrix_;

    SmallMat<space_dimT,hyEdge_dimT, map_float_t> inner_normals_;
    SmallMat<space_dimT,space_dimT-hyEdge_dimT, map_float_t> outer_normals_;

  public:

    Parallelopipedon
    ( 
      const Point<space_dimT,map_float_t>& translation,
      const SmallMat<space_dimT,hyEdge_dimT,map_float_t>& matrix
    )
    : translation_(translation), matrix_(matrix) { }

    map_float_t haussdorff_hyEdge() const
    { return std::sqrt(std::abs(determinant( transposed_mat_times_mat(matrix_, matrix_) ))); }

    map_float_t haussdorfh_hyNode(const unsigned int index) const
    {
      SmallMat<space_dimT,hyEdge_dimT-1,map_float_t> mat_face;
      for (unsigned int i = 0; i < hyEdge_dimT; ++i)  if (i != index)
        mat_face.set_column(i - (i > index), matrix_.get_column(i));
      return std::sqrt(std::abs(determinant( transposed_mat_times_mat(mat_face, mat_face) )));
    }

    Point<space_dimT,map_float_t> inner_normal(const unsigned int index)
    {
      hy_assert( index < hyEdge_dimT ,
                 "The index of the inner normal must not be bigger than their amount." );

      Point<space_dimT,map_float_t> normal = inner_normals_.get_column(index);
      for (unsigned int i = 0; i < normal.size(); ++i)  if (normal[i] != 0.)  return normal;

      SmallMat<space_dimT,space_dimT-1,map_float_t> other_vectors;
      for (unsigned int i = 0; i < space_dimT-hyEdge_dimT; ++i)
        other_vectors.set_column(i, outer_normal(i));
      for (unsigned int i = 0; i < hyEdge_dimT; ++i)  if (i != index)
        other_vectors.set_column(i + space_dimT-hyEdge_dimT - (i > index), matrix_.get_column(i));

      normal = qr_decomp_q(other_vectors).get_column(space_dimT-1);
      map_float_t scalar_pdct = scalar_product(normal, matrix_.get_column(index));
      hy_assert( scalar_pdct != 0. , "Scalar product must not be zero." );
      if (scalar_pdct < 0)  normal *= -1.;

      inner_normals_.set_column(index, normal);
      return normal;
    }
    
    Point<space_dimT,map_float_t> outer_normal(const unsigned int index)
    {
      hy_assert( index < space_dimT-hyEdge_dimT ,
                 "The index of the outer normal must not be bigger than their amount." );

      Point<space_dimT,map_float_t> normal = outer_normals_.get_column(index);
      for (unsigned int i = 0; i < normal.size(); ++i)  if (normal[i] != 0.)  return normal;

      SmallMat<space_dimT,space_dimT,map_float_t> helper = qr_decomp_q(matrix_);
      for (unsigned int i = 0; i < space_dimT - hyEdge_dimT; ++i)
        outer_normals_.set_column(i, helper.get_column(i + hyEdge_dimT));

      normal = outer_normals_.get_column(index);
      return normal;
    }
}; // end class File
