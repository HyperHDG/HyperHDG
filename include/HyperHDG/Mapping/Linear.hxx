#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/DenseLA.hxx>

#include <array>
#include <memory>

namespace Mapping
{

/*!*************************************************************************************************
 * \brief   Mapping of a unit square to a Linear.
 *
 * \todo    Update doxygen in this whole file!!!
 * \todo    Maybe add more template parameters for matrix and vector implementations?!?!
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template 
< unsigned int hyEdge_dimT, unsigned int space_dimT, typename map_float_t >
class Linear
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

    std::shared_ptr< SmallSquareMat<space_dimT,map_float_t> > matrix_q_;
    std::shared_ptr< SmallSquareMat<hyEdge_dimT,map_float_t> > matrix_r_;

    void make_qr_if_needed()
    {
      if (matrix_q_ && matrix_r_)  return;
      hy_assert ( !matrix_q_ && !matrix_r_ ,
                  "Construction of only one matrix has not been implemented, yet!" );
      matrix_q_ = std::make_shared< SmallSquareMat<space_dimT,map_float_t> > ();
      matrix_r_ = std::make_shared< SmallSquareMat<hyEdge_dimT,map_float_t> > ();
      qr_decomp( matrix_ , *matrix_q_, *matrix_r_ );
    }

  public:

    Linear
    ( 
      const Point<space_dimT,map_float_t>& translation,
      const SmallMat<space_dimT,hyEdge_dimT,map_float_t>& matrix
    )
    : translation_(translation), matrix_(matrix) { }

    map_float_t functional_determinant_hyEdge()
    { 
      make_qr_if_needed();
      map_float_t determinant = 1.;
      for (unsigned int i = 0; i < hyEdge_dimT; ++i)  determinant[i] *= matrix_r_->operator()(i,i);
      return determinant;
    }

    map_float_t functional_determinant_hyNode(const unsigned int index) const
    {
      if constexpr ( hyEdge_dimT == 1 )  return 1.;
      SmallMat<space_dimT,hyEdge_dimT-1,map_float_t> mat_face;
      for (unsigned int i = 0; i < hyEdge_dimT; ++i)  if (i != index)
        mat_face.set_column(i - (i > index), matrix_.get_column(i));
      return determinant( mat_face );
    }

    SmallVec<space_dimT,map_float_t> matrix_column(const unsigned int col) const
    { return matrix_.get_column(col); }

    template<unsigned int n_vec>
    SmallMat<space_dimT,n_vec,map_float_t> map_reference_to_physical
    ( const SmallMat<hyEdge_dimT,n_vec,map_float_t>& mat ) const
    { return matrix_ * mat + translation_; }

    template<unsigned int n_vec> // Todo: Implement matrix solution via qr deconposition for this
    SmallMat<hyEdge_dimT,n_vec,map_float_t> map_reference_to_physical
    ( const SmallMat<space_dimT,n_vec,map_float_t>& mat ) const
    { return (mat - translation_) / matrix_; }

    const SmallSquareMat<hyEdge_dimT,map_float_t>& mat_r()
    {
      make_qr_if_needed();
      return *matrix_r_;
    }

    SmallVec<hyEdge_dimT,map_float_t> hyEdge_dim_normal(const unsigned int index)
    {
      hy_assert( index < hyEdge_dimT ,
                 "The index of the searched normal must not be bigger than their amount." );
      if constexpr (hyEdge_dimT == 1)  return SmallVec<hyEdge_dimT,map_float_t>(1.);
      
      make_qr_if_needed();
      SmallMat<hyEdge_dimT,hyEdge_dimT-1,map_float_t> other_vectors;
      for (unsigned int i = 0; i < hyEdge_dimT; ++i)  if (i != index)
        other_vectors.set_column(i - (i > index), matrix_r_->get_column(i));
      
      Point<hyEdge_dimT,map_float_t> normal = qr_decomp_q(other_vectors).get_column(hyEdge_dimT-1);
      map_float_t scalar_pdct = scalar_pdct(normal, matrix_r_->get_column(index));
      hy_assert( scalar_pdct != 0., "Scalar product must not be zero!" );
      if (scalar_pdct < 0.)  normal *= -1.;
      return normal;
    }

    Point<space_dimT,map_float_t> inner_normal(const unsigned int index)
    {
      hy_assert( index < hyEdge_dimT ,
                 "The index of the inner normal must not be bigger than their amount." );
      if constexpr (space_dimT == 1)  return SmallVec<space_dimT,map_float_t>(1.);

      SmallMat<space_dimT,space_dimT-1,map_float_t> other_vectors;
      for (unsigned int i = 0; i < space_dimT-hyEdge_dimT; ++i)
        other_vectors.set_column(i, outer_normal(i));
      for (unsigned int i = 0; i < hyEdge_dimT; ++i)  if (i != index)
        other_vectors.set_column(i + space_dimT-hyEdge_dimT - (i > index), matrix_.get_column(i));

      Point<space_dimT, map_float_t> normal = qr_decomp_q(other_vectors).get_column(space_dimT-1);
      map_float_t scalar_pdct = scalar_product(normal, matrix_.get_column(index));
      hy_assert( scalar_pdct != 0. , "Scalar product must not be zero." );
      if (scalar_pdct < 0)  normal *= -1.;
      return normal;
    }
    
    Point<space_dimT,map_float_t> outer_normal(const unsigned int index)
    {
      hy_assert( index < space_dimT-hyEdge_dimT ,
                 "The index of the outer normal must not be bigger than their amount." );

      make_qr_if_needed();
      return matrix_q_->get_column(index + hyEdge_dimT);;
    }
}; // end class File

} // end of namespace Mapping