#pragma once // Ensure that file is included only once in a single compilation.

/*!*************************************************************************************************
 * \brief   Helper class containing numbers and functions related to hypercubes.
 * 
 * \tparam  dim     Unsigned integer indicating the dimension of the considered hypercube.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int dim>
struct Hypercube
{
  /*!***********************************************************************************************
   * \brief   Number of faces of the hypercube.
   ************************************************************************************************/ 
  static constexpr unsigned int n_faces () { return 2 * dim; }
  /*!***********************************************************************************************
   * \brief   Number of vertices of the hypercube.
   ************************************************************************************************/ 
  static constexpr unsigned int n_vertices () { return 1 << dim; }
  /*!***********************************************************************************************
   * \brief   Return \c n to the power \c dim, which is the size of the dim-dimensional tensor of
   *          dimension n.
   ************************************************************************************************/ 
  static constexpr unsigned int pow (unsigned int n)
  {
    unsigned int result = 1;
    for (unsigned int i = 0; i < dim; ++i)  result *= n;
    return result;
  }
}; // end of struct Hypercube
