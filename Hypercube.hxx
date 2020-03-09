#pragma once // Ensure that file is included only once in a single compilation.

/*!*************************************************************************************************
 * \brief   Helper class containing numbers and functions related to hypercubes.
 * 
 * \todo    All functions related to dim-dimensional tensor can go here.
 * 
 * \tparam  dim     Unsigned integer indicating the dimension of the considered hypercube.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
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
   * \brief   Return `n` to the power `dim`, which is the size of the `dim`-dimensional tensor of
   *          dimension `n`.
   ************************************************************************************************/ 
  static constexpr unsigned int pow(unsigned int n)
  {
    int result = 1;
    for (unsigned int i=0;i<dim;++i)  result *= n;
    return result;
  }
}; // end of struct Hypercube
