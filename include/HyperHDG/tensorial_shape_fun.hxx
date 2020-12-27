#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/shape_fun_1d.hxx>

#include <array>

#ifndef TENSORIAL_INDEX_DECOMPOSE
#define TENSORIAL_INDEX_DECOMPOSE
/*!*************************************************************************************************
 * \brief   Decompose index of tensorial shape function or tensorial point with respect to the
 *          local dimension.
 *
 * The index of the shape function needs to be decomposed into dim indices of one-dimensional shape
 * functions or points.
 *
 * \tparam  dimT          Local dimension of the shape function's domain or the point.
 * \param   index         Local index of the shape function or point.
 * \param   range         Range (maximum) of the 1D indices.
 * \retval  decomposition Array consisting of respective one-dimensional indices.
 **************************************************************************************************/
template <unsigned int dimT>
inline std::array<unsigned int, std::max(dimT, 1U)> tensorial_index_decompose(unsigned int index,
                                                                              unsigned int range)
{
  std::array<unsigned int, std::max(dimT, 1U)> decomposition;
  for (unsigned int dim = 0; dim < decomposition.size(); ++dim)
  {
    decomposition[dim] = index % range;
    index /= range;
  }
  hy_assert(index == 0, "Index initially exceeded given maximum value range^dimT.");
  return decomposition;
}
#endif  // ifndef TENSORIAL_INDEX_DECOMPOSE

/*!*************************************************************************************************
 * \brief   Class that handles the evaluation of tensorial shape functions in tensor product points.
 *
 * \tparam  dimT              Local dimension of the shape function's domain or the point.
 * \tparam  return_t          The floating type that is returned.
 * \tparam  polynomial_degree Maximum polynomial degree that can be evaluated.
 * \tparam  abscissas_sizeT   Size of abszissas array.
 * \tparam  abscissa_float_t  Floating type for the abscissa values.
 * \tparam  shape_t           Type of 1D shapefunctions.
 *
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 * \authors   Simon Schmidt, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename shape_t,
          unsigned int dimT
          unsigned int polynomial_degree>
struct ShapeFunTensorial
{
  typedef shape_t shape_fun_1d_t;
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_shapefun_ = Hypercube<dimT>::pow(polynomial_degree + 1);
  /*!***********************************************************************************************
   * \brief   Evaluate a multidimensional linear combination of tensor product shapefunctions in
   *          tensor product points.
   *
   * \tparam  coefficient_dimension     Amount of linear combination to be evaluated.
   * \param   lin_comb_coeffs           Double array comprising the linear coefficients.
   * \retval  return_values             Evaluation of linear combinations in tensorial points.
   ************************************************************************************************/
  template <unsigned int coefficient_dimension>
  inline std::array<std::array<return_t, n_points_>, coefficient_dimension>
  evaluate_linear_combination_in_all_tensorial_points(
    std::array<return_t, coefficient_dimension * n_shapefun_> lin_comb_coeffs) const
  {
    constexpr unsigned int n_points_ = Hypercube<dimT>::pow(abscissas_sizeT);
    std::array<std::array<return_t, n_points_>, coefficient_dimension> return_values;

    for (unsigned int point = 0; point < n_points_; ++point)
      for (unsigned int d = 0; d < coefficient_dimension; ++d)
        return_values[d][point] = 0;

    double function_value;
    for (unsigned int function = 0; function < n_shapefun_; ++function)
    {
      for (unsigned int point = 0; point < n_points_; ++point)
      {
        function_value = evaluate_in_tensorial_point(function, point);
        for (unsigned int d = 0; d < coefficient_dimension; ++d)
        {
          return_values[d][point] += function_value * lin_comb_coeffs[d * n_shapefun_ + function];
        }
      }
    }
    return return_values;
  }
};
