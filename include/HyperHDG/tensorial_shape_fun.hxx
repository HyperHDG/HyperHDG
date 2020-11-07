#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/shape_fun_1d.hxx>

#include <array>

/*!*************************************************************************************************
 * \brief   Decompose index of tensorial shape function or tensorial point with respect to the
 *          local dimension.
 *
 * The index of the shape function needs to be decomposed into dim indices of one-dimensional shape
 * functions or points.
 *
 * \tparam  dimT          Local dimension of the shape function's domain or the point.
 * \tparam  range         Range (maximum) of the 1D indices.
 * \param   index         Local index of the shape function or point.
 * \retval  decomposition Array consisting of respective one-dimensional indices.
 **************************************************************************************************/
template <unsigned int dimT, unsigned int range>
inline std::array<unsigned int, std::max(dimT, 1U)> index_decompose(unsigned int index)
{
  std::array<unsigned int, std::max(dimT, 1U)> decomposition;
  if (decomposition.size() == 1)
  {
    hy_assert(index < range, "The index must not exceed its possible range.");
    decomposition[0] = index;
  }
  else
  {
    for (unsigned int dim = 0; dim < dimT; ++dim)
    {
      decomposition[dim] = index % range;
      index /= range;
    }
    hy_assert(index == 0, "Index initially exceeded given maximum value range^dimT.");
  }
  return decomposition;
}

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
template <unsigned int dimT,
          typename return_t,
          typename shape_t,
          unsigned int polynomial_degree,
          unsigned int abscissas_sizeT,
          typename abscissa_float_t>
class TensorialShapeFunctionEvaluation
{
 private:
  /*!***********************************************************************************************
   * \brief   Number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_shapefun_ = std::pow(polynomial_degree + 1, dimT);
  /*!***********************************************************************************************
   * \brief   Number of dimT-dimensional abscissa values.
   ************************************************************************************************/
  static constexpr unsigned int n_points_ = std::pow(abscissas_sizeT, dimT);
  /*!***********************************************************************************************
   * \brief   Values of one-dimensional shape-functions at abscissas.
   ************************************************************************************************/
  const std::array<std::array<return_t, abscissas_sizeT>, polynomial_degree + 1> values_1D;
  /*!***********************************************************************************************
   * \brief   Function that fills values1D.
   ************************************************************************************************/
  inline std::array<std::array<return_t, abscissas_sizeT>, polynomial_degree + 1>
  generate_values_for_abscissas(const std::array<abscissa_float_t, abscissas_sizeT> abscissas)
  {
    std::array<unsigned int, polynomial_degree + 1> indices_1D;
    for (unsigned int i = 0; i < polynomial_degree + 1; ++i)
      indices_1D[i] = i;
    return shape_fct_eval<return_t, shape_t>(indices_1D, abscissas);
  }

 public:
  /*!***********************************************************************************************
   * \brief   Return umber of dimT-dimensional abscissa values.
   ************************************************************************************************/
  static constexpr unsigned int n_points() { return n_points_; }
  /*!***********************************************************************************************
   * \brief   Return number of shape functions that span the local polynomial space.
   ************************************************************************************************/
  static constexpr unsigned int n_shapefun() { return n_shapefun_; }
  /*!***********************************************************************************************
   * \brief   Constructor.
   *
   * \param   abscissas                 Array containing one-dimensional abscissa values.
   ************************************************************************************************/
  TensorialShapeFunctionEvaluation(std::array<abscissa_float_t, abscissas_sizeT> abscissas)
  : values_1D(generate_values_for_abscissas(abscissas)){};
  /*!***********************************************************************************************
   * \brief   Return shape function of given index at abscissa of given index.
   *
   * \param   tensorial_function_index  Index of the tensorial function.
   * \param   tensorial_point_index     Index of tensorial point to be evaluated.
   * \retval  function_value            Value of function evaluated at point.
   ************************************************************************************************/
  inline return_t evaluate_in_tensorial_point(unsigned int tensorial_function_index,
                                              unsigned int tensorial_point_index) const
  {
    std::array<unsigned int, std::max(dimT, 1U)> decomposed_function_index =
      index_decompose<dimT, polynomial_degree + 1>(tensorial_function_index);
    std::array<unsigned int, std::max(dimT, 1U)> decomposed_point_index =
      index_decompose<dimT, abscissas_sizeT>(tensorial_point_index);
    return_t function_value = 1.;
    for (unsigned int d = 0; d < std::max(dimT, 1U); ++d)
    {
      function_value *= values_1D[decomposed_function_index[d]][decomposed_point_index[d]];
    }
    return function_value;
  }
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
