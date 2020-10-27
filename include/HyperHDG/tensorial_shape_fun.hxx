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
    decomposition[0] = index;
  else
  {
    for (unsigned int dim = 0; dim < dimT; ++dim)
    {
      decomposition[dim] = index % range;
      index /= range;
    }
  }
  return decomposition;
}

/*!*************************************************************************************************
 * \brief   Class that handles the evaluation of tensorial shape functions in tensor product points.
 *
 * \tparam  dimT              Local dimension of the shape function's domain or the point.
 * \tparam  abscissas_sizeT   size of abszissas array.
 * \tparam  abscissa_float_t  Floating type for the abscissa values.
 * \tparam  shape_t           type of 1D shapefunctions.
 * \tparam  polynomial_degree polynomial degree of 1D shape funcion space.
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
  static constexpr unsigned int shapefunction_count = std::pow(polynomial_degree + 1, dimT);
  static constexpr unsigned int point_count = std::pow(abscissas_sizeT, dimT);
  std::array<unsigned int, polynomial_degree + 1> indices_1D;
  std::array<std::array<return_t, abscissas_sizeT>, polynomial_degree + 1> values_1D;
  inline void generate_values_for_abscissas(
    const std::array<abscissa_float_t, abscissas_sizeT> abscissas)
  {
    values_1D = shape_fct_eval<return_t, shape_t>(indices_1D, abscissas);
  }
  inline void generate_indices()
  {
    for (unsigned int i = 0; i < polynomial_degree + 1; ++i)
      indices_1D[i] = i;
  }

 public:
  TensorialShapeFunctionEvaluation(std::array<abscissa_float_t, abscissas_sizeT> abscissas)
  {
    generate_indices();
    generate_values_for_abscissas(abscissas);
  }
  unsigned int get_point_count() { return point_count; }
  unsigned int get_shapefunction_count() { return shapefunction_count; }
  inline return_t evaluate_in_tensorial_point(unsigned int tensorial_function_index,
                                              unsigned int tensorial_point_index)
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
   ************************************************************************************************/
  template <unsigned int coefficient_dimension>
  inline std::array<std::array<return_t, point_count>, coefficient_dimension>
  evaluate_linear_combination_in_all_tensorial_points(
    std::array<return_t, coefficient_dimension * shapefunction_count>
      linear_combination_coefficients)
  {
    std::array<std::array<return_t, point_count>, coefficient_dimension> return_values;

    for (unsigned int point = 0; point < point_count; ++point)
      for (unsigned int d = 0; d < coefficient_dimension; ++d)
        return_values[d][point] = 0;

    double function_value;
    for (unsigned int function = 0; function < shapefunction_count; ++function)
    {
      for (unsigned int point = 0; point < point_count; ++point)
      {
        function_value = evaluate_in_tensorial_point(function, point);
        for (unsigned int d = 0; d < coefficient_dimension; ++d)
        {
          return_values[d][point] +=
            function_value * linear_combination_coefficients[d * shapefunction_count + function];
        }
      }
    }
    return return_values;
  }
};
