#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

#include <HyperHDG/shape_fun_1d.hxx>
#include <array>

/*!*********************************************************************************************
 * \brief   Decompose index of tensorial shape function or tensorial point with respect to the
 * local dimension.
 *
 * Tthe index of the shape function needs to be decomposed into dim indices of one-
 * dimensional shape functions or points.
 *
 * \tparam  dimT          Local dimension of the shape function's domain or the point.
 * \tparam  range         Range (maximum) of the 1D indices.
 * \param   index         Local index of the shape function or point.
 * \param   decomposition Array consisting of respective one-dimensional indices.
 **********************************************************************************************/
template<unsigned int dimT, unsigned int range>
inline std::array<unsigned int, std::max(dimT,1U)> index_decompose ( unsigned int index )
{
  std::array<unsigned int, std::max(dimT,1U)> decomposition;
  if ( decomposition.size() == 1 )  decomposition[0] = index;
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


/*!*********************************************************************************************
 * \brief   Evaluate tensor product shape function with (multidimensional index) in all tensor product points of abscissas
 *
 *
 * \tparam  dimT          Local dimension of the shape function's domain or the point.
 * \tparam  abscissas_sizeT         size of abszissas array.
 * \tparam  abscissa_float_t         Floating type for the abscissa values.
 * \tparam  shape_t         type of 1D shapefunctions
 * \param   index         index of tensor product shape
 * \param   abscissas Abscissas of the supporting points.
 **********************************************************************************************/
template < unsigned int dimT, unsigned int abscissas_sizeT, typename return_t, typename abscissa_float_t, typename shape_t>
static inline std::array<return_t,(unsigned int) std::pow(abscissas_sizeT, dimT)> evaluate_in_tensorial_points(const unsigned int index, const std::array<abscissa_float_t, abscissas_sizeT>  abscissas)
{
  const unsigned int tensorial_point_count = std::pow(abscissas_sizeT, dimT);
  std::array<return_t,tensorial_point_count> point_values;
  std::array<return_t, abscissas_sizeT>values_1D
      = shape_fct_eval<return_t, shape_t>(index, abscissas);
    for (unsigned int point = 0; point < tensorial_point_count; ++point) {
      std::array<unsigned int, dimT> decomposed_point_index = index_decompose<dimT, abscissas_sizeT>(point);
      return_t function_value = 1.;
      for (unsigned int d = 0; d < dimT; ++d) {
        function_value *= values_1D[decomposed_point_index[d]];
      }
      point_values[point] = function_value;
    }
  return point_values;
}


/*!*********************************************************************************************
 * \brief   Evaluate tensor product shape function with (multidimensional index) in all tensor product points of abscissas
 * thet are on the boundary. We require that abscissas[0] = 0  and abscissas[1] = 1
 *
 *
 * \tparam  dimT          Local dimension of the shape function's domain or the point.
 * \tparam  abscissas_sizeT         size of abszissas array.
 * \tparam  abscissa_float_t         Floating type for the abscissa values.
 * \tparam  shape_t         type of 1D shapefunctions
 * \param   index         index of tensor product shape
 * \param   abscissas Abscissas of the supporting points.
 **********************************************************************************************/
template < typename return_t,  typename abscissa_float_t, typename shape_t, unsigned int dimT, unsigned int abscissas_sizeT>
static inline std::array<return_t,(unsigned int) std::pow(abscissas_sizeT, dimT)*2*dimT> evaluate_in_tensorial_boundary_points(const unsigned int index, const std::array<abscissa_float_t, abscissas_sizeT>  abscissas)
{
  const unsigned int tensorial_point_count = std::pow(abscissas_sizeT, dimT-1); //only count points on single face
  const unsigned int total_point_count = tensorial_point_count*2*dimT;
  std::array<return_t,total_point_count> point_values;

  std::array<return_t, abscissas_sizeT>values_1D
      = shape_fct_eval<return_t, shape_t>(index, abscissas);

  for (unsigned int point = 0; point < total_point_count; ++point) {
    unsigned int face_of_point = point / tensorial_point_count;
    unsigned int local_point_index = point % tensorial_point_count;
    std::array<unsigned int, dimT> decomposed_local_point_index = index_decompose<dimT, abscissas_sizeT>(local_point_index);
    return_t function_value = 1.;
    unsigned int subd = 0;
    for (unsigned int d = 0; d < dimT; ++d) {
      if(d == face_of_point/2)
      {
        if(d%2 == 0)
          function_value *= values_1D[0];
        else
          function_value *= values_1D.back();
      }
      function_value *= values_1D[decomposed_local_point_index[subd++]];
    }
    point_values[point] = function_value;
  }
  return point_values;
}

/*!*********************************************************************************************
 * \brief   Evaluate all tensor product shape functions in all tensor product points of abscissas
 *
 *
 * \tparam  dimT          Local dimension of the shape function's domain or the point.
 * \tparam  abscissas_sizeT         size of abszissas array.
 * \tparam  abscissa_float_t         Floating type for the abscissa values.
 * \tparam  shape_t         type of 1D shapefunctions.
 * \tparam  polynomial_degree         polynomial degree of 1D shape funcion space.
 * \param   abscissas Abscissas of the supporting points.
 **********************************************************************************************/
template < typename return_t, typename abscissa_float_t, typename shape_t,  unsigned int dimT, unsigned int abscissas_sizeT, unsigned int polynomial_degree>
static inline std::array<std::array<return_t, (unsigned int)std::pow(abscissas_sizeT, dimT)>,
    (unsigned int) std::pow(polynomial_degree+1, dimT)>
    evaluate_all_in_tensorial_points(const std::array<abscissa_float_t, abscissas_sizeT>  abscissas)
{
  const unsigned int tensorial_shapefunction_count = std::pow(polynomial_degree+1,dimT);
  const unsigned int tensorial_point_count = std::pow(abscissas_sizeT, dimT);
  std::array<std::array<return_t,tensorial_point_count>, tensorial_shapefunction_count> point_values;

  std::array<unsigned int, polynomial_degree + 1> indices_1D;
  for(unsigned int i = 0; i < polynomial_degree + 1; ++i)
    indices_1D[i] = i;
  std::array<std::array<return_t, abscissas_sizeT>, polynomial_degree+1> values_1D
      = shape_fct_eval<return_t, shape_t>(indices_1D, abscissas);
  for(unsigned int function = 0; function < tensorial_shapefunction_count; ++function) {
    std::array<unsigned int, dimT> decomposed_function_index = index_decompose<dimT, polynomial_degree+1>(function);
    for (unsigned int point = 0; point < tensorial_point_count; ++point) {
      std::array<unsigned int, dimT> decomposed_point_index = index_decompose<dimT, abscissas_sizeT>(point);
      return_t function_value = 1.;
      for (unsigned int d = 0; d < dimT; ++d) {
        function_value *= values_1D[decomposed_function_index[d]][decomposed_point_index[d]];
      }
      point_values[function][point] = function_value;
    }
  }
  return point_values;
}

/*!*********************************************************************************************
 * \brief   Evaluate all tensor product shape functions in all tensor product points of abscissas
 * thet are on the boundary. We require that abscissas[0] = 0  and abscissas[1] = 1
 *
 *
 * \tparam  dimT          Local dimension of the shape function's domain or the point.
 * \tparam  abscissas_sizeT         size of abszissas array.
 * \tparam  abscissa_float_t         Floating type for the abscissa values.
 * \tparam  shape_t         type of 1D shapefunctions
 * \tparam  polynomial_degree         polynomial degree of 1D shape funcion space.
 * \param   abscissas Abscissas of the supporting points.
 **********************************************************************************************/
template < typename return_t,  typename abscissa_float_t, typename shape_t, unsigned int dimT, unsigned int abscissas_sizeT, unsigned int polynomial_degree>
static inline std::array<std::array<return_t,(unsigned int) std::pow(abscissas_sizeT, dimT-1)*2*dimT>,
                         (unsigned int) std::pow(polynomial_degree+1, dimT)>
       evaluate_all_in_tensorial_boundary_points(const std::array<abscissa_float_t, abscissas_sizeT>  abscissas)
{
  const unsigned int tensorial_shapefunction_count = std::pow(polynomial_degree+1,dimT);
  const unsigned int tensorial_point_count = std::pow(abscissas_sizeT, dimT-1); //only count points on single face
  const unsigned int total_point_count = tensorial_point_count*2*dimT;
  std::array<std::array<return_t, total_point_count>, tensorial_shapefunction_count> point_values;

  std::array<unsigned int, polynomial_degree + 1> indices_1D;
  for(unsigned int i = 0; i < polynomial_degree + 1; ++i)
    indices_1D[i] = i;

  std::array<std::array<return_t, abscissas_sizeT>, polynomial_degree+1> values_1D
      = shape_fct_eval<return_t, shape_t>(indices_1D, abscissas);

  for(unsigned int function = 0; function < tensorial_shapefunction_count; ++function) {
    std::array<unsigned int, dimT> decomposed_function_index = index_decompose<dimT, polynomial_degree + 1>(function);
    for (unsigned int point = 0; point < total_point_count; ++point) {
      unsigned int face_of_point = point / tensorial_point_count;
      unsigned int local_point_index = point % tensorial_point_count;
      std::array<unsigned int, dimT>
          decomposed_local_point_index = index_decompose<dimT, abscissas_sizeT>(local_point_index);
      return_t function_value = 1.;
      unsigned int subd = 0;
      for (unsigned int d = 0; d < dimT; ++d) {
        if (d == face_of_point / 2) {
          if (d % 2 == 0)
            function_value *= values_1D[function][0];
          else
            function_value *= values_1D[function].back();
        }
        function_value *= values_1D[function][decomposed_local_point_index[subd++]];
      }
      point_values[function][point] = function_value;
    }
  }
  return point_values;
}

/*!*********************************************************************************************
 * \brief   Evaluate for (possibly) multiple dimensions a linear combination of all tensor product shape functions
 * in all tensor product points of abscissas, each dimension can have a different linear combination
 *
 *
 *
 * \tparam  dimT          Local dimension of the shape function's domain or the point.
 * \tparam  abscissas_sizeT         size of abszissas array.
 * \tparam  abscissa_float_t         Floating type for the abscissa values.
 * \tparam  shape_t         type of 1D shapefunctions.
 * \tparam  polynomial_degree         polynomial degree of 1D shape funcion space.
 * \tparam  polynomial_degree         dimension of function that is evaluated.
 * \param   abscissas Abscissas of the supporting points.
 * \param   coefficients flat array of the coefficients of the linear combination.
 **********************************************************************************************/
template < typename return_t, typename abscissa_float_t, typename shape_t,  unsigned int dimT,
    unsigned int abscissas_sizeT, unsigned int polynomial_degree, unsigned int coefficient_dimension>
static inline std::array<std::array<return_t, (unsigned int)std::pow(abscissas_sizeT, dimT)>, coefficient_dimension>
    sum_all_in_tensorial_points_with_coefficients(const std::array<abscissa_float_t, abscissas_sizeT>  abscissas,
                              std::array<return_t, coefficient_dimension*(unsigned int)std::pow(polynomial_degree+1, dimT)> coefficients )
{
  const unsigned int tensorial_shapefunction_count = std::pow(polynomial_degree+1,dimT);
  const unsigned int tensorial_point_count = std::pow(abscissas_sizeT, dimT);
  std::array<std::array<return_t,tensorial_point_count>, coefficient_dimension> return_values;

  for(unsigned int point = 0; point < tensorial_point_count; ++point)
    for(unsigned int d = 0; d < coefficient_dimension; ++d)
      return_values[d][point] = 0;

  std::array<unsigned int, polynomial_degree + 1> indices_1D; //create indices of all shapefunctions
  for(unsigned int i = 0; i < polynomial_degree + 1; ++i)
    indices_1D[i] = i;

  std::array<std::array<return_t, abscissas_sizeT>, polynomial_degree+1> values_1D
  = shape_fct_eval<return_t, shape_t>(indices_1D, abscissas);

  return_t function_value;
  for(unsigned int function = 0; function < tensorial_shapefunction_count; ++function) {
    std::array<unsigned int, std::max(dimT,1U)> decomposed_function_index = index_decompose<dimT, polynomial_degree+1>(function);
    for (unsigned int point = 0; point < tensorial_point_count; ++point) {
      std::array<unsigned int, std::max(dimT,1U)> decomposed_point_index = index_decompose<dimT, abscissas_sizeT>(point);

      function_value = 1.;
      for (unsigned int d = 0; d < std::max(dimT,1U); ++d) {
        function_value *= values_1D[decomposed_function_index[d]][decomposed_point_index[d]];
      }
      for (unsigned int d = 0; d < coefficient_dimension; ++d) {
        return_values[d][point] += function_value*coefficients[d*tensorial_shapefunction_count+function];
      }
    }
  }
  return return_values;
}
