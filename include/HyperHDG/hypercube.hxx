#pragma once  // Ensure that file is included only once in a single compilation.

/*!*************************************************************************************************
 * \brief   Helper class containing numbers and functions related to hypercubes.
 *
 * \tparam  dim     Unsigned integer indicating the dimension of the considered hypercube.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int dimT>
struct Hypercube
{
  /*!***********************************************************************************************
   * \brief   Number of faces of the hypercube.
   ************************************************************************************************/
  static constexpr unsigned int n_faces() { return 2 * dimT; }
  /*!***********************************************************************************************
   * \brief   Number of vertices of the hypercube.
   ************************************************************************************************/
  static constexpr unsigned int n_vertices() { return 1 << dimT; }
  /*!***********************************************************************************************
   * \brief   Return \c n to the power \c dim.
   ************************************************************************************************/
  static constexpr unsigned int pow(unsigned int base)
  {
    unsigned int result = 1;
    if constexpr (dimT > 0)
      for (unsigned int i = 0; i < dimT; ++i)
        result *= base;
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Decompose global index of tensorial structure with respect to its dimensions.
   *
   * \param   index         Local tensorial index (e.g. of the shape function or point).
   * \param   range         Range (maximum, excluded) of the 1D indices.
   * \retval  decomposition Array consisting of respective one-dimensional indices.
   ************************************************************************************************/
  static inline std::array<unsigned int, std::max(dimT, 1U)> index_decompose(
    unsigned int index,
    const unsigned int range)
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
  /*!***********************************************************************************************
   * \brief   Generate tensorial point whose entries are elements of an array.
   *
   * \param   index         Local tensorial index (e.g. of the shape function or point).
   * \param   points        Array of possible values of point components.
   * \retval  point         Tensorial point of given index.
   ************************************************************************************************/
  template <typename point_t, typename array_t>
  static constexpr point_t tensorial_pt(__attribute__((unused)) const unsigned int index,
                                        const array_t& points)
  {
    static_assert(point_t::size() == dimT, "Point must have appropriate dimension!");
    point_t point;
    if constexpr (dimT > 0)
    {
      std::array<unsigned int, std::max(dimT, 1U)> decomp = index_decompose(index, array_t::size());
      for (unsigned int dim = 0; dim < dimT; ++dim)
        point[dim] = points[decomp[dim]];
    }
    return point;
  }
};  // end of struct Hypercube
