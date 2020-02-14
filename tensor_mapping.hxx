#ifndef TENSOR_GEOMETRY_HXX
#define TENSOR_GEOMETRY_HXX

#include <Point.hxx>

namespace Tensor
{
  /**
   * \brief The mapping of tensor products of point sets.
   *
   * This class represents a mapping from a reference hypercube of
   * dimension domain_dimension() to a hypersurface in space of
   * dimension range_dimension().
   *
   * The reference hypercube is discretized by the tensor product of a
   * one-dimensional point set. For each of the components in range
   * space, a tensor of the mapped values is stored, where the actual
   * format of this tensor is open.
   *
   * The implementation is through an `std::array` of #BaseTensor
   * objects. In order to be able to fill and modify these objects,
   * there is an access function to the individual base
   * tensors. Access to the mapped points should not go through this
   * access function, though.
   *
   * \tparam rdim Dimension of the range of this mapping.
   *
   * \tparam BaseTensor The underlying tensor implementation. The
   * order of this tensor is the dimension of the domain of the
   * mapping, while its dimensions are the numbers of mapped points in
   * each coordinate direction.
   */
  template<int rdim, class BaseTensor>
  class Mapping
  {
    /// The array of base tensors storing the data
    std::array<BaseTensor, rdim> data;
    
  public:
    /// The dimension of the domain of the mapping
    static constexpr unsigned int domain_dimension()
    { return BaseTensor::order(); }

    /// The dimension of the range of the mapping
    static constexpr unsigned int range_dimension()
    { return rdim; }

    /// The number of points in coordinate direction `d`
    static constexpr unsigned int size(unsigned int d)
    { return BaseTensor::dimension(d); }

    /// Access to a mapped point with tensor coordinates `indices`
    Point<rdim> operator() (std::initializer_list<unsigned int> indices) const;

    /**
     * \brief Access to the base tensor.
     *
     * Use this function only for modifying the contents. Access the mapped points always by operator()!
     */
    BaseTensor& base_tensor(unsigned int i)
    { return data[i]; }
};

  template <int rdim, class BaseTensor>
  Point<rdim>
  Mapping<rdim, BaseTensor>::operator() (std::initializer_list<unsigned int> indices) const
  {
    Point<rdim> result;
    for (unsigned int i=0;i<rdim;++i)
      result[i] = data[i](indices);
  }
}

#endif
