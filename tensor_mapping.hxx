#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/HyAssert.hxx>
#include <HyperHDG/SmallVec.hxx>
#include <HyperHDG/Hypercube.hxx>

#include <array>

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
    Point<rdim> operator() (std::array<unsigned int,domain_dimension()> indices) const;

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
  Mapping<rdim, BaseTensor>::operator() (std::array<unsigned int,domain_dimension()> indices) const
  {
    Point<rdim> result;
    for (unsigned int i=0;i<rdim;++i)
      result[i] = data[i](indices);
  }

  /**
   * \brief Mapping by tensor product of linear polynomials
   *
   * \todo Switch the order of rdim and ddim?
   *
   * \tparam rdim: the dimension of the space mapped into
   * \tparam ddim: the dimension of the domain
   * \tparam npts: the number of points in each direction
   * \tparam T: the number type
   */
  template <int rdim, int ddim, std::size_t npts, typename T=double>
  class MappingMultilinear
  {
  public:
    /// The dimension of the domain of the mapping
    static constexpr unsigned int domain_dimension()
    { return ddim; }

    /// The dimension of the range of the mapping
    static constexpr unsigned int range_dimension()
    { return rdim; }

    /// The number of points in coordinate direction `d`
    static constexpr unsigned int size(unsigned int d)
    { return npts; }

    /// The number of points total
    static constexpr unsigned int size()
    { return Hypercube<ddim>::pow(npts); }

    /**
     * Constructor, obtaining the corner points and the one dimensional evaluation points
     */
    template <typename T2>
    MappingMultilinear(const std::array<Point<rdim>, Hypercube<ddim>::n_vertices()>& vertices,
		       const std::array<T2, npts>& points);
    
    /// Access to a mapped point with tensor coordinates `indices`
    Point<rdim> operator() (const std::array<unsigned int, ddim>& indices) const;

    /// Access points in lexicographic order, first index fastest
    Point<rdim> lexicographic (unsigned int index) const;
    
  private:
    /// Array of the corner points of the cell
    std::array<Point<rdim>, Hypercube<ddim>::n_vertices()> vertices;
    /// Array of one-dimensional quadrature points
    std::array<T, npts> points_1d;
  };

  template <int rdim, int ddim, std::size_t npts, typename T>
  template <typename T2>
  MappingMultilinear<rdim,ddim,npts,T>::MappingMultilinear(const std::array<Point<rdim>,Hypercube<ddim>::n_vertices()>& vertices,
							   const std::array<T2, npts>& points)
    : vertices(vertices), points_1d(points)
  {}
  
  template <int rdim, int ddim, std::size_t npts, typename T>
  Point<rdim>
  MappingMultilinear<rdim,ddim,npts,T>::operator() (const std::array<unsigned int, ddim>& ind) const
  {
    static_assert(ind.size() == ddim);
    Point<rdim> result;

    for (unsigned int v=0;v<vertices.size();++v)
      {
	T shape_value = 1.;
	for (unsigned int d=0;d<ddim;++d)
	  shape_value *= ((v & (1<<d)) == 0) ? (1.-points_1d[ind[d]]) : points_1d[ind[d]];
	
	for (unsigned int r=0;r<rdim;++r)
	  result[r] += shape_value * vertices[v][r];
      }
    return result;
  }

  template <int rdim, int ddim, std::size_t npts, typename T>
  Point<rdim>
  MappingMultilinear<rdim,ddim,npts,T>::lexicographic (unsigned i) const
  {
    std::array<unsigned int, ddim> ind;
    for (unsigned int d=0;d<ddim;++d)
      {
	ind[d] = i%npts;
	i /= npts;
      }
    return (*this)(ind);
  }
}
