#pragma once // Ensure that file is included only once in a single compilation.

#include <initializer_list>
#include <array>

/**
 * \brief The tensor objects available in this library.
 *
 * A tensor in this library is an object with at least the signature
 * \code
 * template <typename T>
 * struct Tensor
 * {
 *   typedef T value_type;
 *   static unsigned int order() const;
 *   static unsigned int dimension(unsigned int d) const;
 *   const T operator() (std::initializer_list<unsigned int> indices) const;
 * };
 * \endcode
 *
 * Here, `dimension()` is the dimension of the tensor in direction
 * `d`, that is, the length of fibers in this direction. The argument
 * `indices` to `operator()` is a list of `order()` values of type
 * `unsigned int`, denoting the position of the entry requested.
 *
 * Isotropic tensors implement additionally a function `dimension()`
 * without argument.
 */
namespace Tensor
{
  /**
   * \brief Implementation of a rank one tensor of uniform dimension
   *
   * \tparam order The order of the tensor
   * \tparam dim The dimension of every factor
   * \tparam T A numerical type, the value type of the tensor
   *
   * \author Guido Kanschat, Andreas Rupp, 2020
   */
  template <int ord, int dim, typename T=double>
  class Rank1
  {
  public:
    /// Access to an element of the tensor conforming to the general tensor interface
    const T operator() (std::initializer_list<unsigned int> indices) const;

    /// Access to a factor of the tensor
    std::array<T, dim>& factor(unsigned int i);

    /// The value type of the tensor, namely #T
    typedef T value_type;

    /// The order of the tensor
    static constexpr unsigned int order()
    { return ord; }

    /// The dimension of the factors conforming to general tensor interface
    static constexpr unsigned int dimension(unsigned int)
    { return dim; }

    /// The dimension of the factors conforming to interface for isotropic tensors
    static constexpr unsigned int dimension()
    { return dim; }
  private:
    /// The order one factors of the tensor
    std::array<std::array<T, dim>, ord> data;
  };

  template <int ord, int dim, typename T>
  const T
  Rank1<ord, dim, T>::operator() (std::initializer_list<unsigned int> indices) const
  {
    T result = 1.;
    unsigned int d = 0;
    for (unsigned int i : indices)
      {
	result *= data[d][i];
	++d;
      }
    return result;
  }

  template <int ord, int dim, typename T>
  std::array<T, dim>&
  Rank1<ord, dim, T>::factor(unsigned int i)
  {
    return data[i];
  }
}
