#ifndef TENSOR_RANK1
#define TENSOR_RANK1

#include <initializer_list>
#include <array>

/**
 * \brief The tensor objects available in this library.
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
  template <int order, int dim, typename T=double>
  class Rank1
  {
  public:
    /// Access to an element of the tensor
    const T operator() (std::initializer_list<unsigned int> indices) const;

    /// Access to a factor of the tensor
    std::array<T, dim>& factor(unsigned int i);

    /// The value type of the tensor, namely #T
    typedef T value_type;
  private:
    /// The order one factors of the tensor
    std::array<std::array<T, dim>, order> data;
  };

  template <int order, int dim, typename T>
  const T
  Rank1<order, dim, T>::operator() (std::initializer_list<unsigned int> indices) const
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

  template <int order, int dim, typename T>
  std::array<T, dim>&
  Rank1<order, dim, T>::factor(unsigned int i)
  {
    return data[i];
  }
}

#endif
