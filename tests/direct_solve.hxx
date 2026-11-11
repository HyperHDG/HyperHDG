#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/wrapper/lapack.hxx>

#include <cmath>
#include <cstddef>
#include <vector>

namespace tests
{
/*!*************************************************************************************************
 * \brief   Dense LU solve of the assembled matrix of a GlobalLoop::Generic (with a unit
 *          diagonal on constrained dofs); throws if the matrix is singular. For small meshes only.
 **************************************************************************************************/
template <class LoopT>
std::vector<typename LoopT::value_type> direct_solve(LoopT& loop,
                                                     std::vector<typename LoopT::value_type> rhs)
{
  using real = typename LoopT::value_type;
  const unsigned int size = loop.size_of_system();
  const std::size_t n_triplets = loop.n_triplets(true);
  std::vector<unsigned int> rows(n_triplets), cols(n_triplets);
  std::vector<real> vals(n_triplets), dense((std::size_t)size * size, 0.);
  loop.assemble(typename LoopT::Args(), true, rows, cols, vals);
  for (std::size_t k = 0; k < n_triplets; ++k)
    dense[(std::size_t)cols[k] * size + rows[k]] += vals[k];  // column-major, duplicates summed
  Wrapper::lapack_solve(size, 1, dense.data(), rhs.data());
  return rhs;
}
/*!*************************************************************************************************
 * \brief   Whether the assembled matrix of a GlobalLoop::Generic (with a unit diagonal on
 *          constrained dofs) is positive definite: its dense Cholesky factorization exists. For
 *          small meshes only.
 **************************************************************************************************/
template <class LoopT>
bool positive_definite(LoopT& loop)
{
  using real = typename LoopT::value_type;
  const std::size_t size = loop.size_of_system(), n_triplets = loop.n_triplets(true);
  std::vector<unsigned int> rows(n_triplets), cols(n_triplets);
  std::vector<real> vals(n_triplets), a(size * size, 0.);
  loop.assemble(typename LoopT::Args(), true, rows, cols, vals);
  for (std::size_t k = 0; k < n_triplets; ++k)
    a[rows[k] * size + cols[k]] += vals[k];
  for (std::size_t j = 0; j < size; ++j)  // the lower triangle becomes the factor
  {
    for (std::size_t k = 0; k < j; ++k)
      a[j * size + j] -= a[j * size + k] * a[j * size + k];
    if (!(a[j * size + j] > 0.))
      return false;
    a[j * size + j] = std::sqrt(a[j * size + j]);
    for (std::size_t i = j + 1; i < size; ++i)
    {
      for (std::size_t k = 0; k < j; ++k)
        a[i * size + j] -= a[i * size + k] * a[j * size + k];
      a[i * size + j] /= a[j * size + j];
    }
  }
  return true;
}
}  // end of namespace tests
