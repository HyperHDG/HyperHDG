#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <vector>

#if __has_include(<cholmod.h>)
#include <cholmod.h>
#elif __has_include(<suitesparse/cholmod.h>)
#include <suitesparse/cholmod.h>
#endif

namespace tests
{
/*!*************************************************************************************************
 * \brief   The assembled condensed system of a GlobalLoop::Generic, restricted to its free dofs.
 *
 * Constrained dofs are zero in the solution (Dirichlet data enter through the load), so their rows
 * and columns are dropped; the free block of LDG-H diffusion is symmetric. CSR with duplicates of
 * the per-hyperedge triplets summed; multiply() is y = A x on the free dofs, expand() maps a free
 * vector back to all dofs.
 **************************************************************************************************/
struct FreeSystem
{
  std::vector<unsigned int> free;      // free dof -> global dof
  std::vector<long long> index;        // global dof -> free dof, or -1
  std::vector<std::size_t> row_begin;  // CSR
  std::vector<unsigned int> column;
  std::vector<double> value;
  std::vector<double> rhs;  // free part of -residual(0)

  template <class LoopT>
  explicit FreeSystem(LoopT& loop)
  {
    const unsigned int n = loop.size_of_system();
    index.assign(n, 0);
    for (unsigned int i : loop.constraints())
      index[i] = -1;
    for (unsigned int i = 0; i < n; ++i)
      if (index[i] >= 0)
      {
        index[i] = free.size();
        free.push_back(i);
      }
    const std::size_t m = loop.n_triplets(false);
    std::vector<unsigned int> rows(m), cols(m);
    std::vector<double> vals(m);
    loop.assemble(typename LoopT::Args(), false, rows, cols, vals);
    // CSR of the free block: count, then fill, then sum duplicates per row
    const std::size_t n_free = free.size();
    row_begin.assign(n_free + 1, 0);
    for (std::size_t k = 0; k < m; ++k)
      if (index[rows[k]] >= 0 && index[cols[k]] >= 0)
        ++row_begin[index[rows[k]] + 1];
    std::partial_sum(row_begin.begin(), row_begin.end(), row_begin.begin());
    std::vector<std::size_t> next(row_begin.begin(), row_begin.end() - 1);
    std::vector<std::pair<unsigned int, double>> entries(row_begin.back());
    for (std::size_t k = 0; k < m; ++k)
      if (index[rows[k]] >= 0 && index[cols[k]] >= 0)
        entries[next[index[rows[k]]]++] = {(unsigned int)index[cols[k]], vals[k]};
    std::vector<std::size_t> compressed(n_free + 1, 0);
    for (std::size_t r = 0; r < n_free; ++r)
    {
      std::sort(entries.begin() + row_begin[r], entries.begin() + row_begin[r + 1]);
      for (std::size_t k = row_begin[r]; k < row_begin[r + 1]; ++k)
        if (k > row_begin[r] && column.back() == entries[k].first)  // same row, same column
          value.back() += entries[k].second;
        else
        {
          column.push_back(entries[k].first);
          value.push_back(entries[k].second);
        }
      compressed[r + 1] = column.size();
    }
    row_begin = compressed;
    std::vector<double> x(n, 0.), r(n);
    loop.residual(x, r);
    rhs.resize(n_free);
    for (std::size_t i = 0; i < n_free; ++i)
      rhs[i] = -r[free[i]];
  }
  std::size_t size() const { return free.size(); }
  std::size_t nnz() const { return value.size(); }
  void multiply(const std::vector<double>& x, std::vector<double>& y) const
  {
    for (std::size_t r = 0; r < size(); ++r)
    {
      double sum = 0.;
      for (std::size_t k = row_begin[r]; k < row_begin[r + 1]; ++k)
        sum += value[k] * x[column[k]];
      y[r] = sum;
    }
  }
  std::vector<double> expand(const std::vector<double>& x_free) const
  {
    std::vector<double> x(index.size(), 0.);
    for (std::size_t i = 0; i < size(); ++i)
      x[free[i]] = x_free[i];
    return x;
  }
};

#if defined(CHOLMOD_MAIN_VERSION)
/*!*************************************************************************************************
 * \brief   Sparse Cholesky of a FreeSystem with CHOLMOD (analyze, factorize, solve), RAII.
 *
 * The condensed operator of LDG-H diffusion is negative definite with Generic's sign convention;
 * the factorization is of -A (sign = -1) then, solve returns the solution of A x = b either way.
 **************************************************************************************************/
class Cholmod
{
  cholmod_common common_;
  cholmod_sparse* matrix_ = nullptr;
  cholmod_factor* factor_ = nullptr;
  double sign_ = 1.;

 public:
  explicit Cholmod(const FreeSystem& system)
  {
    cholmod_start(&common_);
    const std::size_t n = system.size();
    sign_ = system.value[0] < 0. ? -1. : 1.;  // the first diagonal entry is the first of row 0
    // upper triangle in CSC: the CSR rows of a symmetric matrix are its CSC columns
    std::size_t nnz_upper = 0;
    for (std::size_t r = 0; r < n; ++r)
      for (std::size_t k = system.row_begin[r]; k < system.row_begin[r + 1]; ++k)
        nnz_upper += system.column[k] <= r;
    matrix_ = cholmod_allocate_sparse(n, n, nnz_upper, 1, 1, 1, CHOLMOD_REAL, &common_);
    auto* p = static_cast<int*>(matrix_->p);
    auto* i = static_cast<int*>(matrix_->i);
    auto* x = static_cast<double*>(matrix_->x);
    std::size_t pos = 0;
    for (std::size_t c = 0; c < n; ++c)
    {
      p[c] = pos;
      for (std::size_t k = system.row_begin[c]; k < system.row_begin[c + 1]; ++k)
        if (system.column[k] <= c)
        {
          i[pos] = system.column[k];
          x[pos++] = sign_ * system.value[k];
        }
    }
    p[n] = pos;
  }
  ~Cholmod()
  {
    cholmod_free_factor(&factor_, &common_);
    cholmod_free_sparse(&matrix_, &common_);
    cholmod_finish(&common_);
  }
  Cholmod(const Cholmod&) = delete;
  Cholmod& operator=(const Cholmod&) = delete;

  void analyze()
  {
    cholmod_free_factor(&factor_, &common_);
    factor_ = cholmod_analyze(matrix_, &common_);
  }
  bool factorize()
  {
    cholmod_factorize(matrix_, factor_, &common_);
    return common_.status == CHOLMOD_OK && factor_->minor == factor_->n;
  }
  std::vector<double> solve(const std::vector<double>& b)
  {
    cholmod_dense* rhs = cholmod_allocate_dense(b.size(), 1, b.size(), CHOLMOD_REAL, &common_);
    for (std::size_t k = 0; k < b.size(); ++k)
      static_cast<double*>(rhs->x)[k] = sign_ * b[k];
    cholmod_dense* x = cholmod_solve(CHOLMOD_A, factor_, rhs, &common_);
    std::vector<double> result(static_cast<double*>(x->x), static_cast<double*>(x->x) + b.size());
    cholmod_free_dense(&x, &common_);
    cholmod_free_dense(&rhs, &common_);
    return result;
  }
  double factor_nnz() const { return common_.lnz; }
  double factor_flops() const { return common_.fl; }
};
#endif

}  // end of namespace tests
