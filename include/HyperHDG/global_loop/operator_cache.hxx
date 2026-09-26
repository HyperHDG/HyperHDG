#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/global_loop/mat_generate.hxx>
#include <HyperHDG/global_loop/prototype.hxx>
#include <HyperHDG/hy_assert.hxx>

#include <vector>

/*!*************************************************************************************************
 * \brief   Optional cache of the condensed operator of a global loop as a sparse matrix.
 *
 * Evaluating \c trace_to_flux matrix-free solves the local problem of every hyperedge in every
 * call. If the cache is enabled, the condensed matrix is assembled once (via \c trace_to_flux_mat)
 * and \c trace_to_flux becomes a sparse matrix--vector product. This pays off if \c trace_to_flux
 *is evaluated many times for the same operator, e.g. by iterative solvers or in time stepping. The
 * results coincide with the matrix-free ones up to rounding errors (the order of summation
 * differs). The price is the memory for the matrix, i.e., \f$(2 d n)^2\f$ values per hyperedge for
 * \f$n\f$ degrees of freedom per face of a \f$d\f$-dimensional hyperedge.
 *
 * The cache is rebuilt if the parameter of \c trace_to_flux (the time for elliptic and parabolic
 * problems, the shift for eigenvalue problems) changes, unless it has been declared to be
 * independent of this parameter. It must be invalidated if the hypergraph changes (refinement).
 *
 * \tparam  LargeVecT   Type of large (global) vectors.
 * \tparam  dof_value_t Floating point type of the degrees of freedom.
 **************************************************************************************************/
template <typename LargeVecT, typename dof_value_t>
class OperatorCache
{
 public:
  /*!***********************************************************************************************
   * \brief   Enable or disable the cache (and drop the cached matrix).
   *
   * \param   enable            Use the cached matrix in \c trace_to_flux.
   * \param   factorize_once    Assemble the matrix with the local matrix computed only once per
   *                            combination of face types (cf. \c prototype_mat_generate_local).
   * \param   param_independent The operator does not depend on the parameter (e.g. time) of
   *                            \c trace_to_flux, i.e., the cached matrix is never rebuilt.
   ************************************************************************************************/
  void configure(const bool enable, const bool factorize_once, const bool param_independent)
  {
    enabled_ = enable;
    factorize_once_ = factorize_once;
    param_independent_ = param_independent;
    invalidate();
  }
  /*!***********************************************************************************************
   * \brief   Drop the cached matrix (e.g. after a refinement of the hypergraph).
   ************************************************************************************************/
  void invalidate()
  {
    valid_ = false;
    row_ptr_.clear();
    col_ind_.clear();
    values_.clear();
  }
  /*!***********************************************************************************************
   * \brief   Whether the cache is enabled.
   ************************************************************************************************/
  bool enabled() const { return enabled_; }
  /*!***********************************************************************************************
   * \brief   Whether the local matrices are only computed once per combination of face types.
   ************************************************************************************************/
  bool factorize_once() const { return factorize_once_; }
  /*!***********************************************************************************************
   * \brief   Evaluate the cached operator, (re)assembling it via \c assemble(param) if necessary.
   *
   * \param   x_vec     Input vector.
   * \param   param     Parameter of \c trace_to_flux (time or shift).
   * \param   assemble  Callable that returns the \c sparse_mat of the operator for a parameter.
   * \retval  vec_Ax    Result of the matrix--vector product.
   ************************************************************************************************/
  template <typename AssembleT>
  LargeVecT apply(const LargeVecT& x_vec, const dof_value_t param, AssembleT&& assemble)
  {
    if (!valid_ || (!param_independent_ && param != param_))
    {
      auto triplets = assemble(param);
      build(triplets, x_vec.size());
      param_ = param;
    }
    hy_assert(x_vec.size() + 1 == row_ptr_.size(), "Vector does not fit to cached operator!");

    LargeVecT vec_Ax(x_vec.size(), 0.);
    const long long n_rows = (long long)row_ptr_.size() - 1;
    HYPERHDG_OMP_PARALLEL_FOR
    for (long long row = 0; row < n_rows; ++row)
    {
      dof_value_t sum = 0.;
      for (unsigned int k = row_ptr_[row]; k < row_ptr_[row + 1]; ++k)
        sum += values_[k] * x_vec[col_ind_[k]];
      vec_Ax[row] = sum;
    }
    return vec_Ax;
  }

 private:
  /*!***********************************************************************************************
   * \brief   Convert triplets into compressed sparse row format (duplicates are kept).
   ************************************************************************************************/
  template <typename SparseMatT>
  void build(SparseMatT& triplets, const unsigned int n_rows)
  {
    const auto& rows = triplets.get_rows();
    const auto& cols = triplets.get_cols();
    const auto& vals = triplets.get_values();
    row_ptr_.assign(n_rows + 1, 0);
    for (unsigned int k = 0; k < rows.size(); ++k)
      ++row_ptr_[rows[k] + 1];
    for (unsigned int row = 0; row < n_rows; ++row)
      row_ptr_[row + 1] += row_ptr_[row];
    col_ind_.resize(rows.size());
    values_.resize(rows.size());
    std::vector<unsigned int> next(row_ptr_.begin(), row_ptr_.end() - 1);
    for (unsigned int k = 0; k < rows.size(); ++k)
    {
      col_ind_[next[rows[k]]] = cols[k];
      values_[next[rows[k]]++] = vals[k];
    }
    valid_ = true;
  }

  bool enabled_ = false, factorize_once_ = false, param_independent_ = false, valid_ = false;
  dof_value_t param_ = 0.;
  std::vector<unsigned int> row_ptr_, col_ind_;
  std::vector<dof_value_t> values_;
};
