#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/global_loop/mat_generate.hxx>
#include <HyperHDG/global_loop/prototype.hxx>
#include <HyperHDG/hy_assert.hxx>

#include <algorithm>
#include <cstddef>
#include <vector>

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Optional cache of the condensed operator of a global loop as a sparse matrix.
 *
 * Evaluating \c trace_to_flux matrix-free solves the local problem of every hyperedge in every
 * call. If the cache is enabled, the condensed matrix is assembled once and \c trace_to_flux
 * becomes a sparse matrix--vector product. This pays off if \c trace_to_flux is evaluated many
 * times for the same operator, e.g. by iterative solvers or in time stepping. The results coincide
 * with the matrix-free ones up to rounding errors (the order of summation differs).
 *
 * The matrix is stored in compressed sparse row (CSR) format, where the contributions of all
 * hyperedges to an entry are summed up (in increasing order of the hyperedges, i.e., independently
 * of the number of OpenMP threads). Its memory is about \f$(2 d n)^2\f$ values and column indices
 * per hyperedge for \f$n\f$ degrees of freedom per face of a \f$d\f$-dimensional hyperedge, minus
 * the entries shared by neighboring hyperedges. It is assembled directly into this format, i.e.,
 * no intermediate triplets (as returned by \c trace_to_flux_mat) are stored.
 *
 * The cache is rebuilt if the parameter of \c trace_to_flux (the time for elliptic and parabolic
 * problems, the shift for eigenvalue problems) changes, unless it has been declared to be
 * independent of this parameter. It must be invalidated if the hypergraph changes (refinement).
 *
 * \warning The cache assumes that the condensed operator is determined by the hypergraph, the local
 *          solver, and the parameter of \c trace_to_flux only. This holds for the local solvers of
 *          HyperHDG, whose \c trace_to_flux does not use the data stored at the hyperedges (such
 *          as the solution of the previous time step, which only enters \c residual_flux). A local
 *          solver whose \c trace_to_flux depends on these data (e.g. the linearization of a
 *          nonlinear problem) would silently be evaluated with a stale matrix. In this case, the
 *          cache must not be used (or it has to be reconfigured whenever the data change).
 *
 * \tparam  LargeVecT   Type of large (global) vectors.
 **************************************************************************************************/
template <typename LargeVecT>
class OperatorCache
{
  /*!***********************************************************************************************
   * \brief   Floating type is determined by floating type of large vector's entries.
   ************************************************************************************************/
  using dof_value_t = typename LargeVecT::value_type;

 public:
  /*!***********************************************************************************************
   * \brief   Enable or disable the cache (and drop the cached matrix).
   *
   * \param   enable            Use the cached matrix in \c trace_to_flux.
   * \param   factorize_once    Assemble the matrix with the local matrix computed only once per
   *                            combination of face types (cf. \c for_each_local_matrix).
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
    std::vector<std::size_t>().swap(row_ptr_);
    std::vector<unsigned int>().swap(col_ind_);
    std::vector<dof_value_t>().swap(values_);
  }
  /*!***********************************************************************************************
   * \brief   Whether the cache is enabled.
   ************************************************************************************************/
  bool enabled() const { return enabled_; }
  /*!***********************************************************************************************
   * \brief   Evaluate the cached operator, (re)assembling it if necessary.
   *
   * \tparam  hyNode_index_t  Typename of the hypernode index.
   * \param   x_vec           Input vector.
   * \param   param           Parameter of \c trace_to_flux (time or shift).
   * \param   hyper_graph     The hypergraph.
   * \param   local_solver    The local solver.
   * \retval  vec_Ax          Result of the matrix--vector product.
   ************************************************************************************************/
  template <typename hyNode_index_t, typename HyperGraphT, typename LocalSolverT>
  LargeVecT apply(const LargeVecT& x_vec,
                  const dof_value_t param,
                  HyperGraphT& hyper_graph,
                  const LocalSolverT& local_solver)
  {
    if (!valid_ || (!param_independent_ && param != param_))
    {
      assemble<hyNode_index_t>(hyper_graph, local_solver, param);
      param_ = param;
    }
    hy_assert(x_vec.size() + 1 == row_ptr_.size(), "Vector does not fit to cached operator!");

    LargeVecT vec_Ax(x_vec.size(), 0.);
    const long long n_rows = (long long)row_ptr_.size() - 1;
    HYPERHDG_OMP_PARALLEL_FOR
    for (long long row = 0; row < n_rows; ++row)
    {
      dof_value_t sum = 0.;
      for (std::size_t k = row_ptr_[row]; k < row_ptr_[row + 1]; ++k)
        sum += values_[k] * x_vec[col_ind_[k]];
      vec_Ax[row] = sum;
    }
    return vec_Ax;
  }

 private:
  /*!***********************************************************************************************
   * \brief   Assemble the condensed matrix in CSR format (with summed up duplicates).
   *
   * Hypernode \c node has the degrees of freedom \c node * n + \c k for \c k < n (cf.
   * HyperNodeFactory). Thus, all rows of a hypernode have the same columns, namely the degrees of
   * freedom of all hypernodes that share a hyperedge with it (its neighbors, in increasing order).
   ************************************************************************************************/
  template <typename hyNode_index_t, typename HyperGraphT, typename LocalSolverT>
  void assemble(HyperGraphT& hyper_graph, const LocalSolverT& local_solver, const dof_value_t param)
  {
    constexpr unsigned int n_faces = 2 * HyperGraphT::hyEdge_dim();
    constexpr unsigned int n_dofs = LocalSolverT::n_glob_dofs_per_node();
    constexpr unsigned int n_loc = n_faces * n_dofs;
    MatGenerate::check_dof_index_range<n_dofs>(hyper_graph);
    invalidate();  // Free the memory of the old matrix first.

    const std::size_t n_hyEdges = hyper_graph.n_hyEdges();
    const std::size_t n_hyNodes = hyper_graph.hyNode_factory().n_hyNodes();

    // Hypernodes of all hyperedges, and hyperedges of all hypernodes.
    std::vector<hyNode_index_t> edge_nodes(n_hyEdges * n_faces);
    for (std::size_t edge = 0; edge < n_hyEdges; ++edge)
    {
      const SmallVec<n_faces, hyNode_index_t> hyNodes =
        hyper_graph[edge].topology.get_hyNode_indices();
      for (unsigned int face = 0; face < n_faces; ++face)
        edge_nodes[edge * n_faces + face] = hyNodes[face];
    }
    std::vector<std::size_t> node_edges_ptr(n_hyNodes + 1, 0);
    for (const hyNode_index_t node : edge_nodes)
      ++node_edges_ptr[node + 1];
    for (std::size_t node = 0; node < n_hyNodes; ++node)
      node_edges_ptr[node + 1] += node_edges_ptr[node];
    std::vector<std::size_t> node_edges(edge_nodes.size());
    {
      std::vector<std::size_t> next(node_edges_ptr.begin(), node_edges_ptr.end() - 1);
      for (std::size_t k = 0; k < edge_nodes.size(); ++k)
        node_edges[next[edge_nodes[k]]++] = k / n_faces;
    }

    // Neighbors of all hypernodes (sorted, including the hypernode itself).
    std::vector<std::size_t> neighbors_ptr(n_hyNodes + 1, 0);
    std::vector<hyNode_index_t> neighbors;
    {
      std::vector<hyNode_index_t> helper;
      for (std::size_t node = 0; node < n_hyNodes; ++node)
      {
        helper.clear();
        for (std::size_t k = node_edges_ptr[node]; k < node_edges_ptr[node + 1]; ++k)
          for (unsigned int face = 0; face < n_faces; ++face)
            helper.push_back(edge_nodes[node_edges[k] * n_faces + face]);
        std::sort(helper.begin(), helper.end());
        helper.erase(std::unique(helper.begin(), helper.end()), helper.end());
        neighbors.insert(neighbors.end(), helper.begin(), helper.end());
        neighbors_ptr[node + 1] = neighbors.size();
      }
    }
    std::vector<std::size_t>().swap(node_edges);
    std::vector<std::size_t>().swap(node_edges_ptr);
    std::vector<hyNode_index_t>().swap(edge_nodes);

    // CSR structure: Row node * n_dofs + k has the columns neighbor * n_dofs + l of all neighbors.
    row_ptr_.resize(n_hyNodes * n_dofs + 1);
    for (std::size_t node = 0; node < n_hyNodes; ++node)
    {
      const std::size_t n_neighbors = neighbors_ptr[node + 1] - neighbors_ptr[node];
      for (unsigned int k = 0; k < n_dofs; ++k)
        row_ptr_[node * n_dofs + k] = n_dofs * (n_dofs * neighbors_ptr[node] + k * n_neighbors);
    }
    row_ptr_[n_hyNodes * n_dofs] = n_dofs * n_dofs * neighbors_ptr[n_hyNodes];
    col_ind_.resize(row_ptr_.back());
    values_.assign(row_ptr_.back(), 0.);
    for (std::size_t node = 0; node < n_hyNodes; ++node)
      for (unsigned int k = 0; k < n_dofs; ++k)
      {
        std::size_t index = row_ptr_[node * n_dofs + k];
        for (std::size_t p = neighbors_ptr[node]; p < neighbors_ptr[node + 1]; ++p)
          for (unsigned int l = 0; l < n_dofs; ++l)
            col_ind_[index++] = neighbors[p] * n_dofs + l;
      }

    // Sum up the local matrices (serially in increasing order of the hyperedges).
    MatGenerate::for_each_local_matrix<hyNode_index_t, false>(
      hyper_graph, local_solver, param, factorize_once_,
      [&](const std::size_t, const auto& hyNodes, const dof_value_t* loc_mat)
      {
        for (unsigned int i = 0; i < n_faces; ++i)
        {
          const auto first = neighbors.begin() + neighbors_ptr[hyNodes[i]],
                     last = neighbors.begin() + neighbors_ptr[hyNodes[i] + 1];
          for (unsigned int j = 0; j < n_faces; ++j)
          {
            const std::size_t block = std::lower_bound(first, last, hyNodes[j]) - first;
            // Entry (k, l) of the block is at row_start + k * row_length + l.
            dof_value_t* const block_start =
              &values_[row_ptr_[hyNodes[i] * n_dofs] + block * n_dofs];
            const std::size_t row_length = n_dofs * (last - first);
            for (unsigned int l = 0; l < n_dofs; ++l)
            {
              const dof_value_t* column = &loc_mat[(j * n_dofs + l) * n_loc + i * n_dofs];
              for (unsigned int k = 0; k < n_dofs; ++k)
                block_start[k * row_length + l] += column[k];
            }
          }
        }
      });
    valid_ = true;
  }

  bool enabled_ = false, factorize_once_ = false, param_independent_ = false, valid_ = false;
  dof_value_t param_ = 0.;
  std::vector<std::size_t> row_ptr_;
  std::vector<unsigned int> col_ind_;
  std::vector<dof_value_t> values_;
};

}  // end of namespace GlobalLoop
