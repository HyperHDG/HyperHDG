#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/global_loop/prototype.hxx>
#include <HyperHDG/hy_assert.hxx>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

/*!*************************************************************************************************
 * \brief   True if HyperHDG is compiled in debug mode (i.e. without \c NDEBUG).
 **************************************************************************************************/
#ifdef NDEBUG
inline constexpr bool hy_debug_mode = false;
#else
inline constexpr bool hy_debug_mode = true;
#endif

/*!*************************************************************************************************
 * \brief   Parallelize the following loop with OpenMP (if HyperHDG is compiled with OpenMP).
 **************************************************************************************************/
#ifdef _OPENMP
#define HYPERHDG_OMP_PARALLEL_FOR _Pragma("omp parallel for schedule(static)")
#define HYPERHDG_OMP_CRITICAL _Pragma("omp critical(hyperhdg_exceptions)")
#else
#define HYPERHDG_OMP_PARALLEL_FOR
#define HYPERHDG_OMP_CRITICAL
#endif

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Forward exceptions thrown in the body of a parallel loop to the calling thread.
 *
 * An exception that leaves an OpenMP parallel region terminates the program. Therefore, the body of
 * a parallel loop is executed via \c run, which catches all exceptions. Once an exception has been
 * caught, the remaining iterations are skipped, and \c rethrow (called after the loop) rethrows the
 * first caught exception. Without OpenMP, this results in the same behavior as a serial loop.
 **************************************************************************************************/
class ParallelExceptions
{
 public:
  /*!***********************************************************************************************
   * \brief   Execute \c body unless an exception has been thrown before, and catch exceptions.
   ************************************************************************************************/
  template <typename BodyT>
  void run(BodyT&& body) noexcept
  {
    if (failed_.load(std::memory_order_relaxed))
      return;
    try
    {
      body();
    }
    catch (...)
    {
      HYPERHDG_OMP_CRITICAL
      {
        if (!exception_)
          exception_ = std::current_exception();
      }
      failed_.store(true, std::memory_order_relaxed);
    }
  }
  /*!***********************************************************************************************
   * \brief   Rethrow the first caught exception (if any).
   ************************************************************************************************/
  void rethrow()
  {
    if (exception_)
      std::rethrow_exception(exception_);
  }

 private:
  std::atomic<bool> failed_{false};
  std::exception_ptr exception_;
};

/*!*************************************************************************************************
 * \brief   Number of threads that a parallel loop uses.
 **************************************************************************************************/
inline unsigned int hyperhdg_n_threads()
{
#ifdef _OPENMP
  return omp_get_max_threads();
#else
  return 1;
#endif
}

namespace MatGenerate
{
/*!*************************************************************************************************
 * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
 **************************************************************************************************/
HAS_MEMBER_FUNCTION(trace_to_flux, has_trace_to_flux);
/*!*************************************************************************************************
 * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
 **************************************************************************************************/
HAS_MEMBER_FUNCTION(trace_to_flux_mat, has_trace_to_flux_mat);

/*!*************************************************************************************************
 * \brief   Throw if the degrees of freedom of a hypergraph cannot be indexed by \c unsigned \c int.
 *
 * The global loops index degrees of freedom by \c unsigned \c int. Without this check, too large
 * hypergraphs would result in silently wrapped indices.
 **************************************************************************************************/
template <unsigned int n_dofs_per_node, typename HyperGraphT>
void check_dof_index_range(const HyperGraphT& hyper_graph)
{
  const std::size_t n_dofs =
    (std::size_t)hyper_graph.hyNode_factory().n_hyNodes() * n_dofs_per_node;
  if (n_dofs > std::numeric_limits<unsigned int>::max())
  {
    std::stringstream msg;
    msg << "The hypergraph has " << n_dofs << " degrees of freedom, which exceeds the maximum "
        << "index " << std::numeric_limits<unsigned int>::max() << " of the global vectors.";
    throw std::length_error(msg.str());
  }
}

/*!*************************************************************************************************
 * \brief   Local condensed matrix of a hyperedge, i.e., the matrix representation of trace_to_flux.
 *
 * If the local solver provides \c trace_to_flux_mat(hyper_edge, time), which returns the local
 * condensed matrix as a \c SmallSquareMat, this function is used. Local solvers can implement it by
 * factorizing their local matrix once and reusing the factorization for all unit vectors (cf.
 * local_solver/trace_to_flux_mat.hxx). Otherwise, the local condensed matrix is obtained by
 * applying \c trace_to_flux to all local unit vectors.
 *
 * \param   local_solver  The local solver.
 * \param   hyper_edge    The hyperedge.
 * \param   time          Parameter of \c trace_to_flux (time or shift).
 * \param   loc_mat       Array of \f$n^2\f$ entries that is filled with the local matrix in
 *                        column-major order, i.e., entry (row, col) is \c loc_mat[col * n + row].
 **************************************************************************************************/
template <typename LocalSolverT, typename hyEdgeT, typename dof_value_t>
void local_condensed_matrix(const LocalSolverT& local_solver,
                            hyEdgeT& hyper_edge,
                            const dof_value_t time,
                            dof_value_t* loc_mat)
{
  constexpr unsigned int hyEdge_dim = LocalSolverT::hyEdge_dim();
  constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();
  constexpr unsigned int n_loc = 2 * hyEdge_dim * n_dofs_per_node;
  using dofs_t = std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>;

  if constexpr (has_trace_to_flux_mat<LocalSolverT, SmallSquareMat<n_loc, dof_value_t>(
                                                      hyEdgeT&, dof_value_t)>::value)
  {
    const SmallSquareMat<n_loc, dof_value_t> helper =
      local_solver.trace_to_flux_mat(hyper_edge, time);
    for (unsigned int col = 0; col < n_loc; ++col)
      for (unsigned int row = 0; row < n_loc; ++row)
        loc_mat[col * n_loc + row] = helper(row, col);
  }
  else
  {
    dofs_t dofs_old, dofs_new;
    for (unsigned int node_j = 0; node_j < 2 * hyEdge_dim; ++node_j)
      for (unsigned int dof_j = 0; dof_j < n_dofs_per_node; ++dof_j)
      {
        for (unsigned int node = 0; node < 2 * hyEdge_dim; ++node)
        {
          dofs_old[node].fill(0.);
          dofs_new[node].fill(0.);
        }
        dofs_old[node_j][dof_j] = 1.;
        if constexpr (has_trace_to_flux<LocalSolverT,
                                        dofs_t&(dofs_t&, dofs_t&, dof_value_t)>::value)
          local_solver.trace_to_flux(dofs_old, dofs_new, time);
        else if constexpr (has_trace_to_flux<LocalSolverT, dofs_t&(dofs_t&, dofs_t&, hyEdgeT&,
                                                                   dof_value_t)>::value)
          local_solver.trace_to_flux(dofs_old, dofs_new, hyper_edge, time);
        else
          hy_assert(false, "Function seems not to be implemented!");

        const unsigned int col = node_j * n_dofs_per_node + dof_j;
        for (unsigned int node_i = 0; node_i < 2 * hyEdge_dim; ++node_i)
          for (unsigned int dof_i = 0; dof_i < n_dofs_per_node; ++dof_i)
            loc_mat[col * n_loc + node_i * n_dofs_per_node + dof_i] = dofs_new[node_i][dof_i];
      }
  }
}

/*!*************************************************************************************************
 * \brief   Throw if a local matrix differs from the local matrix reused by \c factorize_once.
 **************************************************************************************************/
template <typename dof_value_t>
void check_reused_local_matrix(const dof_value_t* computed,
                               const dof_value_t* reused,
                               const std::size_t size,
                               const std::size_t edge)
{
  dof_value_t max_entry = 0., max_diff = 0.;
  for (std::size_t k = 0; k < size; ++k)
  {
    max_entry = std::max(max_entry, (dof_value_t)std::abs(reused[k]));
    max_diff = std::max(max_diff, (dof_value_t)std::abs(computed[k] - reused[k]));
  }
  const dof_value_t tolerance =
    std::max((dof_value_t)1e-10, (dof_value_t)1e3 * std::numeric_limits<dof_value_t>::epsilon());
  if (!(max_diff <= tolerance * max_entry))  // Also true if max_diff is NaN.
  {
    std::stringstream msg;
    msg << "Option factorize_once requires all hyperedges with the same face types to have the "
        << "same local matrix (e.g. uniform meshes and coefficients that are constant in space), "
        << "but the local matrix of hyperedge " << edge << " differs by " << max_diff
        << " (maximum entry " << max_entry << "). Use factorize_once = false instead!";
    throw std::invalid_argument(msg.str());
  }
}

/*!*************************************************************************************************
 * \brief   Compute the local condensed matrix of every hyperedge and hand it to \c sink.
 *
 * For every hyperedge \c edge, \c sink(edge, hyNodes, loc_mat) is called, where \c hyNodes are the
 * indices of the hypernodes of the hyperedge and \c loc_mat points to the local condensed matrix
 * (cf. \c local_condensed_matrix for its layout).
 *
 * If \c factorize_once is \c true, the local matrix is only computed once for every combination of
 * face (hypernode) types, i.e., for every combination of boundary conditions, and reused for all
 * other hyperedges with the same face types. This is only correct if all these hyperedges have the
 * same local matrix, e.g. for uniform meshes with constant (in space) coefficients. It is wrong if,
 * for example, the diffusion tensor is space dependent or the hyperedges differ in size or shape.
 * Therefore, the local matrices of some hyperedges (the last one of every combination of face types
 * and 16 hyperedges evenly spread over all hyperedges) are recomputed and compared to the reused
 * ones, and an \c std::invalid_argument is thrown if they differ. This costs at most one local
 * matrix per combination of face types plus 16, which is negligible unless the mesh is small.
 * If the code is compiled without \c NDEBUG, this check is done for every hyperedge (which is as
 * expensive as \c factorize_once being false).
 *
 * The local matrices are computed in parallel if HyperHDG is compiled with OpenMP. Exceptions
 * thrown by the local solvers are forwarded to the caller.
 *
 * \tparam  hyNode_index_t  Typename of the hypernode index.
 * \tparam  parallel_sink   If true, \c sink is called in parallel (it must be safe to call it for
 *                          different hyperedges concurrently). If false, \c sink is called by one
 *                          thread in increasing order of the hyperedges, which allows to sum up
 *                          contributions of several hyperedges in a deterministic order.
 * \param   hyper_graph     The hypergraph.
 * \param   local_solver    The local solver.
 * \param   time            Parameter of \c trace_to_flux (time or shift).
 * \param   factorize_once  Compute the local matrix only once per combination of face types.
 * \param   sink            Callable that receives the local matrices.
 **************************************************************************************************/
template <typename hyNode_index_t,
          bool parallel_sink,
          typename HyperGraphT,
          typename LocalSolverT,
          typename dof_value_t,
          typename SinkT>
void for_each_local_matrix(HyperGraphT& hyper_graph,
                           const LocalSolverT& local_solver,
                           const dof_value_t time,
                           const bool factorize_once,
                           SinkT&& sink)
{
  constexpr unsigned int hyEdge_dim = HyperGraphT::hyEdge_dim();
  constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();
  constexpr std::size_t loc_size =
    (2 * hyEdge_dim * n_dofs_per_node) * (2 * hyEdge_dim * n_dofs_per_node);
  using face_types_t = std::array<unsigned int, 2 * hyEdge_dim>;
  using hyNodes_t = SmallVec<2 * hyEdge_dim, hyNode_index_t>;
  const std::size_t n_hyEdges = hyper_graph.n_hyEdges();

  // factorize_once: Compute the local matrix for every combination of face types (serially), and
  // check that it can be reused for a sample of hyperedges.
  std::map<face_types_t, std::vector<dof_value_t> > loc_mats_by_face_types;
  std::vector<const dof_value_t*> loc_mat_of_edge(factorize_once ? n_hyEdges : 0, nullptr);
  if (factorize_once)
  {
    std::vector<std::size_t> samples, first_edges;
    std::map<face_types_t, std::size_t> last_edge_of_face_types;
    for (std::size_t edge = 0; edge < n_hyEdges; ++edge)
    {
      auto hyper_edge = hyper_graph[edge];
      face_types_t face_types;
      for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
        face_types[face] = hyper_edge.node_descriptor[face];
      auto it = loc_mats_by_face_types.find(face_types);
      if (it == loc_mats_by_face_types.end())
      {
        it = loc_mats_by_face_types.emplace(face_types, std::vector<dof_value_t>(loc_size)).first;
        local_condensed_matrix(local_solver, hyper_edge, time, it->second.data());
        first_edges.push_back(edge);
      }
      loc_mat_of_edge[edge] = it->second.data();
      last_edge_of_face_types[face_types] = edge;
    }

    if constexpr (!hy_debug_mode)
    {
      for (const auto& face_types_and_edge : last_edge_of_face_types)
        samples.push_back(face_types_and_edge.second);
      const std::size_t n_spread = std::min<std::size_t>(n_hyEdges, 16);
      for (std::size_t k = 0; k < n_spread; ++k)
        samples.push_back(n_spread == 1 ? 0 : k * (n_hyEdges - 1) / (n_spread - 1));
      std::sort(samples.begin(), samples.end());
      samples.erase(std::unique(samples.begin(), samples.end()), samples.end());
      std::vector<dof_value_t> computed(loc_size);
      for (const std::size_t edge : samples)
      {
        // The local matrix of the first hyperedge of its face types is the reused one.
        if (std::binary_search(first_edges.begin(), first_edges.end(), edge))
          continue;
        auto hyper_edge = hyper_graph[edge];
        local_condensed_matrix(local_solver, hyper_edge, time, computed.data());
        check_reused_local_matrix(computed.data(), loc_mat_of_edge[edge], loc_size, edge);
      }
    }
  }

  // Local matrix of a hyperedge. It is computed into buffer (if necessary).
  auto local_matrix = [&](const std::size_t edge, auto& hyper_edge,
                          dof_value_t* buffer) -> const dof_value_t*
  {
    if (!factorize_once)
    {
      local_condensed_matrix(local_solver, hyper_edge, time, buffer);
      return buffer;
    }
    if constexpr (hy_debug_mode)
    {
      local_condensed_matrix(local_solver, hyper_edge, time, buffer);
      check_reused_local_matrix(buffer, loc_mat_of_edge[edge], loc_size, edge);
    }
    return loc_mat_of_edge[edge];
  };

  ParallelExceptions exceptions;
  if constexpr (parallel_sink)
  {
    HYPERHDG_OMP_PARALLEL_FOR
    for (long long edge = 0; edge < (long long)n_hyEdges; ++edge)
      exceptions.run(
        [&]()
        {
          auto hyper_edge = hyper_graph[edge];
          std::vector<dof_value_t> buffer(!factorize_once || hy_debug_mode ? loc_size : 0);
          const hyNodes_t hyNodes = hyper_edge.topology.get_hyNode_indices();
          sink((std::size_t)edge, hyNodes, local_matrix(edge, hyper_edge, buffer.data()));
        });
    exceptions.rethrow();
  }
  else
  {
    // Compute the local matrices of a chunk of hyperedges in parallel, and hand them to the sink
    // serially. The chunk size limits the memory of the buffer (to about 2^20 values, unless there
    // are many threads).
    const std::size_t chunk = std::min<std::size_t>(
      std::max<std::size_t>(n_hyEdges, 1),
      std::max<std::size_t>(8 * hyperhdg_n_threads(), (std::size_t(1) << 20) / loc_size));
    std::vector<dof_value_t> buffer(!factorize_once || hy_debug_mode ? chunk * loc_size : 0);
    std::vector<const dof_value_t*> loc_mats(chunk);
    std::vector<hyNodes_t> hyNodes(chunk);
    for (std::size_t begin = 0; begin < n_hyEdges; begin += chunk)
    {
      const std::size_t end = std::min(n_hyEdges, begin + chunk);
      HYPERHDG_OMP_PARALLEL_FOR
      for (long long edge = begin; edge < (long long)end; ++edge)
        exceptions.run(
          [&]()
          {
            const std::size_t k = edge - begin;
            auto hyper_edge = hyper_graph[edge];
            hyNodes[k] = hyper_edge.topology.get_hyNode_indices();
            loc_mats[k] =
              local_matrix(edge, hyper_edge, buffer.empty() ? nullptr : &buffer[k * loc_size]);
          });
      exceptions.rethrow();
      for (std::size_t edge = begin; edge < end; ++edge)
        sink(edge, hyNodes[edge - begin], loc_mats[edge - begin]);
    }
  }
}

/*!*************************************************************************************************
 * \brief   Assemble the condensed matrix in coordinate (triplet) format.
 *
 * This returns the same matrix as \c prototype_mat_generate (with the same ordering of entries),
 * but it obtains the local condensed matrix of each hyperedge at once (cf. \c for_each_local_matrix
 * for details and the meaning of \c factorize_once). Entries with the same row and column index
 * have to be summed up (as done by SciPy's \c csr_matrix). Every hyperedge writes its own block of
 * triplets, so the result does not depend on the number of OpenMP threads.
 *
 * \tparam  hyNode_index_t  Typename of the hypernode index.
 * \tparam  LargeVecT       Typename of the large vector.
 * \param   hyper_graph     The hypergraph.
 * \param   local_solver    The local solver.
 * \param   time            Parameter of \c trace_to_flux (time or shift).
 * \param   factorize_once  Compute the local matrix only once per combination of face types.
 * \retval  sparse_mat      Condensed matrix in coordinate format.
 **************************************************************************************************/
template <typename hyNode_index_t, typename LargeVecT, typename HyperGraphT, typename LocalSolverT>
sparse_mat<LargeVecT> assemble_triplets(HyperGraphT& hyper_graph,
                                        const LocalSolverT& local_solver,
                                        const typename LargeVecT::value_type time,
                                        const bool factorize_once)
{
  constexpr unsigned int hyEdge_dim = HyperGraphT::hyEdge_dim();
  constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();
  constexpr unsigned int n_loc = 2 * hyEdge_dim * n_dofs_per_node;
  check_dof_index_range<n_dofs_per_node>(hyper_graph);

  // The number of triplets may exceed the range of unsigned int (the type of the indices).
  const std::size_t n_triplets = (std::size_t)hyper_graph.n_hyEdges() * n_loc * n_loc;
  sparse_mat<LargeVecT> result_mat;
  result_mat.value_vec.resize(n_triplets);
  result_mat.row_vec.resize(n_triplets);
  result_mat.col_vec.resize(n_triplets);

  for_each_local_matrix<hyNode_index_t, true>(
    hyper_graph, local_solver, time, factorize_once,
    [&](const std::size_t edge, const auto& hyNodes, const auto* loc_mat)
    {
      std::array<std::array<unsigned int, n_dofs_per_node>, 2 * hyEdge_dim> dof_indices;
      for (unsigned int node = 0; node < hyNodes.size(); ++node)
        hyper_graph.hyNode_factory().get_dof_indices(hyNodes[node], dof_indices[node]);

      std::size_t index = edge * n_loc * n_loc;
      for (unsigned int node_j = 0; node_j < hyNodes.size(); ++node_j)
        for (unsigned int dof_j = 0; dof_j < n_dofs_per_node; ++dof_j)
          for (unsigned int node_i = 0; node_i < hyNodes.size(); ++node_i)
            for (unsigned int dof_i = 0; dof_i < n_dofs_per_node; ++dof_i, ++index)
            {
              result_mat.row_vec[index] = dof_indices[node_i][dof_i];
              result_mat.col_vec[index] = dof_indices[node_j][dof_j];
              result_mat.value_vec[index] = loc_mat[(node_j * n_dofs_per_node + dof_j) * n_loc +
                                                    node_i * n_dofs_per_node + dof_i];
            }
    });

  return result_mat;
}

}  // end of namespace MatGenerate
}  // end of namespace GlobalLoop
