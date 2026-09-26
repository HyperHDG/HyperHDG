#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/global_loop/prototype.hxx>
#include <HyperHDG/hy_assert.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <vector>

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
#else
#define HYPERHDG_OMP_PARALLEL_FOR
#endif

/*!*************************************************************************************************
 * \brief   Macro that allows to use an implemented matrix generation based on local matrices.
 *
 * This macro implements the generation of a condensed matrix commonly used in HyperHDG. It returns
 * the same matrix as \c prototype_mat_generate (with the same ordering of entries), but it obtains
 * the local condensed matrix of each hyperedge at once:
 * - If the local solver provides the function \c mat_name(hyper_edge, time), which returns the
 *   local condensed matrix as a \c SmallSquareMat, this function is used. Local solvers can
 *   implement it by factorizing their local matrix once and reusing the factorization for all
 *   unit vectors.
 * - Otherwise, the local condensed matrix is obtained by applying \c fun_name to all unit vectors
 *   (which is what \c prototype_mat_generate does).
 *
 * If \c factorize_once is \c true, the local condensed matrix is only computed once for every
 * combination of face (hypernode) types, i.e., for every combination of boundary conditions, and
 * reused for all other hyperedges with the same face types. This is only correct if all these
 * hyperedges have the same local matrix, e.g. for uniform meshes with constant (in space)
 * coefficients. It is wrong if, for example, the diffusion tensor is space dependent or the
 * hyperedges differ in size or shape. If the code is compiled without \c NDEBUG, this assumption is
 * checked for every hyperedge (which is as expensive as \c factorize_once being \c false).
 *
 * Several requirements have to be meet for the macro to be used. That is:
 * - LargeVecT has to define the typename of the large vector,
 * - dof_value_t is the typename of a dof value,
 * - SmallVec of file dense_la has to be available,
 * - hyEdge_dim and n_dofs_per_node need to be constexpr variables,
 * - hyper_graph_ is an instance of HDGHyperGraph,
 * - has_fun_name and has_mat_name are constructed via the macro HAS_MEMBER_FUNCTION,
 * - local_solver_ is a LocalSolverT,
 * - time is the time at which the matrix is generated.
 **************************************************************************************************/
#define prototype_mat_generate_local(fun_name, has_fun_name, mat_name, has_mat_name,               \
                                     factorize_once)                                               \
  [&]()                                                                                            \
  {                                                                                                \
    constexpr unsigned int n_loc = 2 * hyEdge_dim * n_dofs_per_node;                               \
    typedef std::vector<dof_value_t>                                                               \
      loc_mat_t; /* column-major: (row, col) -> col * n_loc + row */                               \
                                                                                                   \
    const auto n_hyEdges = hyper_graph_.n_hyEdges();                                               \
    sparse_mat<LargeVecT> result_mat(n_hyEdges* n_loc* n_loc);                                     \
                                                                                                   \
    auto compute_loc_mat = [&](auto& hyper_edge, auto& solver, loc_mat_t& mat)                     \
    {                                                                                              \
      if constexpr (has_mat_name<LocalSolverT, SmallSquareMat<n_loc, dof_value_t>(                 \
                                                 decltype(hyper_edge)&, dof_value_t)>::value)      \
      {                                                                                            \
        SmallSquareMat<n_loc, dof_value_t> helper = solver.mat_name(hyper_edge, time);             \
        for (unsigned int col = 0; col < n_loc; ++col)                                             \
          for (unsigned int row = 0; row < n_loc; ++row)                                           \
            mat[col * n_loc + row] = helper(row, col);                                             \
      }                                                                                            \
      else                                                                                         \
      {                                                                                            \
        std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> dofs_old, dofs_new;   \
        for (unsigned int node_j = 0; node_j < 2 * hyEdge_dim; ++node_j)                           \
          for (unsigned int dof_j = 0; dof_j < n_dofs_per_node; ++dof_j)                           \
          {                                                                                        \
            for (unsigned int node = 0; node < 2 * hyEdge_dim; ++node)                             \
            {                                                                                      \
              dofs_old[node].fill(0.);                                                             \
              dofs_new[node].fill(0.);                                                             \
            }                                                                                      \
            dofs_old[node_j][dof_j] = 1.;                                                          \
            if constexpr (                                                                         \
              has_fun_name<                                                                        \
                LocalSolverT,                                                                      \
                std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&(             \
                  std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,           \
                  std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,           \
                  dof_value_t)>::value)                                                            \
              solver.fun_name(dofs_old, dofs_new, time);                                           \
            else if constexpr (                                                                    \
              has_fun_name<                                                                        \
                LocalSolverT,                                                                      \
                std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&(             \
                  std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,           \
                  std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,           \
                  decltype(hyper_edge)&, dof_value_t)>::value)                                     \
              solver.fun_name(dofs_old, dofs_new, hyper_edge, time);                               \
            else                                                                                   \
              hy_assert(false, "Function seems not to be implemented!");                           \
                                                                                                   \
            const unsigned int col = node_j * n_dofs_per_node + dof_j;                             \
            for (unsigned int node_i = 0; node_i < 2 * hyEdge_dim; ++node_i)                       \
              for (unsigned int dof_i = 0; dof_i < n_dofs_per_node; ++dof_i)                       \
                mat[col * n_loc + node_i * n_dofs_per_node + dof_i] = dofs_new[node_i][dof_i];     \
          }                                                                                        \
      }                                                                                            \
    };                                                                                             \
                                                                                                   \
    /* factorize_once: compute the local matrix for every combination of face types (serially). */ \
    std::map<std::array<unsigned int, 2 * hyEdge_dim>, loc_mat_t> loc_mats_by_face_types;          \
    std::vector<const loc_mat_t*> loc_mat_of_edge(factorize_once ? n_hyEdges : 0, nullptr);        \
    if (factorize_once)                                                                            \
      for (std::size_t edge = 0; edge < n_hyEdges; ++edge)                                         \
      {                                                                                            \
        auto hyper_edge = hyper_graph_[edge];                                                      \
        std::array<unsigned int, 2 * hyEdge_dim> face_types;                                       \
        for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)                                 \
          face_types[face] = hyper_edge.node_descriptor[face];                                     \
        auto it = loc_mats_by_face_types.find(face_types);                                         \
        if (it == loc_mats_by_face_types.end())                                                    \
        {                                                                                          \
          it = loc_mats_by_face_types.emplace(face_types, loc_mat_t(n_loc * n_loc)).first;         \
          compute_loc_mat(hyper_edge, local_solver_, it->second);                                  \
        }                                                                                          \
        loc_mat_of_edge[edge] = &it->second;                                                       \
      }                                                                                            \
                                                                                                   \
    /* Every hyperedge writes to its own part of result_mat, which allows for parallelization. */  \
    HYPERHDG_OMP_PARALLEL_FOR                                                                      \
    for (long long edge = 0; edge < (long long)n_hyEdges; ++edge)                                  \
    {                                                                                              \
      auto hyper_edge = hyper_graph_[edge];                                                        \
      SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes = hyper_edge.topology.get_hyNode_indices(); \
      std::array<std::array<unsigned int, n_dofs_per_node>, 2 * hyEdge_dim> dof_indices;           \
      for (unsigned int node = 0; node < hyNodes.size(); ++node)                                   \
        hyper_graph_.hyNode_factory().get_dof_indices(hyNodes[node], dof_indices[node]);           \
                                                                                                   \
      loc_mat_t loc_mat;                                                                           \
      const loc_mat_t* mat_ptr = factorize_once ? loc_mat_of_edge[edge] : &loc_mat;                \
      if (!factorize_once || hy_debug_mode)                                                        \
      {                                                                                            \
        loc_mat.resize(n_loc* n_loc);                                                              \
        compute_loc_mat(hyper_edge, local_solver_, loc_mat);                                       \
      }                                                                                            \
      if (factorize_once && hy_debug_mode)                                                         \
        for (unsigned int k = 0; k < n_loc * n_loc; ++k)                                           \
          hy_assert(                                                                               \
            std::abs(loc_mat[k] - (*mat_ptr)[k]) <= 1e-10 * (1. + std::abs((*mat_ptr)[k])),        \
            "Option factorize_once requires all hyperedges with the same face types "              \
              << "to have the same local matrix (e.g. uniform meshes and coefficients "            \
              << "that are constant in space). Use factorize_once = false instead!");              \
                                                                                                   \
      std::size_t index = (std::size_t)edge * n_loc * n_loc;                                       \
      for (unsigned int node_j = 0; node_j < hyNodes.size(); ++node_j)                             \
        for (unsigned int dof_j = 0; dof_j < n_dofs_per_node; ++dof_j)                             \
          for (unsigned int node_i = 0; node_i < hyNodes.size(); ++node_i)                         \
            for (unsigned int dof_i = 0; dof_i < n_dofs_per_node; ++dof_i, ++index)                \
            {                                                                                      \
              result_mat.row_vec[index] = dof_indices[node_i][dof_i];                              \
              result_mat.col_vec[index] = dof_indices[node_j][dof_j];                              \
              result_mat.value_vec[index] =                                                        \
                (*mat_ptr)[(node_j * n_dofs_per_node + dof_j) * n_loc + node_i * n_dofs_per_node + \
                           dof_i];                                                                 \
            }                                                                                      \
    }                                                                                              \
                                                                                                   \
    return result_mat;                                                                             \
  }()
