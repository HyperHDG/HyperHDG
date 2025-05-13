#pragma once  // Ensure that file is included only once in a single compilation.

/*!*************************************************************************************************
 * \brief   Macro that allows to use an implemented a matrix--vector multpilication.
 *
 * This macro implements a condensed matrix--vector multiplication commonly used in HyperHDG.
 * However, several requirements have to be meet for the macro to be used. That is:
 * - LargeVecT has to define the typename of the large vector,
 * - dof_value_t is the typename of a dof value,
 * - SmallVec of file dense_la has to be available,
 * - hyEdge_dim and n_dofs_per_node need to be constexpr variables,
 * - the algortihm library needs to be included,
 * - hyper_graph_ is an instance of HDGHyperGraph,
 * - has_fun_name is constructed via the macro HAS_MEMBER_FUNCTION,
 * - local_solver_ is a LocalSolverT,
 * - ... .
 **************************************************************************************************/
#define prototype_mat_vec_multiply(fun_name, has_fun_name)                                        \
  [&]()                                                                                           \
  {                                                                                               \
    LargeVecT vec_Ax(x_vec.size(), 0.);                                                           \
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;                                             \
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> dofs_old, dofs_new;      \
                                                                                                  \
    std::for_each(                                                                                \
      hyper_graph_.begin(), hyper_graph_.end(),                                                   \
      [&](auto hyper_edge)                                                                        \
      {                                                                                           \
        hyNodes = hyper_edge.topology.get_hyNode_indices();                                       \
        for (unsigned int node = 0; node < hyNodes.size(); ++node)                                \
        {                                                                                         \
          hyper_graph_.hyNode_factory().get_dof_values(hyNodes[node], x_vec, dofs_old[node]);     \
          dofs_new[node].fill(0.);                                                                \
        }                                                                                         \
                                                                                                  \
        if constexpr (has_fun_name<                                                               \
                        LocalSolverT,                                                             \
                        std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&(    \
                          std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,  \
                          std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,  \
                          param_time_t)>::value)                                                   \
          local_solver_.fun_name(dofs_old, dofs_new, time);                                       \
        else if constexpr (                                                                       \
          has_fun_name<LocalSolverT,                                                              \
                       std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&(     \
                         std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,   \
                         std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,   \
                         decltype(hyper_edge)&, param_time_t)>::value)                             \
          local_solver_.fun_name(dofs_old, dofs_new, hyper_edge, time);                           \
        else                                                                                      \
          hy_assert(false, "Function seems not to be implemented!");                              \
                                                                                                  \
        for (unsigned int node = 0; node < hyNodes.size(); ++node)                                \
          hyper_graph_.hyNode_factory().add_to_dof_values(hyNodes[node], vec_Ax, dofs_new[node]); \
      });                                                                                         \
                                                                                                  \
    return vec_Ax;                                                                                \
  }()

template <typename LargeVecT>
struct sparse_mat
{
  LargeVecT value_vec;
  std::vector<unsigned int> col_vec, row_vec;
  sparse_mat() {}
  sparse_mat(const unsigned int length) : value_vec(length), col_vec(length), row_vec(length) {}
  LargeVecT& get_values() { return value_vec; }
  std::vector<unsigned int>& get_cols() { return col_vec; }
  std::vector<unsigned int>& get_rows() { return row_vec; }
};

/*!*************************************************************************************************
 * \brief   Macro that allows to use an implemented a matrix generation.
 *
 * This macro implements the generation of a condensed matrix commonly used in HyperHDG.
 * However, several requirements have to be meet for the macro to be used. That is:
 * - LargeVecT has to define the typename of the large vector,
 * - dof_value_t is the typename of a dof value,
 * - SmallVec of file dense_la has to be available,
 * - hyEdge_dim and n_dofs_per_node need to be constexpr variables,
 * - the algortihm library needs to be included,
 * - hyper_graph_ is an instance of HDGHyperGraph,
 * - has_fun_name is constructed via the macro HAS_MEMBER_FUNCTION,
 * - local_solver_ is a LocalSolverT,
 * - ... .
 **************************************************************************************************/
#define prototype_mat_generate(fun_name, has_fun_name)                                        \
  [&]()                                                                                       \
  {                                                                                           \
    sparse_mat<LargeVecT> result_mat(hyper_graph_.n_hyEdges() * 4 * hyEdge_dim * hyEdge_dim * \
                                     n_dofs_per_node * n_dofs_per_node);                      \
                                                                                              \
    auto value_it = result_mat.value_vec.begin();                                             \
    auto col_it = result_mat.col_vec.begin(), row_it = result_mat.row_vec.begin();            \
                                                                                              \
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;                                         \
    std::array<std::array<unsigned int, n_dofs_per_node>, 2 * hyEdge_dim> dof_indices;        \
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> dofs_old, dofs_new;  \
                                                                                              \
    std::for_each(                                                                            \
      hyper_graph_.begin(), hyper_graph_.end(),                                               \
      [&](auto hyper_edge)                                                                    \
      {                                                                                       \
        hyNodes = hyper_edge.topology.get_hyNode_indices();                                   \
        for (unsigned int node = 0; node < hyNodes.size(); ++node)                            \
          hyper_graph_.hyNode_factory().get_dof_indices(hyNodes[node], dof_indices[node]);    \
                                                                                              \
        for (unsigned int node_j = 0; node_j < hyNodes.size(); ++node_j)                      \
          for (unsigned int dof_j = 0; dof_j < n_dofs_per_node; ++dof_j)                      \
          {                                                                                   \
            for (unsigned int node = 0; node < hyNodes.size(); ++node)                        \
            {                                                                                 \
              dofs_old[node].fill(0.);                                                        \
              dofs_new[node].fill(0.);                                                        \
            }                                                                                 \
            dofs_old[node_j][dof_j] = 1.;                                                     \
            if constexpr (                                                                    \
              has_fun_name<                                                                   \
                LocalSolverT,                                                                 \
                std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&(        \
                  std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,      \
                  std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,      \
                  param_time_t)>::value)                                                       \
              local_solver_.fun_name(dofs_old, dofs_new, time);                               \
            else if constexpr (                                                               \
              has_fun_name<                                                                   \
                LocalSolverT,                                                                 \
                std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&(        \
                  std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,      \
                  std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,      \
                  decltype(hyper_edge)&, param_time_t)>::value)                                \
              local_solver_.fun_name(dofs_old, dofs_new, hyper_edge, time);                   \
            else                                                                              \
              hy_assert(false, "Function seems not to be implemented!");                      \
                                                                                              \
            for (unsigned int node_i = 0; node_i < hyNodes.size(); ++node_i)                  \
              for (unsigned int dof_i = 0; dof_i < n_dofs_per_node; ++dof_i)                  \
              {                                                                               \
                *(row_it++) = dof_indices[node_i][dof_i];                                     \
                *(col_it++) = dof_indices[node_j][dof_j];                                     \
                *(value_it++) = dofs_new[node_i][dof_i];                                      \
              }                                                                               \
          }                                                                                   \
      });                                                                                     \
                                                                                              \
    return result_mat;                                                                        \
  }()

/*!*************************************************************************************************
 * \brief   Macro that allows to use an implemented error evaluation.
 *
 * This macro implements a condensed matrix--vector multiplication commonly used in HyperHDG.
 * However, several requirements have to be meet for the macro to be used. That is:
 * - LargeVecT has to define the typename of the large vector,
 * - dof_value_t is the typename of a dof value,
 * - SmallVec of file dense_la has to be available,
 * - hyEdge_dim and n_dofs_per_node need to be constexpr variables,
 * - the algortihm library needs to be included,
 * - hyper_graph_ is an instance of HDGHyperGraph,
 * - has_fun_name is constructed via the macro HAS_MEMBER_FUNCTION,
 * - local_solver_ is a LocalSolverT,
 * - ... .
 **************************************************************************************************/
#define prototype_errors(fun_name, has_fun_name)                                                   \
  [&]()                                                                                            \
  {                                                                                                \
    typedef typename LocalSolverT::error_def::error_t error_t;                                     \
    error_t result = LocalSolverT::error_def::initial_error();                                     \
                                                                                                   \
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;                                              \
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> dofs;                     \
                                                                                                   \
    std::for_each(                                                                                 \
      hyper_graph_.begin(), hyper_graph_.end(),                                                    \
      [&](auto hyper_edge)                                                                         \
      {                                                                                            \
        hyNodes = hyper_edge.topology.get_hyNode_indices();                                        \
        for (unsigned int node = 0; node < hyNodes.size(); ++node)                                 \
          hyper_graph_.hyNode_factory().get_dof_values(hyNodes[node], x_vec, dofs[node]);          \
                                                                                                   \
        if constexpr (has_fun_name<LocalSolverT,                                                   \
                                   error_t(std::array<std::array<dof_value_t, n_dofs_per_node>,    \
                                                      2 * hyEdge_dim>&,                            \
                                           param_time_t)>::value)                                   \
          result = LocalSolverT::error_def::sum_error(result, local_solver_.fun_name(dofs, time)); \
        else if constexpr (has_fun_name<LocalSolverT,                                              \
                                        error_t(                                                   \
                                          std::array<std::array<dof_value_t, n_dofs_per_node>,     \
                                                     2 * hyEdge_dim>&,                             \
                                          decltype(hyper_edge)&, param_time_t)>::value)             \
          result = LocalSolverT::error_def::sum_error(                                             \
            result, local_solver_.fun_name(dofs, hyper_edge, time));                               \
        else                                                                                       \
          hy_assert(false, "Function seems not to be ímplemented");                                \
      });                                                                                          \
                                                                                                   \
    return LocalSolverT::error_def::postprocess_error(result);                                     \
  }()
  
/*!*************************************************************************************************
 * \brief   Macro that allows to use an implemented mean evaluation.
 *
 * This macro implements a condensed matrix--vector multiplication commonly used in HyperHDG.
 * However, several requirements have to be meet for the macro to be used. That is:
 * - LargeVecT has to define the typename of the large vector,
 * - dof_value_t is the typename of a dof value,
 * - SmallVec of file dense_la has to be available,
 * - hyEdge_dim and n_dofs_per_node need to be constexpr variables,
 * - the algortihm library needs to be included,
 * - hyper_graph_ is an instance of HDGHyperGraph,
 * - has_fun_name is constructed via the macro HAS_MEMBER_FUNCTION,
 * - local_solver_ is a LocalSolverT,
 * - ... .
 **************************************************************************************************/
#define prototype_mean(fun_name, has_fun_name)                                                   \
  [&]()                                                                                            \
  {                                                                                                \
    std::vector<dof_value_t> result(hyEdge_dim + 1);                                               \
                                                                                                   \
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;                                              \
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> dofs;                     \
                                                                                                   \
    std::for_each(                                                                                 \
      hyper_graph_.begin(), hyper_graph_.end(),                                                    \
      [&](auto hyper_edge)                                                                         \
      {                                                                                            \
        hyNodes = hyper_edge.topology.get_hyNode_indices();                                        \
        for (unsigned int node = 0; node < hyNodes.size(); ++node)                                 \
          hyper_graph_.hyNode_factory().get_dof_values(hyNodes[node], x_vec, dofs[node]);          \
                                                                                                   \
        if constexpr (has_fun_name<LocalSolverT,                                                   \
                                   std::array<dof_value_t, hyEdge_dim + 1>(                        \
                                        std::array<std::array<dof_value_t, n_dofs_per_node>,       \
                                                      2 * hyEdge_dim>&,                            \
                                                      param_time_t)>::value) {                     \
          std::array<dof_value_t, hyEdge_dim + 1> loc_contr = local_solver_.fun_name(dofs, time);  \
          for(unsigned int i = 0; i <= hyEdge_dim; i++)                                            \
            result[i] += loc_contr[i];                                                             \
          }                                                                                        \
        else if constexpr (has_fun_name<LocalSolverT,                                              \
                                        std::array<dof_value_t, hyEdge_dim + 1>(                   \
                                          std::array<std::array<dof_value_t, n_dofs_per_node>,     \
                                                     2 * hyEdge_dim>&,                             \
                                          decltype(hyper_edge)&, param_time_t)>::value) {          \
          std::array<dof_value_t, hyEdge_dim + 1> loc_contr                                        \
            = local_solver_.fun_name(dofs, hyper_edge, time);                                      \
          for(unsigned int i = 0; i <= hyEdge_dim; i++)                                            \
            result[i] += loc_contr[i];                                                             \
          }                                                                                        \
        else                                                                                       \
          hy_assert(false, "Function seems not to be ímplemented");                                \
      });                                                                                          \
                                                                                                   \
    return result;                                                                                 \
  }()
