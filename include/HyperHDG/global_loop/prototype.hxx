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
#define prototype_mat_vec_multiply(fun_name, has_fun_name)                                         \
  [&]() {                                                                                          \
    LargeVecT vec_Ax(x_vec.size(), 0.);                                                            \
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;                                              \
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> dofs_old, dofs_new;       \
                                                                                                   \
    _Pragma("acc parallel loop copyin(x_vec) copyout(vec_Ax) create(hyNodes, dofs_old, dofs_new)") \
    for (unsigned int index = 0; index < hyper_graph_.n_hyEdges(); ++index)                        \
    {                                                                                              \
      auto hyper_edge = hyper_graph_[index];                                                       \
      hyNodes = hyper_edge.topology.get_hyNode_indices();                                          \
      for (unsigned int node = 0; node < hyNodes.size(); ++node)                                   \
      {                                                                                            \
        hyper_graph_.hyNode_factory().get_dof_values(hyNodes[node], x_vec, dofs_old[node]);        \
        dofs_new[node].fill(0.);                                                                   \
      }                                                                                            \
                                                                                                   \
      if constexpr (has_fun_name<                                                                  \
                      LocalSolverT,                                                                \
                      std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&(       \
                        std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,     \
                        std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,     \
                        dof_value_t)>::value)                                                      \
        local_solver_.fun_name(dofs_old, dofs_new, time);                                          \
      else if constexpr (                                                                          \
        has_fun_name<LocalSolverT,                                                                 \
                     std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&(        \
                       std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,      \
                       std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>&,      \
                       decltype(hyper_edge)&, dof_value_t)>::value)                                \
        local_solver_.fun_name(dofs_old, dofs_new, hyper_edge, time);                              \
      else                                                                                         \
        hy_assert(false, "Function seems not to be implemented!");                                 \
                                                                                                   \
      for (unsigned int node = 0; node < hyNodes.size(); ++node)                                   \
        hyper_graph_.hyNode_factory().add_to_dof_values(hyNodes[node], vec_Ax, dofs_new[node]);    \
    };                                                                                             \
                                                                                                   \
    return vec_Ax;                                                                                 \
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
  [&]() {                                                                                          \
    typedef typename LocalSolverT::error_def::error_t error_t;                                     \
    error_t result = LocalSolverT::error_def::initial_error();                                     \
                                                                                                   \
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;                                              \
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> dofs;                     \
                                                                                                   \
    std::for_each(hyper_graph_.begin(), hyper_graph_.end(), [&](auto hyper_edge) {                 \
      hyNodes = hyper_edge.topology.get_hyNode_indices();                                          \
      for (unsigned int node = 0; node < hyNodes.size(); ++node)                                   \
        hyper_graph_.hyNode_factory().get_dof_values(hyNodes[node], x_vec, dofs[node]);            \
                                                                                                   \
      if constexpr (has_fun_name<LocalSolverT,                                                     \
                                 error_t(std::array<std::array<dof_value_t, n_dofs_per_node>,      \
                                                    2 * hyEdge_dim>&,                              \
                                         dof_value_t)>::value)                                     \
        result = LocalSolverT::error_def::sum_error(result, local_solver_.fun_name(dofs, time));   \
      else if constexpr (has_fun_name<LocalSolverT,                                                \
                                      error_t(std::array<std::array<dof_value_t, n_dofs_per_node>, \
                                                         2 * hyEdge_dim>&,                         \
                                              decltype(hyper_edge)&, dof_value_t)>::value)         \
        result = LocalSolverT::error_def::sum_error(                                               \
          result, local_solver_.fun_name(dofs, hyper_edge, time));                                 \
      else                                                                                         \
        hy_assert(false, "Function seems not to be ímplemented");                                  \
    });                                                                                            \
                                                                                                   \
    return LocalSolverT::error_def::postprocess_error(result);                                     \
  }()
