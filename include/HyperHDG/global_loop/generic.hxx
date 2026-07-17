#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/global_loop/prototype.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/plot.hxx>
#include <algorithm>
#include <array>
#include <vector>

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Consolidated global loop: one class for every problem type.
 *
 * Lifecycle follows the FESTUNG generic problem framework: configure (runtime setup, zero node
 * types) -> initialize (initial state) -> per step { residual, jacobian, external solve,
 * postprocess } -> output (errors / plot). The elliptic/parabolic/hyperbolic loops are subsets of
 * this protocol; they remain in place until their users migrate.
 *
 * Every operator entry is templated on a generic argument type ArgT that is threaded through to
 * the local solver untouched: a plain scalar carries the time, a Gauss::StageTime carries
 * (time, stage), applications may define their own. The loop never inspects the argument.
 *
 * Zero handling: configure takes a runtime list of node types (e.g. the ternary Cubic encoding);
 * every node whose descriptor matches one of them exactly has its dofs collected into the zero
 * index sets. jacobian/residual zero those entries of the output span, jacobian_mat zeroes the
 * matching rows AND columns and puts a single 1 on each owned diagonal. Per-dof control (partial
 * constraints, static-only bits) stays in the local solver (is_dirichlet et al.); both mechanisms
 * compose.
 *
 * \tparam  TopologyT       Class type containing topological information.
 * \tparam  GeometryT       Class type containing geometrical information.
 * \tparam  NodeDescriptorT Class type containing the information of nodes of hyperedges.
 * \tparam  LocalSolverT    Class type of the local solver.
 * \tparam  LargeVecT       Class type of large, global vector.
 * \tparam  dof_index_t     Index type of hyperedges. Default is \c unsigned \c int.
 **************************************************************************************************/
template <class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename LargeVecT = std::vector<double>,
          typename dof_index_t = unsigned int>
class Generic
{
 public:
  /*!***********************************************************************************************
   * \brief   Dimension of a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
  /*!***********************************************************************************************
   * \brief   Global degrees of freedom per hypernode.
   ************************************************************************************************/
  static constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();
  /*!***********************************************************************************************
   * \brief   Floating type is determined by floating type of large vector's entries.
   ************************************************************************************************/
  using dof_value_t = typename LargeVecT::value_type;

 private:
  /*!***********************************************************************************************
   * \brief   Instantiation of a hypergraph.
   ************************************************************************************************/
  HDGHyperGraph<LocalSolverT::n_glob_dofs_per_node(),
                TopologyT,
                GeometryT,
                NodeDescriptorT,
                typename LocalSolverT::data_type>
    hyper_graph_;
  /*!***********************************************************************************************
   * \brief   Instantiation of a local solver. Non-const so configure can reach it.
   ************************************************************************************************/
  LocalSolverT local_solver_;
  /*!***********************************************************************************************
   * \brief   Struct encoding the options for plotting.
   ************************************************************************************************/
  PlotOptions plot_options;
  /*!***********************************************************************************************
   * \brief   Zero dof indices in local numbering (matches the n_local_dofs-sized spans).
   ************************************************************************************************/
  std::vector<dof_index_t> zero_indices_local_;
  /*!***********************************************************************************************
   * \brief   Zero dof indices in global numbering (matches the COO rows/cols of jacobian_mat).
   ************************************************************************************************/
  std::vector<dof_index_t> zero_indices_global_;
  /*!***********************************************************************************************
   * \brief   Owned subset of \c zero_indices_global_: exactly one diagonal 1 per zero dof
   *          globally, since PETSc-style COO assembly sums duplicate triplets.
   ************************************************************************************************/
  std::vector<dof_index_t> zero_indices_diag_;

  /*!***********************************************************************************************
   * \brief   Collect the dof indices of every node whose descriptor is in the given type list.
   *
   * Exact-match semantics on \c hyper_edge.node_descriptor[node]; assumes both edges incident to
   * a node see the same type (also across ranks: every rank touching the node zeroes its local
   * entries, so the cross-rank add stays consistent). Ownership of a dof is "local index below
   * n_owned_dofs" -- hypernodes are numbered owned-first (see distribute_domain.hxx).
   ************************************************************************************************/
  void build_zero_indices(const std::vector<unsigned int>& zero_node_types)
  {
    zero_indices_local_.clear();
    zero_indices_global_.clear();
    zero_indices_diag_.clear();

    std::vector<unsigned int> types = zero_node_types;
    std::sort(types.begin(), types.end());

    const dof_index_t n_owned = hyper_graph_.n_owned_dofs();
    std::array<dof_index_t, n_dofs_per_node> idx_loc, idx_glob;

    const dof_index_t n_e = hyper_graph_.n_hyEdges();
    for (dof_index_t iedge = 0; iedge < n_e; ++iedge)
    {
      auto hyper_edge = hyper_graph_[iedge];
      const auto hyNodes = hyper_edge.topology.get_hyNode_indices();
      for (unsigned int node = 0; node < hyNodes.size(); ++node)
      {
        if (!std::binary_search(types.begin(), types.end(),
                                static_cast<unsigned int>(hyper_edge.node_descriptor[node])))
          continue;
        hyper_graph_.hyNode_factory().get_dof_indices(hyNodes[node], idx_loc);
        hyper_graph_.hyNode_factory().get_global_dof_indices(hyNodes[node], idx_glob);
        for (unsigned int d = 0; d < n_dofs_per_node; ++d)
        {
          zero_indices_local_.push_back(idx_loc[d]);
          zero_indices_global_.push_back(idx_glob[d]);
          if (idx_loc[d] < n_owned)
            zero_indices_diag_.push_back(idx_glob[d]);
        }
      }
    }

    for (auto* vec : {&zero_indices_local_, &zero_indices_global_, &zero_indices_diag_})
    {
      std::sort(vec->begin(), vec->end());
      vec->erase(std::unique(vec->begin(), vec->end()), vec->end());
    }
  }
  /*!***********************************************************************************************
   * \brief   Zero the configured entries of an output span.
   ************************************************************************************************/
  template <typename SpanT>
  void zero_span(SpanT& vec) const
  {
    for (const dof_index_t idx : zero_indices_local_)
    {
      hy_assert(idx < vec.size(), "zero index " << idx << " exceeds vector size " << vec.size());
      vec[idx] = 0.;
    }
  }
  /*!***********************************************************************************************
   * \brief   Zero the configured rows and columns of an assembled COO matrix, diagonal 1.
   *
   * Blanking by row OR column membership also wipes contributions of neighboring edges (which the
   * per-edge is_dirichlet shortcut of the assembly cannot see) and any per-edge diagonal entries,
   * so the appended owned diagonal 1 is exact -- no double counting under COO duplicate summing.
   ************************************************************************************************/
  template <typename MatVecT>
  void zero_mat(sparse_mat<MatVecT>& mat) const
  {
    if (zero_indices_global_.empty())
      return;
    const auto is_zero = [&](const unsigned int idx)
    { return std::binary_search(zero_indices_global_.begin(), zero_indices_global_.end(), idx); };
    for (std::size_t i = 0; i < mat.value_vec.size(); ++i)
      if (is_zero(mat.row_vec[i]) || is_zero(mat.col_vec[i]))
        mat.value_vec[i] = 0.;
    for (const dof_index_t idx : zero_indices_diag_)
    {
      mat.row_vec.push_back(idx);
      mat.col_vec.push_back(idx);
      mat.value_vec.push_back(1.);
    }
  }

 public:
  /*!***********************************************************************************************
   * \brief   Abstract problem constructor.
   *
   * \param   construct_topo    Information to construct a topology.
   * \param   construct_geom    Information to construct a geometry.
   * \param   construct_loc_sol Information to construct a local solver.
   ************************************************************************************************/
  Generic(const typename TopologyT::constructor_value_type& construct_topo,
          const typename GeometryT::constructor_value_type& construct_geom,
          const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo, construct_geom), local_solver_(construct_loc_sol)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Abstract problem constructor.
   *
   * \param   construct_topo    Information to construct a topology.
   * \param   construct_loc_sol Information to construct a local solver.
   ************************************************************************************************/
  Generic(const typename TopologyT::constructor_value_type& construct_topo,
          const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo), local_solver_(construct_loc_sol)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Abstract problem constructor.
   *
   * \param   construct_topo    Information to construct a topology.
   ************************************************************************************************/
  Generic(const typename TopologyT::constructor_value_type& construct_topo)
  : hyper_graph_(construct_topo)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }

  // -----------------------------------------------------------------------------------------------
  // Lifecycle: configure.
  // -----------------------------------------------------------------------------------------------
  /*!***********************************************************************************************
   * \brief   Runtime problem configuration.
   *
   * Builds the zero index sets from the given node-type list and forwards any further arguments
   * verbatim to \c LocalSolverT::configure. Passing local config to a solver without a matching
   * configure is a compile error (silently dropping runtime configuration would be a footgun).
   *
   * Must be re-run after \c set_refinement (which clears the zero sets: the dof layout changes).
   ************************************************************************************************/
  template <typename... LocalCfgTs>
  void configure(const std::vector<unsigned int>& zero_node_types, LocalCfgTs&&... local_cfg)
  {
    build_zero_indices(zero_node_types);
    if constexpr (sizeof...(LocalCfgTs) > 0)
    {
      static_assert(requires { local_solver_.configure(std::forward<LocalCfgTs>(local_cfg)...); },
                    "LocalSolverT has no configure() matching the given arguments!");
      local_solver_.configure(std::forward<LocalCfgTs>(local_cfg)...);
    }
  }
  /*!***********************************************************************************************
   * \brief   Zero dof indices in local numbering (introspection / cross-validation).
   ************************************************************************************************/
  const std::vector<dof_index_t>& zero_indices_local() const { return zero_indices_local_; }
  /*!***********************************************************************************************
   * \brief   Zero dof indices in global numbering (introspection / cross-validation).
   ************************************************************************************************/
  const std::vector<dof_index_t>& zero_indices_global() const { return zero_indices_global_; }

  // -----------------------------------------------------------------------------------------------
  // Lifecycle: initialize. Optional capability, detected via constexpr: a solver without
  // make_initial gets a no-op (the user opted into this semantics); a solver whose make_initial
  // exists but does not accept the given ArgT is a compile error, not a silent skip.
  // -----------------------------------------------------------------------------------------------
  /*!***********************************************************************************************
   * \brief   Write the initial trace into the caller's span (zero-filled first) and set up the
   *          local solvers' initial state.
   ************************************************************************************************/
  template <typename SpanT, typename ArgT = dof_value_t, typename hyNode_index_t = dof_index_t>
  void initialize(SpanT& x_vec, const ArgT arg = 0.)
  {
    using span_value_t = typename std::decay_t<SpanT>::value_type;
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;
    std::array<std::array<span_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

    std::fill(x_vec.begin(), x_vec.end(), 0.);
    std::for_each(
      hyper_graph_.begin(), hyper_graph_.end(),
      [&](auto hyper_edge)
      {
        hyNodes = hyper_edge.topology.get_hyNode_indices();
        for (unsigned int node = 0; node < hyNodes.size(); ++node)
          hyEdge_dofs[node].fill(0.);

        if constexpr (requires { local_solver_.make_initial(hyEdge_dofs, prvalue_of(arg)); })
          local_solver_.make_initial(hyEdge_dofs, arg);
        else if constexpr (requires {
                             local_solver_.make_initial(hyEdge_dofs, hyper_edge, prvalue_of(arg));
                           })
          local_solver_.make_initial(hyEdge_dofs, hyper_edge, arg);
        else if constexpr (requires { local_solver_.make_initial(hyEdge_dofs, 0.); } ||
                           requires { local_solver_.make_initial(hyEdge_dofs, hyper_edge, 0.); })
          static_assert(always_false_v<LocalSolverT, decltype(hyper_edge)>,
                        "LocalSolverT has make_initial, but no overload for this ArgT!");
        else
          return;  // no make_initial: initialize is a no-op

        for (unsigned int node = 0; node < hyNodes.size(); ++node)
          hyper_graph_.hyNode_factory().set_dof_values(hyNodes[node], x_vec, hyEdge_dofs[node]);
      });
  }
  /*!***********************************************************************************************
   * \brief   Initialize the state by solving the static problem for the given boundary trace
   *          (reads and rewrites the caller's span). No-op without local support.
   ************************************************************************************************/
  template <typename SpanT, typename ArgT = dof_value_t, typename hyNode_index_t = dof_index_t>
  void initialize_from_static(SpanT& x_vec, const ArgT arg = 0.)
  {
    using span_value_t = typename std::decay_t<SpanT>::value_type;
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;
    std::array<std::array<span_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

    std::for_each(
      hyper_graph_.begin(), hyper_graph_.end(),
      [&](auto hyper_edge)
      {
        hyNodes = hyper_edge.topology.get_hyNode_indices();
        for (unsigned int node = 0; node < hyNodes.size(); ++node)
          hyper_graph_.hyNode_factory().get_dof_values(hyNodes[node], x_vec, hyEdge_dofs[node]);

        if constexpr (requires {
                        local_solver_.make_initial_from_static(hyEdge_dofs, prvalue_of(arg));
                      })
          local_solver_.make_initial_from_static(hyEdge_dofs, arg);
        else if constexpr (requires {
                             local_solver_.make_initial_from_static(hyEdge_dofs, hyper_edge,
                                                                    prvalue_of(arg));
                           })
          local_solver_.make_initial_from_static(hyEdge_dofs, hyper_edge, arg);
        else if constexpr (requires { local_solver_.make_initial_from_static(hyEdge_dofs, 0.); } ||
                           requires {
                             local_solver_.make_initial_from_static(hyEdge_dofs, hyper_edge, 0.);
                           })
          static_assert(always_false_v<LocalSolverT, decltype(hyper_edge)>,
                        "LocalSolverT has make_initial_from_static, but no overload for this "
                        "ArgT!");
        else
          return;  // no make_initial_from_static: initialize_from_static is a no-op

        for (unsigned int node = 0; node < hyNodes.size(); ++node)
          hyper_graph_.hyNode_factory().set_dof_values(hyNodes[node], x_vec, hyEdge_dofs[node]);
      });
  }
  /*!***********************************************************************************************
   * \brief   Post-process step: hand the solved trace to the local solvers (FESTUNG's
   *          "Post-process Step"; the former set_data).
   *
   * For Gauss stage solvers the stage index rides in the argument (Gauss::StageTime); the
   * sentinel stage == -1 means "all stages are in" and turns this call into the step
   * recombination (the former finalize_step -- not a separate API here). Solvers without
   * set_data get a no-op.
   ************************************************************************************************/
  template <typename SpanT, typename ArgT = dof_value_t, typename hyNode_index_t = dof_index_t>
  void postprocess(const SpanT& x_vec, const ArgT arg = 0.)
  {
    using span_value_t = typename std::decay_t<SpanT>::value_type;
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyNodes;
    std::array<std::array<span_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

    std::for_each(
      hyper_graph_.begin(), hyper_graph_.end(),
      [&](auto hyper_edge)
      {
        hyNodes = hyper_edge.topology.get_hyNode_indices();
        for (unsigned int node = 0; node < hyNodes.size(); ++node)
          hyper_graph_.hyNode_factory().get_dof_values(hyNodes[node], x_vec, hyEdge_dofs[node]);

        if constexpr (requires { local_solver_.set_data(hyEdge_dofs, prvalue_of(arg)); })
          local_solver_.set_data(hyEdge_dofs, arg);
        else if constexpr (requires {
                             local_solver_.set_data(hyEdge_dofs, hyper_edge, prvalue_of(arg));
                           })
          local_solver_.set_data(hyEdge_dofs, hyper_edge, arg);
        else if constexpr (requires { local_solver_.set_data(hyEdge_dofs, 0.); } ||
                           requires { local_solver_.set_data(hyEdge_dofs, hyper_edge, 0.); })
          static_assert(always_false_v<LocalSolverT, decltype(hyper_edge)>,
                        "LocalSolverT has set_data, but no overload for this ArgT!");
        // else: no set_data -- postprocess is a no-op (steady problems)
      });
  }

  // -----------------------------------------------------------------------------------------------
  // Core operators. No allocation: outputs are caller-provided spans / sparse_mat buffers.
  // -----------------------------------------------------------------------------------------------
  /*!***********************************************************************************************
   * \brief   Apply the condensed system operator: vec_Ax = A(arg) x_vec.
   *
   * The homogeneous part of the problem (the matrix the solver inverts); dispatches to the local
   * solver's trace_to_flux. Zeroes the configured zero entries of the output.
   ************************************************************************************************/
  template <typename XSpanT, typename YSpanT, typename ArgT = dof_value_t,
            typename hyNode_index_t = dof_index_t>
  void jacobian(const XSpanT& x_vec, YSpanT& vec_Ax, const ArgT arg = 0.)
  {
    hy_assert(x_vec.size() == vec_Ax.size(), "x_vec and vec_Ax need to be of same size");
    const ArgT& time = arg;
    prototype_mat_vec_multiply_span(trace_to_flux);
    zero_span(vec_Ax);
  }
  /*!***********************************************************************************************
   * \brief   Residual of the condensed system: operator action plus data terms (loads, boundary
   *          values, old-step history).
   *
   * Dispatches to the local solver's residual_flux. Zeroes the configured zero entries of the
   * output.
   ************************************************************************************************/
  template <typename XSpanT, typename YSpanT, typename ArgT = dof_value_t,
            typename hyNode_index_t = dof_index_t>
  void residual(const XSpanT& x_vec, YSpanT& vec_Ax, const ArgT arg = 0.)
  {
    hy_assert(x_vec.size() == vec_Ax.size(), "x_vec and vec_Ax need to be of same size");
    const ArgT& time = arg;
    prototype_mat_vec_multiply_span(residual_flux);
    zero_span(vec_Ax);
  }
  /*!***********************************************************************************************
   * \brief   Assemble the condensed system operator as a COO matrix into the caller's buffer.
   *
   * Reusing one \c sparse_mat across repeated assemblies avoids the allocation. Configured zero
   * rows and columns are blanked with a single 1 on each owned diagonal.
   ************************************************************************************************/
  template <typename MatVecT = LargeVecT, typename ArgT = dof_value_t,
            typename hyNode_index_t = dof_index_t>
  void jacobian_mat(sparse_mat<MatVecT>& mat, const ArgT arg = 0.)
  {
    const ArgT& time = arg;
    prototype_mat_generate_into(trace_to_flux, mat);
    zero_mat(mat);
  }

  // -----------------------------------------------------------------------------------------------
  // Output: errors, norms, plotting.
  // -----------------------------------------------------------------------------------------------
  /*!***********************************************************************************************
   * \brief   Returns a zero vector of system size (explicit factory; the only allocating entry).
   ************************************************************************************************/
  LargeVecT zero_vector() const { return LargeVecT(hyper_graph_.n_local_dofs(), 0.); }
  /*!***********************************************************************************************
   * \brief   Calculate L2 error at the given time.
   ************************************************************************************************/
  template <typename SpanT, typename hyNode_index_t = dof_index_t>
  auto errors(const SpanT& x_vec, const double time = 0.)
  {
    auto result = prototype_errors(errors);
    return std::vector<typename decltype(result)::value_type>(result.begin(), result.end());
  }
  /*!***********************************************************************************************
   * \brief   Calculate L2 norm of the analytic solution at the given time.
   ************************************************************************************************/
  template <typename SpanT, typename hyNode_index_t = dof_index_t>
  auto norms(const SpanT& x_vec, const double time = 0.)
  {
    auto result = prototype_errors(norms);
    return std::vector<typename decltype(result)::value_type>(result.begin(), result.end());
  }
  /*!***********************************************************************************************
   * \brief   Set plot option and return old plot option.
   ************************************************************************************************/
  std::string plot_option(const std::string& option, std::string value = "")
  {
    return set_plot_option(plot_options, option, value);
  }
  /*!***********************************************************************************************
   * \brief   Plot solution.
   ************************************************************************************************/
  template <typename SpanT>
  void plot_solution(const SpanT& lambda, const double time = 0.)
  {
    plot(hyper_graph_, local_solver_, lambda, plot_options, time);
  }

  // -----------------------------------------------------------------------------------------------
  // Gauss stage introspection (stage counts only; stage solves ride through the generic entries
  // with a Gauss::StageTime argument, stage == -1 turns postprocess into the recombination).
  // -----------------------------------------------------------------------------------------------
  static constexpr bool gauss_stage_capable()
  {
    return requires { LocalSolverT::n_gauss_stages(); };
  }
  static constexpr unsigned int n_gauss_stages()
  {
    if constexpr (gauss_stage_capable())
      return LocalSolverT::n_gauss_stages();
    else
      return 1;
  }
  static constexpr unsigned int n_gauss_reps()
  {
    if constexpr (requires { LocalSolverT::n_gauss_reps(); })
      return LocalSolverT::n_gauss_reps();
    else
      return 1;
  }

  // -----------------------------------------------------------------------------------------------
  // Size queries and distributed-memory information.
  // -----------------------------------------------------------------------------------------------
  dof_index_t size_of_system() const { return hyper_graph_.n_global_dofs(); }
  dof_index_t n_owned_dofs() const { return hyper_graph_.n_owned_dofs(); }
  dof_index_t n_local_dofs() const { return hyper_graph_.n_local_dofs(); }
  dof_index_t n_edges() const { return hyper_graph_.n_hyEdges(); }
  std::vector<dof_index_t> local_to_global_dofs() const
  {
    return hyper_graph_.local_to_global_dofs();
  }
  static constexpr unsigned int space_dim() { return TopologyT::space_dim(); }
  /*!***********************************************************************************************
   * \brief   Flat coordinates of this rank's owned hypernodes; see
   *          \c HDGHyperGraph::owned_point_coords().
   ************************************************************************************************/
  std::vector<double> owned_point_coords() const { return hyper_graph_.owned_point_coords(); }
  /*!***********************************************************************************************
   * \brief   This rank's owned hyperedges as global hypernode index pairs; see
   *          \c HDGHyperGraph::owned_edges_global().
   ************************************************************************************************/
  std::vector<dof_index_t> owned_edges_global() const { return hyper_graph_.owned_edges_global(); }
  /*!***********************************************************************************************
   * \brief   Return refinement level.
   ************************************************************************************************/
  unsigned int get_refinement() { return hyper_graph_.get_refinement(); }
  /*!***********************************************************************************************
   * \brief   Set refinement level. Clears the zero index sets (the dof layout changes);
   *          re-run configure afterwards.
   ************************************************************************************************/
  void set_refinement(unsigned int level)
  {
    hyper_graph_.set_refinement(level);
    zero_indices_local_.clear();
    zero_indices_global_.clear();
    zero_indices_diag_.clear();
  }
};  // end of class Generic

}  // end of namespace GlobalLoop
