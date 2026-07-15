#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/global_loop/prototype.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/plot.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Combine local solver and global information for parabolic problems.
 *
 * \tparam  TopologyT       Class type containing topological information.
 * \tparam  GeometryT       Class type containing geometrical information.
 * \tparam  NodeDescriptorT Class type containing the information of nodes of hyperedges.
 * \tparam  LocalSolverT    Class type of the local solver.
 * \tparam  LargeVecT       Clas type of large, global vector.
 * \tparam  dof_index_t     Index type of hyperedges. Default is \c unsigned \c int.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename LargeVecT = std::vector<double>,
          typename dof_index_t = unsigned int>
class Hyperbolic
{
 public:
  /*!***********************************************************************************************
   * \brief   Some constant variable that might be helpful.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
  /*!***********************************************************************************************
   * \brief   Some constant variable that might be helpful.
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
   * \brief   Instantiation of a local solver.
   ************************************************************************************************/
  const LocalSolverT local_solver_;
  /*!***********************************************************************************************
   * \brief   Struct encoding the options for plotting.
   ************************************************************************************************/
  PlotOptions plot_options;

 public:
  /*!***********************************************************************************************
   * \brief   Abstract problem constructor.
   *
   * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
   * hyperedge.
   *
   * \param   construct_topo    Information to construct a topology.
   * \param   construct_geom    Information to construct a geometry.
   * \param   construct_loc_sol Information to construct a local solver.
   ************************************************************************************************/
  Hyperbolic(const typename TopologyT::constructor_value_type& construct_topo,
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
   * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
   * hyperedge.
   *
   * \param   construct_topo    Information to construct a topology.
   * \param   construct_loc_sol Information to construct a local solver.
   ************************************************************************************************/
  Hyperbolic(const typename TopologyT::constructor_value_type& construct_topo,
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
   * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
   * hyperedge.
   *
   * \param   construct_topo    Information to construct a topology.
   ************************************************************************************************/
  Hyperbolic(const typename TopologyT::constructor_value_type& construct_topo)
  : hyper_graph_(construct_topo)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Returns vector of appropriate size for the predefined problem.
   *
   * Returns a vector containing only the value zero, but of the size \f$n\f$ which is also the
   * number which is returned if \c size_of_system() is evaluated.
   *
   * \retval  zero          A vector of the correct size for the unknowns of the given problem.
   ************************************************************************************************/
  LargeVecT zero_vector() const { return LargeVecT(hyper_graph_.n_local_dofs(), 0.); }
  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product.
   *
   * Function that evaluates the condensed, matrix-free version of the matrix-vector product
   * \f$A x = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method that needs to be
   * inverted for a time step, \f$x\f$ is the vector of parameters to define the skeletal variable
   * \f$\lambda\f$, and \f$y\f$ is the resulting vector, which has the same size as the input vector
   * \f$x\f$.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time at which the new time step will end.
   * \retval  y_vec         A vector containing the product \f$y = Ax\f$.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT trace_to_flux(const LargeVecT& x_vec, const dof_value_t time = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    LargeVecT vec_Ax(x_vec.size(), 0.);
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs_old,
      hyEdge_dofs_new;

    // Do matrix--vector multiplication by iterating over all hyperedges.
    std::for_each(
      hyper_graph_.begin(), hyper_graph_.end(),
      [&](auto hyper_edge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        {
          hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec,
                                                       hyEdge_dofs_old[hyNode]);
          hyEdge_dofs_new[hyNode].fill(0.);
        }

        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr (requires { local_solver_.trace_to_flux(hyEdge_dofs_old, hyEdge_dofs_new, prvalue_of(time)); })
          local_solver_.trace_to_flux(hyEdge_dofs_old, hyEdge_dofs_new, time);
        else if constexpr (requires { local_solver_.trace_to_flux(hyEdge_dofs_old, hyEdge_dofs_new, hyper_edge, prvalue_of(time)); })
          local_solver_.trace_to_flux(hyEdge_dofs_old, hyEdge_dofs_new, hyper_edge, time);
        else
          static_assert(always_false_v<LocalSolverT, decltype(hyper_edge)>,
                        "LocalSolverT implements no usable overload of trace_to_flux!");

        // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
          hyper_graph_.hyNode_factory().add_to_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                          hyEdge_dofs_new[hyNode]);
      });

    return vec_Ax;
  }


  template <typename hyNode_index_t = dof_index_t,
            typename MatVecT = LargeVecT,
            typename time_t = dof_value_t>
  sparse_mat<MatVecT> trace_to_flux_mat(const time_t time = 0.)
  {
    return prototype_mat_generate(trace_to_flux, MatVecT);
  }
 
  template <typename hyNode_index_t = dof_index_t, typename SpanT, typename time_t = dof_value_t>
  void residual_flux2(const SpanT& x_vec, SpanT& vec_Ax, time_t time = 0.) {
    hy_assert(x_vec.size() == vec_Ax.size(), "x_vec and vec_Ax need to be of same size");
    prototype_mat_vec_multiply_span(residual_flux);
  }


  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product.
   *
   * Function that evaluates the condensed, matrix-free version of the matrix-vector product
   * \f$A x = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method that needs to be
   * inverted to do a time step, \f$x\f$ is the vector of parameters to define the skeletal variable
   * \f$\lambda\f$, and \f$y\f$ is the resulting vector, which has the same size as the input vector
   * \f$x\f$.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time at which the time step ends.
   * \retval  y_vec         A vector containing the product \f$y = Ax\f$.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  std::vector<dof_value_t> residual_flux(const std::vector<dof_value_t>& x_vec,
                                         const dof_value_t time = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    LargeVecT vec_Ax(x_vec.size(), 0.);
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs_old,
      hyEdge_dofs_new;

    // Do matrix--vector multiplication by iterating over all hyperedges.
    std::for_each(
      hyper_graph_.begin(), hyper_graph_.end(),
      [&](auto hyper_edge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        {
          hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec,
                                                       hyEdge_dofs_old[hyNode]);
          hyEdge_dofs_new[hyNode].fill(0.);
        }

        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr (requires { local_solver_.residual_flux(hyEdge_dofs_old, hyEdge_dofs_new, prvalue_of(time)); })
          local_solver_.residual_flux(hyEdge_dofs_old, hyEdge_dofs_new, time);
        else if constexpr (requires { local_solver_.residual_flux(hyEdge_dofs_old, hyEdge_dofs_new, hyper_edge, prvalue_of(time)); })
          local_solver_.residual_flux(hyEdge_dofs_old, hyEdge_dofs_new, hyper_edge, time);
        else
          static_assert(always_false_v<LocalSolverT, decltype(hyper_edge)>,
                        "LocalSolverT implements no usable overload of residual_flux!");

        // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
          hyper_graph_.hyNode_factory().add_to_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                          hyEdge_dofs_new[hyNode]);
      });

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Set data using the result of the old time step.
   *
   * \param   x_vec         A \c std::vector containing the input vector \f$x\f$.
   * \param   time          Time at which the old time step ended.
   ************************************************************************************************/
  template <typename SpanT, typename hyNode_index_t = dof_index_t, typename time_t = dof_value_t>
  void set_data(const SpanT& x_vec, const time_t time = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();
    using span_value_t = typename std::decay_t<SpanT>::value_type;

    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<span_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

    // Do matrix--vector multiplication by iterating over all hyperedges.
    std::for_each(
      hyper_graph_.begin(), hyper_graph_.end(),
      [&](auto hyper_edge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
          hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec,
                                                       hyEdge_dofs[hyNode]);

        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr (requires { local_solver_.set_data(hyEdge_dofs, prvalue_of(time)); })
          local_solver_.set_data(hyEdge_dofs, time);
        else if constexpr (requires { local_solver_.set_data(hyEdge_dofs, hyper_edge, prvalue_of(time)); })
          local_solver_.set_data(hyEdge_dofs, hyper_edge, time);
        else
          static_assert(always_false_v<LocalSolverT, decltype(hyper_edge)>,
                        "LocalSolverT implements no usable overload of set_data!");
      });
  }
  /*!***********************************************************************************************
   * \brief   Gauss step protocol: combine stashed stage data into the new state on every edge.
   *
   * Explicit driver-called completion of a Gauss time step (after all per-stage set_data calls).
   * Optional capability: solvers without it fail loudly at runtime, not at compile time.
   ************************************************************************************************/
  void finalize_step()
  {
    std::for_each(hyper_graph_.begin(), hyper_graph_.end(),
                  [&](auto hyper_edge)
                  {
                    if constexpr (requires { local_solver_.finalize_step(hyper_edge); })
                      local_solver_.finalize_step(hyper_edge);
                    else
                      hy_check(false, "LocalSolverT implements no finalize_step!");
                  });
  }
  // -----------------------------------------------------------------------------------------------
  // Gauss stage introspection: stage counts and endpoint-recombination weights. The stage SOLVES
  // reuse the generic entries above (trace_to_flux_mat / residual_flux2 / set_data), instantiated
  // with complex vectors and a Gauss::StageTime as the time argument (the local solver unpacks
  // time and stage index from it).
  // -----------------------------------------------------------------------------------------------
  /*!***********************************************************************************************
   * \brief   Number of Gauss collocation stages / solved stage representatives of the solver.
   ************************************************************************************************/
  static constexpr unsigned int n_gauss_stages()
  {
    if constexpr (requires { LocalSolverT::n_gauss_stages(); })
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
  /*!***********************************************************************************************
   * \brief   Endpoint-update weights for the driver's trace recombination.
   ************************************************************************************************/
  void stage_weights(const unsigned int rep,
                     dof_value_t& affine,
                     dof_value_t& mult,
                     dof_value_t& w_re,
                     dof_value_t& w_im) const
  {
    if constexpr (requires {
                    local_solver_.stage_affine();
                    local_solver_.stage_mult(rep);
                    local_solver_.stage_w(rep);
                  })
    {
      affine = local_solver_.stage_affine();
      mult = local_solver_.stage_mult(rep);
      w_re = local_solver_.stage_w(rep).real();
      w_im = local_solver_.stage_w(rep).imag();
    }
    else
      hy_check(false, "LocalSolverT exposes no stage weights!");
  }
  /*!***********************************************************************************************
   * \brief   Evaluate the initial flux of the problem.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time for initial data.
   * \retval  y_vec         A vector containing the initial fluxes.
   ************************************************************************************************/
  template <typename SpanT, typename hyNode_index_t = dof_index_t>
  void make_initial_from_static(const SpanT& x_vec, const dof_index_t time = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

    // Do matrix--vector multiplication by iterating over all hyperedges.
    std::for_each(
      hyper_graph_.begin(), hyper_graph_.end(),
      [&](auto hyper_edge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
          for (unsigned int d = 0; d < n_dofs_per_node; ++d)
            hyEdge_dofs[hyNode][d] = x_vec[hyEdge_hyNodes[hyNode] * n_dofs_per_node + d];

        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr (requires { local_solver_.make_initial_from_static(hyEdge_dofs, prvalue_of(time)); })
          local_solver_.make_initial_from_static(hyEdge_dofs, time);
        else if constexpr (requires { local_solver_.make_initial_from_static(hyEdge_dofs, hyper_edge, prvalue_of(time)); })
          local_solver_.make_initial_from_static(hyEdge_dofs, hyper_edge, time);
        else
          // optional capability (static-init hack): fail loudly at runtime instead of at compile
          // time, so solvers without it still instantiate this loop
          hy_check(false, "LocalSolverT implements no usable overload of make_initial_from_static!");

        // write per-edge dofs back to the global trace vector (mirrors make_initial)
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
          hyper_graph_.hyNode_factory().set_dof_values(hyEdge_hyNodes[hyNode], x_vec,
                                                       hyEdge_dofs[hyNode]);
      });
  }
  /*!***********************************************************************************************
   * \brief   Evaluate the initial flux of the problem.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time for initial data.
   * \retval  y_vec         A vector containing the initial fluxes.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT make_initial(const LargeVecT& x_vec, const dof_index_t time = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    std::vector<dof_value_t> vec_Ax(x_vec.size(), 0.);
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

    // Do matrix--vector multiplication by iterating over all hyperedges.
    std::for_each(
      hyper_graph_.begin(), hyper_graph_.end(),
      [&](auto hyper_edge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
          hyEdge_dofs[hyNode].fill(0.);

        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr (requires { local_solver_.make_initial(hyEdge_dofs, prvalue_of(time)); })
          local_solver_.make_initial(hyEdge_dofs, time);
        else if constexpr (requires { local_solver_.make_initial(hyEdge_dofs, hyper_edge, prvalue_of(time)); })
          local_solver_.make_initial(hyEdge_dofs, hyper_edge, time);
        else
          static_assert(always_false_v<LocalSolverT, decltype(hyper_edge)>,
                        "LocalSolverT implements no usable overload of make_initial!");

        // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
        for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
          hyper_graph_.hyNode_factory().set_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                       hyEdge_dofs[hyNode]);
      });

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Calculate L2 error.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time at which error is evaluated.
   * \retval  error         L2 error.
   ************************************************************************************************/
  template <typename SpanT, typename hyNode_index_t = dof_index_t>
  std::vector<dof_value_t> errors(const SpanT& x_vec, const dof_value_t time = 0.)
  {
    auto result = prototype_errors(errors);
    return std::vector<dof_value_t>(result.begin(), result.end());
  }
  /*!***********************************************************************************************
   * \brief   Calculate L2 norm of the analytic solution.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time at which norm is evaluated.
   * \retval  norm          L2 norm of the analytic solution.
   ************************************************************************************************/
  template <typename SpanT, typename hyNode_index_t = dof_index_t>
  std::vector<dof_value_t> norms(const SpanT& x_vec, const dof_value_t time = 0.)
  {
    auto result = prototype_errors(norms);
    return std::vector<dof_value_t>(result.begin(), result.end());
  }
  /*!***********************************************************************************************
   * \brief   Number of hyperedges in the underlying hypergraph.
   ************************************************************************************************/
  dof_index_t n_edges() const { return hyper_graph_.n_hyEdges(); }
  /*!***********************************************************************************************
   * \brief   Number of energy components reported per edge by the local solver.
   ************************************************************************************************/
  static constexpr unsigned int n_energy_components()
  {
    return LocalSolverT::n_energy_components();
  }
  /*!***********************************************************************************************
   * \brief   Fill caller-provided span with per-edge energy components.
   *
   * Does not allocate. The span must hold exactly \c n_edges() * \c n_energy_components()
   * entries; the local solver's array for edge \c e is written into
   * \c out[e * n_energy_components() .. (e+1) * n_energy_components()).
   ************************************************************************************************/
  template <typename SpanT, typename LambdaSpanT>
  void energy(const LambdaSpanT& x_vec, SpanT out, const dof_value_t time = 0.)
  {
    constexpr unsigned int n_comp = n_energy_components();
    hy_check(out.size() == static_cast<std::size_t>(n_edges()) * n_comp,
             "energy span size mismatch: got " << out.size()
             << " expected " << static_cast<std::size_t>(n_edges()) * n_comp);

    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> dofs;

    const dof_index_t n_e = n_edges();
    for (dof_index_t e = 0; e < n_e; ++e)
    {
      auto hyper_edge = hyper_graph_[e];
      const auto hyNodes = hyper_edge.topology.get_hyNode_indices();
      for (unsigned int node = 0; node < hyNodes.size(); ++node)
        hyper_graph_.hyNode_factory().get_dof_values(hyNodes[node], x_vec, dofs[node]);

      const auto local = local_solver_.energy(dofs, hyper_edge, time);
      const std::size_t base = static_cast<std::size_t>(e) * n_comp;
      for (unsigned int c = 0; c < n_comp; ++c)
        out[base + c] = local[c];
    }
  }
  /*!***********************************************************************************************
   * \brief   Determine size of condensed system for the skeletal unknowns.
   *
   * Function that returns the size \f$n\f$ of the \f$n \times n\f$ linear, sparse system
   * \f$Ax = b\f$ that is solved by the program in a matrix-free fashion.
   *
   * This function is needed to define a \c LinearOperator from Python's \c scipy.sparse.linalg
   * package which can be used to define iterative solvers for sparse systems.
   *
   * \retval  n             Size of condensed system of equations.
   ************************************************************************************************/
  dof_index_t size_of_system() const { return hyper_graph_.n_global_dofs(); }

  dof_index_t n_owned_dofs() const { return hyper_graph_.n_owned_dofs(); }

  dof_index_t n_local_dofs() const { return hyper_graph_.n_local_dofs(); }

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
   * \brief   Set plot option and return old plot option.
   *
   * Function to set and / or read the current plot option.
   *
   * \param   option        A \c std::string containing the plot option to be considered.
   * \param   value         A \c std::string containing the new value of the considered option.
   *                        If empty, the old value is kept.
   * \retval  opt_value     A \c std::string containing the value of the plot option.
   ************************************************************************************************/
  std::string plot_option(const std::string& option, std::string value = "")
  {
    return set_plot_option(plot_options, option, value);
  }
  /*!***********************************************************************************************
   * \brief   Plot solution in vtu format.
   *
   * Function that plots the solution of the problem to a predefined file.
   *
   * \param   lambda        A vector of unknowns containing the data vector.
   * \param   time          Time at which analytic functions are evaluated.
   * \retval  file          A file in the output directory.
   ************************************************************************************************/
  template<typename SpanT>
  void plot_solution(const SpanT& lambda, const dof_value_t time = 0.)
  {
    plot(hyper_graph_, local_solver_, lambda, plot_options, time);
  }
  /*!***********************************************************************************************
   * \brief   Return refinement level.
   ************************************************************************************************/
  unsigned int get_refinement() { return hyper_graph_.get_refinement(); }
  /*!***********************************************************************************************
   * \brief   Set refinement level.
   ************************************************************************************************/
  void set_refinement(unsigned int level) { hyper_graph_.set_refinement(level); }
};  // end of class Hyperbolic

}  // end of namespace GlobalLoop
