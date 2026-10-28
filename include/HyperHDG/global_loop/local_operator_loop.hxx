#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/global_loop/kinds.hxx>
#include <HyperHDG/global_loop/prototype.hxx>  // sparse_mat
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/plot.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Prototype of the single global loop built on local operators.
 *
 * The loop is a scatter--gather engine over hyperedges. For every hyperedge it gathers the local
 * trace from a global vector, hands it to the hyperedge's local operator, and scatters the result.
 * It never sees block matrices, factorizations, or boundary conditions of the local solver.
 *
 * Local solver protocol (cf. local_operator_design_v2.md):
 * - \c Args: runtime parameters (time, ...) forwarded to every call,
 * - \c local_operator(edge, Stiffness, args): builds the immutable per-edge operator with
 *   \c apply<k>(edge, in, args, out) and \c residual(edge, lambda, load, args, r),
 * - \c local_load(edge, args): the local right-hand side (loads, Dirichlet data) used by residual,
 * - optional \c operator_key(edge, Stiffness, args): share one operator among equal edges,
 * - optional \c constrained_node<space_dim>(node_type): default predicate for constrain_dirichlet.
 *
 * Store: after \c prepare the operators live in a contiguous vector (one per distinct key) until
 * \c release. Without prepare every sweep streams: build, use, drop per hyperedge.
 *
 * Constraints: a sorted index set of global dofs owned by the loop. apply is the identity there,
 * residual returns x_i (prescribed values are zero in this prototype, Dirichlet data enters via
 * the local load), assemble replaces the rows by zero or identity rows.
 *
 * Parallelism: with OpenMP, sweeps run in parallel over hyperedges with per-thread scatter
 * buffers summed in thread order, so results are deterministic for a fixed thread count.
 **************************************************************************************************/
template <class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename LargeVecT = std::vector<double>,
          typename dof_index_t = unsigned int>
class LocalOperatorLoop
{
 public:
  static constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
  static constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();
  static constexpr unsigned int n_loc = 2 * hyEdge_dim * n_dofs_per_node;
  using dof_value_t = typename LargeVecT::value_type;
  using Args = typename LocalSolverT::Args;
  using HyperGraphT = HDGHyperGraph<n_dofs_per_node,
                                    TopologyT,
                                    GeometryT,
                                    NodeDescriptorT,
                                    typename LocalSolverT::data_type>;
  // The hypergraph keeps its hyEdge type private; name it through operator[].
  using Edge = decltype(std::declval<HyperGraphT&>()[(unsigned int)0]);
  using Operator = std::remove_cvref_t<decltype(std::declval<const LocalSolverT&>().local_operator(
    std::declval<Edge&>(), Stiffness{}, std::declval<const Args&>()))>;
  template <unsigned int k>
  using Traces = SmallMat<n_loc, k, dof_value_t>;

 private:
  HyperGraphT hyper_graph_;
  const LocalSolverT local_solver_;
  PlotOptions plot_options_;

  std::vector<Operator> store_;             // one operator per distinct key (prepared only)
  std::vector<std::size_t> slot_;           // hyperedge -> index into store_
  std::vector<std::size_t> representative_; // slot -> a hyperedge with that key
  bool prepared_ = false;

  std::vector<dof_index_t> constraints_;  // sorted global dof indices
  std::vector<dof_value_t> free_mask_;     // 1 on free dofs, 0 on constrained dofs

 public:
  LocalOperatorLoop(const typename TopologyT::constructor_value_type& construct_topo,
                    const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo), local_solver_(construct_loc_sol)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
    static_assert(std::is_trivially_copyable_v<Operator>,
                  "Local operators must be trivially copyable (contiguous store, device copy).");
  }
  LocalOperatorLoop(const typename TopologyT::constructor_value_type& construct_topo)
  : LocalOperatorLoop(construct_topo, typename LocalSolverT::constructor_value_type())
  {
  }

  // ---- sizes, refinement, plotting (as in the other loops) --------------------------------------

  dof_index_t size_of_system() const { return hyper_graph_.n_global_dofs(); }
  LargeVecT zero_vector() const { return LargeVecT(hyper_graph_.n_global_dofs(), 0.); }
  unsigned int get_refinement() { return hyper_graph_.get_refinement(); }
  void set_refinement(unsigned int level)
  {
    hyper_graph_.set_refinement(level);
    release();
    constraints_.clear();
    free_mask_.clear();
  }
  std::string plot_option(const std::string& option, std::string value = "")
  {
    return set_plot_option(plot_options_, option, value);
  }
  void plot_solution(const LargeVecT& lambda, const dof_value_t time = 0.)
  {
    plot(hyper_graph_, local_solver_, lambda, plot_options_, time);
  }

  // ---- constraints -------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Constrain every dof of every hypernode whose node type satisfies \c pred.
   ************************************************************************************************/
  template <class Pred>
  void constrain_nodes(Pred pred)
  {
    free_mask_.assign(hyper_graph_.n_global_dofs(), 1.);
    std::array<unsigned int, n_dofs_per_node> dofs;
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    for (std::size_t e = 0; e < n_edges; ++e)
    {
      Edge edge = hyper_graph_[e];
      const auto& hyNodes = edge.topology.get_hyNode_indices();
      for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
        if (pred(edge.node_descriptor[face]))
        {
          hyper_graph_.hyNode_factory().get_dof_indices(hyNodes[face], dofs);
          for (unsigned int d = 0; d < n_dofs_per_node; ++d)
            free_mask_[dofs[d]] = 0.;
        }
    }
    constraints_.clear();
    for (dof_index_t i = 0; i < free_mask_.size(); ++i)
      if (free_mask_[i] == 0.)
        constraints_.push_back(i);
  }
  /*!***********************************************************************************************
   * \brief   Constrain the Dirichlet nodes as the local solver's parameters define them.
   ************************************************************************************************/
  void constrain_dirichlet_nodes()
  {
    constrain_nodes([](unsigned int type)
                    { return LocalSolverT::template constrained_node<GeometryT::space_dim()>(type); });
  }
  void set_constraints(const std::vector<unsigned int>& indices)
  {
    free_mask_.assign(hyper_graph_.n_global_dofs(), 1.);
    for (unsigned int i : indices)
    {
      hy_assert(i < free_mask_.size(), "Constraint index out of range.");
      free_mask_[i] = 0.;
    }
    constraints_.clear();
    for (dof_index_t i = 0; i < free_mask_.size(); ++i)
      if (free_mask_[i] == 0.)
        constraints_.push_back(i);
  }
  const std::vector<dof_index_t>& constraints() const { return constraints_; }

  // ---- store -------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Build the operators of all hyperedges (one per distinct key) and keep them.
   *
   * Pass 1 (serial): keys and deduplication. Pass 2 (parallel over slots): factorizations. Slots
   * are distinct, so no scatter buffers are needed, which is why this is a parallel_for and not a
   * sweep.
   ************************************************************************************************/
  void prepare(const Args& args = Args())
  {
    release();
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    slot_.resize(n_edges);
    std::vector<std::size_t>& representative = representative_;

    if constexpr (requires(Edge & e) { local_solver_.operator_key(e, Stiffness{}, args); })
    {
      using Key = std::remove_cvref_t<decltype(local_solver_.operator_key(
        std::declval<Edge&>(), Stiffness{}, args))>;
      std::map<Key, std::size_t> slot_of_key;
      for (std::size_t e = 0; e < n_edges; ++e)
      {
        Edge edge = hyper_graph_[e];
        auto [it, inserted] = slot_of_key.emplace(
          local_solver_.operator_key(edge, Stiffness{}, args), representative.size());
        if (inserted)
          representative.push_back(e);
        slot_[e] = it->second;
      }
    }
    else
    {
      representative.resize(n_edges);
      for (std::size_t e = 0; e < n_edges; ++e)
        slot_[e] = representative[e] = e;
    }

    store_.resize(representative.size());
    try
    {
      parallel_for(representative.size(),
                   [&](std::size_t s)
                   {
                     Edge edge = hyper_graph_[representative[s]];
                     store_[s] = local_solver_.local_operator(edge, Stiffness{}, args);
                   });
    }
    catch (...)
    {
      release();
      throw;
    }
    prepared_ = true;
  }
  void release()
  {
    store_.clear();
    store_.shrink_to_fit();
    slot_.clear();
    representative_.clear();
    prepared_ = false;
  }
  bool prepared() const { return prepared_; }
  std::size_t n_operators() const { return store_.size(); }

  // ---- sweeps ------------------------------------------------------------------------------------
  // Every public sweep decides once whether the store is used, and instantiates a branch-free
  // per-edge body for that mode (cf. dispatch / visit).

  /*!***********************************************************************************************
   * \brief   y = J x with identity rows on constrained dofs.
   ************************************************************************************************/
  LargeVecT apply(const LargeVecT& x, const Args& args = Args())
  {
    hy_assert(x.size() == hyper_graph_.n_global_dofs(), "Vector has wrong size.");
    LargeVecT y(x.size(), 0.);
    dispatch(
      [&](auto mode)
      {
        sweep<mode.value>(y, args,
                          [&](Edge& edge, const Operator& op, const Dofs& dofs, LargeVecT& out)
                          {
                            Traces<1> in, res;
                            gather(dofs, x, in);
                            op.apply(edge, in, args, res);
                            scatter(dofs, res, out);
                          });
      });
    for (dof_index_t i : constraints_)
      y[i] = x[i];
    return y;
  }
  /*!***********************************************************************************************
   * \brief   r = F(x): flux balance on free dofs, x_i on constrained dofs (prescribed value zero).
   ************************************************************************************************/
  LargeVecT residual(const LargeVecT& x, const Args& args = Args())
  {
    hy_assert(x.size() == hyper_graph_.n_global_dofs(), "Vector has wrong size.");
    LargeVecT r(x.size(), 0.);
    dispatch(
      [&](auto mode)
      {
        sweep<mode.value>(r, args,
                          [&](Edge& edge, const Operator& op, const Dofs& dofs, LargeVecT& out)
                          {
                            Traces<1> in, res;
                            gather(dofs, x, in);
                            op.residual(edge, in, local_solver_.local_load(edge, args), args, res);
                            scatter(dofs, res, out);
                          });
      });
    for (dof_index_t i : constraints_)
      r[i] = x[i];
    return r;
  }
  /*!***********************************************************************************************
   * \brief   Triplets of J: every hyperedge's block is apply<n_loc> of the identity.
   *
   * Entries in constrained rows are multiplied by zero (branch-free mask); with \c identity_rows a
   * unit diagonal entry is appended for every constrained dof. Every hyperedge writes its own
   * block, so the result does not depend on the number of threads.
   ************************************************************************************************/
  sparse_mat<LargeVecT> assemble(const Args& args = Args(), const bool identity_rows = false)
  {
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    const std::size_t block = (std::size_t)n_loc * n_loc;
    sparse_mat<LargeVecT> result;
    result.value_vec.resize(n_edges * block + (identity_rows ? constraints_.size() : 0));
    result.row_vec.resize(result.value_vec.size());
    result.col_vec.resize(result.value_vec.size());
    if (free_mask_.size() != hyper_graph_.n_global_dofs())
      free_mask_.assign(hyper_graph_.n_global_dofs(), 1.);

    Traces<n_loc> identity;
    for (unsigned int i = 0; i < n_loc; ++i)
      identity(i, i) = 1.;

    auto write_block = [&](std::size_t e, const Dofs& dofs, const Traces<n_loc>& loc_mat)
    {
      std::size_t index = e * block;
      for (unsigned int col = 0; col < n_loc; ++col)
        for (unsigned int row = 0; row < n_loc; ++row, ++index)
        {
          result.row_vec[index] = dofs[row];
          result.col_vec[index] = dofs[col];
          result.value_vec[index] = loc_mat(row, col) * free_mask_[dofs[row]];
        }
    };

    if (prepared_ && store_.size() < n_edges)
    {
      // Shared operators: one dense block per slot, then a copy per hyperedge.
      std::vector<Traces<n_loc> > blocks(store_.size());
      parallel_for(store_.size(),
                   [&](std::size_t s)
                   {
                     Edge edge = hyper_graph_[representative_[s]];
                     store_[s].apply(edge, identity, args, blocks[s]);
                   });
      parallel_for(n_edges,
                   [&](std::size_t e)
                   {
                     Dofs dofs;
                     local_dof_indices(hyper_graph_[e], dofs);
                     write_block(e, dofs, blocks[slot_[e]]);
                   });
    }
    else
      dispatch(
        [&](auto mode)
        {
          parallel_for(n_edges,
                       [&](std::size_t e)
                       {
                         Edge edge = hyper_graph_[e];
                         Dofs dofs;
                         local_dof_indices(edge, dofs);
                         Traces<n_loc> loc_mat;
                         visit<mode.value>(e, edge, args, [&](const Operator& op)
                                           { op.apply(edge, identity, args, loc_mat); });
                         write_block(e, dofs, loc_mat);
                       });
        });

    if (identity_rows)
    {
      std::size_t index = n_edges * block;
      for (dof_index_t i : constraints_)
      {
        result.row_vec[index] = i;
        result.col_vec[index] = i;
        result.value_vec[index++] = 1.;
      }
    }
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Fold f(edge, operator, local trace) over all hyperedges (serial).
   ************************************************************************************************/
  template <class F, class R, class Combine = std::plus<> >
  R reduce(const LargeVecT& x, const Args& args, F f, R result, Combine combine = Combine())
  {
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    dispatch(
      [&](auto mode)
      {
        for (std::size_t e = 0; e < n_edges; ++e)
        {
          Edge edge = hyper_graph_[e];
          Dofs dofs;
          local_dof_indices(edge, dofs);
          Traces<1> in;
          gather(dofs, x, in);
          visit<mode.value>(e, edge, args, [&](const Operator& op)
                            { result = combine(result, f(edge, op, in)); });
        }
      });
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Errors as the local solver defines them (reduce + square root), cf. Elliptic::errors.
   ************************************************************************************************/
  std::vector<dof_value_t> errors(const LargeVecT& x, const Args& args = Args())
  {
    using error_t = std::remove_cvref_t<decltype(local_solver_.error_sq(
      std::declval<Edge&>(), std::declval<const Operator&>(), std::declval<const Traces<1>&>(),
      args))>;
    error_t init;
    init.fill(0.);
    error_t sum = reduce(
      x, args, [&](Edge& edge, const Operator& op, const Traces<1>& lam)
      { return local_solver_.error_sq(edge, op, lam, args); }, init,
      [](error_t a, const error_t& b)
      {
        for (std::size_t k = 0; k < a.size(); ++k)
          a[k] += b[k];
        return a;
      });
    std::vector<dof_value_t> result(sum.begin(), sum.end());
    for (auto& v : result)
      v = std::sqrt(v);
    return result;
  }

 private:
  using Dofs = std::array<dof_index_t, n_loc>;

  void local_dof_indices(const Edge& edge, Dofs& dofs) const
  {
    const auto& hyNodes = edge.topology.get_hyNode_indices();
    std::array<dof_index_t, n_dofs_per_node> node_dofs;
    for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
    {
      hyper_graph_.hyNode_factory().get_dof_indices(hyNodes[face], node_dofs);
      for (unsigned int d = 0; d < n_dofs_per_node; ++d)
        dofs[face * n_dofs_per_node + d] = node_dofs[d];
    }
  }
  static void gather(const Dofs& dofs, const LargeVecT& x, Traces<1>& loc)
  {
    for (unsigned int i = 0; i < n_loc; ++i)
      loc[i] = x[dofs[i]];
  }
  static void scatter(const Dofs& dofs, const Traces<1>& loc, LargeVecT& y)
  {
    for (unsigned int i = 0; i < n_loc; ++i)
      y[dofs[i]] += loc[i];
  }
  /*!***********************************************************************************************
   * \brief   Call f with a compile-time bool (std::true_type / std::false_type) telling whether
   *          the store is used. The runtime decision is taken once per public call.
   ************************************************************************************************/
  template <class F>
  void dispatch(F&& f) const
  {
    if (prepared_)
      f(std::true_type{});
    else
      f(std::false_type{});
  }
  /*!***********************************************************************************************
   * \brief   Hand the operator of hyperedge e to f: from the store, or built on the fly.
   ************************************************************************************************/
  template <bool use_store, class F>
  void visit(const std::size_t e, Edge& edge, const Args& args, F&& f) const
  {
    if constexpr (use_store)
      f(store_[slot_[e]]);
    else
    {
      const Operator op = local_solver_.local_operator(edge, Stiffness{}, args);
      f(op);
    }
  }
  /*!***********************************************************************************************
   * \brief   Run body(i) for i in [0, n), in parallel with OpenMP. The first exception thrown by a
   *          body is rethrown after the loop (an exception leaving an OpenMP region would
   *          terminate the program).
   ************************************************************************************************/
  template <class Body>
  static void parallel_for(const std::size_t n, Body body)
  {
    std::exception_ptr failure;
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (long long i = 0; i < (long long)n; ++i)
    {
      try
      {
        body((std::size_t)i);
      }
      catch (...)
      {
#ifdef _OPENMP
#pragma omp critical(hyperhdg_local_operator_loop)
#endif
        if (!failure)
          failure = std::current_exception();
      }
    }
    if (failure)
      std::rethrow_exception(failure);
  }
  /*!***********************************************************************************************
   * \brief   Run body(edge, operator, dofs, out) over all hyperedges, scattering into y.
   *
   * With OpenMP every thread scatters into its own buffer; the buffers are summed in thread order
   * afterwards (deterministic for a fixed thread count). Memory: n_threads * n_dofs, which is
   * fine for a prototype and not for very large graphs (coloring or node-centric scatter then).
   ************************************************************************************************/
  template <bool use_store, class Body>
  void sweep(LargeVecT& y, const Args& args, Body body)
  {
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
#ifdef _OPENMP
    std::vector<LargeVecT> buffers(omp_get_max_threads(), LargeVecT(y.size(), 0.));
    auto out_of = [&](std::size_t) -> LargeVecT& { return buffers[omp_get_thread_num()]; };
#else
    auto out_of = [&](std::size_t) -> LargeVecT& { return y; };
#endif
    parallel_for(n_edges,
                 [&](std::size_t e)
                 {
                   Edge edge = hyper_graph_[e];
                   Dofs dofs;
                   local_dof_indices(edge, dofs);
                   visit<use_store>(e, edge, args, [&](const Operator& op)
                                    { body(edge, op, dofs, out_of(e)); });
                 });
#ifdef _OPENMP
    for (const LargeVecT& buffer : buffers)
      for (std::size_t i = 0; i < y.size(); ++i)
        y[i] += buffer[i];
#endif
  }
};  // end of class LocalOperatorLoop

}  // end of namespace GlobalLoop
