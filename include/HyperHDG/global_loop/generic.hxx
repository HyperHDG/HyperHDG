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
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace GlobalLoop
{
namespace detail
{
// The kinds of operators of a local solver: its Kinds (a std::tuple), else Stiffness only.
template <class SolverT>
struct kinds
{
  using type = std::tuple<Stiffness>;
};
template <class SolverT>
  requires requires { typename SolverT::Kinds; }
struct kinds<SolverT>
{
  using type = typename SolverT::Kinds;
};
template <class SolverT>
using kinds_of = typename kinds<SolverT>::type;
}  // end of namespace detail

/*!*************************************************************************************************
 * \brief   The single global loop built on local operators.
 *
 * The loop is a scatter--gather engine over hyperedges. For every hyperedge it gathers the local
 * trace from a global vector, hands it to the hyperedge's local operator, and scatters the result.
 * It never sees block matrices, factorizations, or boundary conditions of the local solver, and
 * it never allocates global vectors: every entry point writes into caller-provided vectors (any
 * type with size() and operator[], e.g. std::vector or std::span).
 *
 * Local solver protocol (cf. doc/local_operator_design.md):
 * - \c Args: runtime parameters (time, ...) forwarded to every call,
 * - \c local_operator(edge, Kind, args): builds the immutable operator of that kind with
 *   \c apply<k>(solver, edge, in, args, out) and \c residual(solver, edge, lambda, args, r); what
 *   it computes and keeps is the solver's business (e.g. local_solver/local_operator.hxx). The
 *   loop hands the solver and the hyperedge to every call of the operator,
 * - optional \c n_operator_slots(edge, Kind, args) and \c operator_slot(edge, Kind, args): the
 *   solver numbers its operators of that kind, and hyperedges with the same number share one,
 * - optional \c Kinds: a std::tuple of the kinds of operators the solver builds (kinds.hxx), if
 *   there are others than Stiffness,
 * - optional \c constrained_node<space_dim>(node_type): default predicate for constrain_dirichlet.
 *
 * Kinds: \c prepare and the sweeps are templates on the kind of operator, Stiffness by default
 * (loop.apply<Mass>(x, y)); each kind has its own store, prepared on its own.
 *
 * Store: after \c prepare the operators live in a contiguous vector until \c release, one per
 * slot; without operator_slot, the slot of a hyperedge is its index. A slot names the operator
 * that is stored, not the local problem: hyperedges of one slot read the same operator, whose
 * results may still depend on the hyperedge it is applied to. Operators without state (an empty
 * type) are all the same, so there is one slot for all hyperedges. The loop keeps nothing per
 * hyperedge beyond the operators.
 * Without prepare every sweep builds the operator of a hyperedge, uses it and drops it. Every
 * public sweep takes that decision once and runs a branch-free per-edge body for the chosen mode.
 *
 * Constraints: a sorted index set of global dofs owned by the loop. Constrained dofs are not part
 * of the system (prescribed values are zero in this prototype, Dirichlet data enters via the
 * local load): every sweep reads x as zero there, and apply and residual return zero there. So
 * apply is P J P with the projection P onto the free dofs, symmetric if J is, and singular;
 * assemble gives the same matrix, with a unit diagonal on the constrained dofs on request.
 *
 * Serial. Parallel sweeps (scatter strategy) are a separate step.
 **************************************************************************************************/
template <class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename dof_value_t = double,
          typename dof_index_t = unsigned int>
class Generic
{
 public:
  static constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
  static constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();
  static constexpr unsigned int n_loc = 2 * hyEdge_dim * n_dofs_per_node;
  using value_type = dof_value_t;
  using index_type = dof_index_t;
  using Args = typename LocalSolverT::Args;
  using HyperGraphT = HDGHyperGraph<n_dofs_per_node,
                                    TopologyT,
                                    GeometryT,
                                    NodeDescriptorT,
                                    typename LocalSolverT::data_type>;
  // The hypergraph keeps its hyEdge type private; name it through operator[].
  using Edge = decltype(std::declval<HyperGraphT&>()[(unsigned int)0]);
  template <class Kind>
  using OperatorOf = std::remove_cvref_t<
    decltype(std::declval<const LocalSolverT&>()
               .local_operator(std::declval<Edge&>(), Kind{}, std::declval<const Args&>()))>;
  using Operator = OperatorOf<Stiffness>;
  template <unsigned int k>
  using Traces = SmallMat<n_loc, k, dof_value_t>;

 private:
  HyperGraphT hyper_graph_;
  const LocalSolverT local_solver_;

  template <class Kind>
  struct Store
  {
    static_assert(std::is_trivially_copyable_v<OperatorOf<Kind> >,
                  "Local operators must be trivially copyable (contiguous store, device copy).");
    std::vector<OperatorOf<Kind> > operators;  // one per slot (prepared only)
    bool prepared = false;
  };
  template <class... Kinds>
  static std::tuple<Store<Kinds>...> stores_of(std::tuple<Kinds...>*);
  decltype(stores_of((detail::kinds_of<LocalSolverT>*)nullptr)) stores_;  // one per kind

  PlotOptions plot_options_;

  std::vector<dof_index_t> constraints_;  // sorted global dof indices
  std::vector<dof_value_t> free_mask_;    // 1 on free dofs, 0 on constrained dofs

  static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                "Hyperedge dimension of topology and geometry must be equal!");
  static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                "Space dimension of topology and geometry must be equal!");
  static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                "Hyperedge dimension of hypergraph and local solver must be equal!");

 public:
  Generic(const typename TopologyT::constructor_value_type& construct_topo,
          const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo), local_solver_(construct_loc_sol)
  {
    free_mask_.assign(hyper_graph_.n_global_dofs(), 1.);
  }
  // The local solver's own defaults (e.g. tau = 1 for Diffusion), not a value-initialized
  // constructor argument (tau = 0, a singular local problem).
  Generic(const typename TopologyT::constructor_value_type& construct_topo)
  : hyper_graph_(construct_topo), local_solver_()
  {
    free_mask_.assign(hyper_graph_.n_global_dofs(), 1.);
  }

  // ---- sizes, refinement
  // ---------------------------------------------------------------------------

  const LocalSolverT& local_solver() const { return local_solver_; }
  dof_index_t size_of_system() const { return hyper_graph_.n_global_dofs(); }
  std::size_t n_hyEdges() const { return hyper_graph_.n_hyEdges(); }
  unsigned int get_refinement() { return hyper_graph_.get_refinement(); }
  void set_refinement(unsigned int level)
  {
    hyper_graph_.set_refinement(level);
    release();
    constraints_.clear();
    free_mask_.assign(hyper_graph_.n_global_dofs(), 1.);
  }

  // ---- constraints
  // -------------------------------------------------------------------------------

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
    collect_constraints();
  }
  /*!***********************************************************************************************
   * \brief   Constrain the Dirichlet nodes as the local solver's parameters define them.
   ************************************************************************************************/
  void constrain_dirichlet_nodes()
  {
    constrain_nodes(
      [](unsigned int type)
      { return LocalSolverT::template constrained_node<GeometryT::space_dim()>(type); });
  }
  template <class IndexVecT>
  void set_constraints(const IndexVecT& indices)
  {
    free_mask_.assign(hyper_graph_.n_global_dofs(), 1.);
    for (std::size_t k = 0; k < indices.size(); ++k)
    {
      hy_assert(indices[k] < free_mask_.size(), "Constraint index out of range.");
      free_mask_[indices[k]] = 0.;
    }
    collect_constraints();
  }
  const std::vector<dof_index_t>& constraints() const { return constraints_; }

  // ---- store
  // -------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Build the operators of a kind (one per slot) and keep them.
   ************************************************************************************************/
  template <class Kind = Stiffness>
  void prepare(const Args& args = Args())
  {
    release<Kind>();
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    std::vector<OperatorOf<Kind> >& store_ = store<Kind>().operators;
    store_.resize(n_operator_slots<Kind>(args));
    std::vector<bool> built(store_.size(), false);
    try
    {
      for (std::size_t e = 0; e < n_edges; ++e)
      {
        Edge edge = hyper_graph_[e];
        const std::size_t slot = operator_slot<Kind>(e, edge, args);
        hy_assert(slot < store_.size(), "The local solver has only "
                                          << store_.size() << " operator slots, but named slot "
                                          << slot << ".");
        if (!built[slot])
        {
          store_[slot] = local_solver_.local_operator(edge, Kind{}, args);
          built[slot] = true;
        }
      }
    }
    catch (...)
    {
      release<Kind>();
      throw;
    }
    store<Kind>().prepared = true;
  }
  /*!***********************************************************************************************
   * \brief   Drop the stored operators of a kind, or of all kinds.
   ************************************************************************************************/
  template <class Kind>
  void release()
  {
    store<Kind>() = Store<Kind>();
  }
  void release() { stores_ = {}; }
  template <class Kind = Stiffness>
  bool prepared() const
  {
    return store<Kind>().prepared;
  }
  template <class Kind = Stiffness>
  std::size_t n_operators() const
  {
    return store<Kind>().operators.size();
  }

  // ---- sweeps
  // ------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   y = J x on free dofs, with x read as zero on constrained dofs; zero there. y is
   *          overwritten.
   *
   * x is the direction of the directional derivative. Where a nonlinear problem is linearized is
   * not an argument: args reaches the local solver as it is, but a global vector does not (only
   * the loop can gather its local traces), so the local solver keeps the point, e.g. in the data
   * of the hyperedges.
   ************************************************************************************************/
  template <class Kind = Stiffness, class VecInT, class VecOutT>
  void apply(const VecInT& x, VecOutT& y, const Args& args = Args())
  {
    check_sizes(x, y);
    zero(y);
    dispatch<Kind>(
      [&](auto mode)
      {
        sweep<Kind, mode.value>(args,
                                [&](std::size_t, Edge& edge, const auto& op, const Dofs& dofs)
                                {
                                  Traces<1> in, res;
                                  gather(dofs, x, in);
                                  op.apply(local_solver_, edge, in, args, res);
                                  scatter(dofs, res, y);
                                });
      });
    for (dof_index_t i : constraints_)
      y[i] = 0.;
  }
  /*!***********************************************************************************************
   * \brief   r = F(x): flux balance on free dofs, with x read as zero on constrained dofs; zero
   *          there. r is overwritten.
   *
   * x is the point at which the residual is evaluated.
   ************************************************************************************************/
  template <class Kind = Stiffness, class VecInT, class VecOutT>
  void residual(const VecInT& x, VecOutT& r, const Args& args = Args())
  {
    check_sizes(x, r);
    zero(r);
    dispatch<Kind>(
      [&](auto mode)
      {
        sweep<Kind, mode.value>(args,
                                [&](std::size_t, Edge& edge, const auto& op, const Dofs& dofs)
                                {
                                  Traces<1> in, res;
                                  gather(dofs, x, in);
                                  op.residual(local_solver_, edge, in, args, res);
                                  scatter(dofs, res, r);
                                });
      });
    for (dof_index_t i : constraints_)
      r[i] = 0.;
  }
  /*!***********************************************************************************************
   * \brief   Number of triplets assemble() writes: one block per hyperedge plus the identity rows.
   ************************************************************************************************/
  std::size_t n_triplets(const bool identity_rows = false) const
  {
    return hyper_graph_.n_hyEdges() * (std::size_t)n_loc * n_loc +
           (identity_rows ? constraints_.size() : 0);
  }
  /*!***********************************************************************************************
   * \brief   Triplets of J into caller-provided arrays of size n_triplets(identity_rows).
   *
   * Every hyperedge's block is apply<n_loc> of the identity (so the matrix is the matrix-free
   * operator by construction). Entries in constrained rows and columns are multiplied by zero
   * (branch-free mask), so the matrix is symmetric if J is; with \c identity_rows a unit diagonal
   * entry is appended for every constrained dof, which makes it non-singular.
   * Duplicate (row, col) pairs have to be summed by the consumer (as scipy's coo_matrix does).
   ************************************************************************************************/
  template <class Kind = Stiffness, class IndexVecT, class ValueVecT>
  void assemble(const Args& args,
                const bool identity_rows,
                IndexVecT& rows,
                IndexVecT& cols,
                ValueVecT& vals)
  {
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    const std::size_t block = (std::size_t)n_loc * n_loc;
    hy_assert(rows.size() == n_triplets(identity_rows) && cols.size() == rows.size() &&
                vals.size() == rows.size(),
              "Triplet arrays must have size n_triplets(identity_rows).");
    Traces<n_loc> identity;
    for (unsigned int i = 0; i < n_loc; ++i)
      identity(i, i) = 1.;

    auto write_block = [&](std::size_t e, const Dofs& dofs, const Traces<n_loc>& loc_mat)
    {
      std::size_t index = e * block;
      for (unsigned int col = 0; col < n_loc; ++col)
        for (unsigned int row = 0; row < n_loc; ++row, ++index)
        {
          rows[index] = dofs[row];
          cols[index] = dofs[col];
          vals[index] = loc_mat(row, col) * free_mask_[dofs[row]] * free_mask_[dofs[col]];
        }
    };

    dispatch<Kind>(
      [&](auto mode)
      {
        sweep<Kind, mode.value>(args,
                                [&](std::size_t e, Edge& edge, const auto& op, const Dofs& dofs)
                                {
                                  Traces<n_loc> loc_mat;
                                  op.apply(local_solver_, edge, identity, args, loc_mat);
                                  write_block(e, dofs, loc_mat);
                                });
      });

    if (identity_rows)
    {
      std::size_t index = n_edges * block;
      for (dof_index_t i : constraints_)
      {
        rows[index] = i;
        cols[index] = i;
        vals[index++] = 1.;
      }
    }
  }
  /*!***********************************************************************************************
   * \brief   Fold f(edge, operator, local trace) over all hyperedges: result = combine(result,
   *          f(...)), by default with operator+.
   *
   * \c combine is associative and \c result its neutral element (0 for +), so that partial
   * results can be combined in any grouping.
   ************************************************************************************************/
  template <class Kind = Stiffness, class VecT, class F, class R, class Combine = std::plus<> >
  R reduce(const VecT& x, const Args& args, F f, R result, Combine combine = Combine())
  {
    dispatch<Kind>(
      [&](auto mode)
      {
        sweep<Kind, mode.value>(args,
                                [&](std::size_t, Edge& edge, const auto& op, const Dofs& dofs)
                                {
                                  Traces<1> in;
                                  gather(dofs, x, in);
                                  result = combine(std::move(result), f(edge, op, in));
                                });
      });
    return result;
  }

  // ---- plot
  // --------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Set an option of the plot (cf. PlotOptions in plot.hxx, e.g. "fileName", or
   *          "fileEnding" as "vtu" or "vtkhdf") and return its value; without a value, only that.
   ************************************************************************************************/
  std::string plot_option(const std::string& option, const std::string& value = "")
  {
    return set_plot_option(plot_options_, option, value);
  }
  /*!***********************************************************************************************
   * \brief   Write the local solutions of x to the file the plot options name (plot.hxx; the
   *          local solver's bulk_values). The format vtkhdf needs HyperHDG configured with HDF5.
   ************************************************************************************************/
  template <class VecT>
  void plot_solution(const VecT& x, const dof_value_t time = 0.)
  {
    plot(hyper_graph_, local_solver_, x, plot_options_, time);
  }

 private:
  using Dofs = std::array<dof_index_t, n_loc>;

  void collect_constraints()
  {
    constraints_.clear();
    for (dof_index_t i = 0; i < free_mask_.size(); ++i)
      if (free_mask_[i] == 0.)
        constraints_.push_back(i);
  }
  template <class VecInT, class VecOutT>
  void check_sizes(const VecInT& x, const VecOutT& y) const
  {
    hy_assert(x.size() == hyper_graph_.n_global_dofs() && y.size() == x.size(),
              "Vectors must have size_of_system() entries.");
  }
  template <class VecT>
  static void zero(VecT& y)
  {
    for (std::size_t i = 0; i < y.size(); ++i)
      y[i] = 0.;
  }
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
  // The local trace of x, zero on constrained dofs (branch-free mask).
  template <class VecT>
  void gather(const Dofs& dofs, const VecT& x, Traces<1>& loc) const
  {
    for (unsigned int i = 0; i < n_loc; ++i)
      loc[i] = x[dofs[i]] * free_mask_[dofs[i]];
  }
  template <class VecT>
  static void scatter(const Dofs& dofs, const Traces<1>& loc, VecT& y)
  {
    for (unsigned int i = 0; i < n_loc; ++i)
      y[dofs[i]] += loc[i];
  }
  /*!***********************************************************************************************
   * \brief   The store of a kind of operators.
   ************************************************************************************************/
  template <class Kind>
  Store<Kind>& store()
  {
    return std::get<Store<Kind> >(stores_);
  }
  template <class Kind>
  const Store<Kind>& store() const
  {
    return std::get<Store<Kind> >(stores_);
  }
  /*!***********************************************************************************************
   * \brief   Call f with std::true_type / std::false_type telling whether the store is used. The
   *          runtime decision is taken once per public call.
   ************************************************************************************************/
  template <class Kind, class F>
  void dispatch(F&& f) const
  {
    if (store<Kind>().prepared)
      f(std::true_type{});
    else
      f(std::false_type{});
  }
  /*!***********************************************************************************************
   * \brief   The number of slots and the slot of a hyperedge: a single one if operators have no
   *          state, else as the local solver numbers them, else one per hyperedge.
   ************************************************************************************************/
  template <class Kind>
  std::size_t n_operator_slots(const Args& args)
  {
    if constexpr (std::is_empty_v<OperatorOf<Kind> >)
      return 1;
    else if constexpr (requires(Edge& edge) { local_solver_.operator_slot(edge, Kind{}, args); })
    {
      if (hyper_graph_.n_hyEdges() == 0)
        return 0;
      Edge edge = hyper_graph_[0];  // any hyperedge can tell
      return local_solver_.n_operator_slots(edge, Kind{}, args);
    }
    else
      return hyper_graph_.n_hyEdges();
  }
  template <class Kind>
  std::size_t operator_slot(const std::size_t e, Edge& edge, const Args& args) const
  {
    if constexpr (std::is_empty_v<OperatorOf<Kind> >)
      return 0;
    else if constexpr (requires { local_solver_.operator_slot(edge, Kind{}, args); })
      return local_solver_.operator_slot(edge, Kind{}, args);
    else
      return e;
  }
  /*!***********************************************************************************************
   * \brief   Run body(e, edge, operator, dofs) over all hyperedges; the operator comes from the
   *          store or is built on the fly, decided at compile time.
   ************************************************************************************************/
  template <class Kind, bool use_store, class Body>
  void sweep(const Args& args, Body body)
  {
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    for (std::size_t e = 0; e < n_edges; ++e)
    {
      Edge edge = hyper_graph_[e];
      Dofs dofs;
      local_dof_indices(edge, dofs);
      if constexpr (use_store)
        body(e, edge, store<Kind>().operators[operator_slot<Kind>(e, edge, args)], dofs);
      else
      {
        const OperatorOf<Kind> op = local_solver_.local_operator(edge, Kind{}, args);
        body(e, edge, op, dofs);
      }
    }
  }
};  // end of class Generic

/*!*************************************************************************************************
 * \brief   Errors of a solution: the sum of the local solver's error_sq(edge, operator, lambda,
 *          args) over all hyperedges (a reduction of the loop), square root taken componentwise.
 **************************************************************************************************/
template <class LoopT, class VecT>
auto errors(LoopT& loop, const VecT& x, const typename LoopT::Args& args = typename LoopT::Args())
{
  using error_t = std::remove_cvref_t<decltype(loop.local_solver().error_sq(
    std::declval<typename LoopT::Edge&>(), std::declval<const typename LoopT::Operator&>(),
    std::declval<const typename LoopT::template Traces<1>&>(), args))>;
  error_t sum = loop.reduce(
    x, args, [&](auto& edge, const auto& op, const auto& lambda)
    { return loop.local_solver().error_sq(edge, op, lambda, args); }, error_t());
  for (unsigned int k = 0; k < sum.size(); ++k)
    sum[k] = std::sqrt(sum[k]);
  return sum;
}

}  // end of namespace GlobalLoop
