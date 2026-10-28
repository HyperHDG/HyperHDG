#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/global_loop/kinds.hxx>
#include <HyperHDG/global_loop/prototype.hxx>  // sparse_mat
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <map>
#include <type_traits>
#include <utility>
#include <vector>

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Prototype of the single global loop built on local operators.
 *
 * The loop is a scatter--gather engine over hyperedges. For every hyperedge it gathers the local
 * trace from a global vector, hands it to the hyperedge's local operator, and scatters the result.
 * It never sees block matrices, factorizations, or boundary conditions of the local solver, and
 * it never allocates global vectors: every entry point writes into caller-provided vectors (any
 * type with size() and operator[], e.g. std::vector or std::span).
 *
 * Local solver protocol (cf. local_operator_design_v2.md):
 * - \c Args: runtime parameters (time, ...) forwarded to every call,
 * - \c local_operator(edge, Stiffness, args): builds the immutable per-edge operator with
 *   \c apply<k>(edge, in, args, out) and \c residual(edge, lambda, load, args, r),
 * - \c local_load(edge, args): the local right-hand side (loads, Dirichlet data) used by residual,
 * - \c error_sq(edge, op, lambda, args): per-edge error contribution; its type needs operator+,
 * - optional \c operator_key(edge, Stiffness, args): share one operator among equal edges,
 * - optional \c constrained_node<space_dim>(node_type): default predicate for constrain_dirichlet.
 *
 * Store: after \c prepare the operators live in a contiguous vector (one per distinct key) until
 * \c release. Without prepare every sweep streams: build, use, drop per hyperedge. Every public
 * sweep takes that decision once and runs a branch-free per-edge body for the chosen mode.
 *
 * Constraints: a sorted index set of global dofs owned by the loop. apply is the identity there,
 * residual returns x_i (prescribed values are zero in this prototype, Dirichlet data enters via
 * the local load), assemble replaces the rows by zero or identity rows.
 *
 * Serial. Parallel sweeps (scatter strategy) are a separate step.
 **************************************************************************************************/
template <class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename dof_value_t = double,
          typename dof_index_t = unsigned int>
class LocalOperatorLoop
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
  using Operator = std::remove_cvref_t<
    decltype(std::declval<const LocalSolverT&>()
               .local_operator(std::declval<Edge&>(), Stiffness{}, std::declval<const Args&>()))>;
  template <unsigned int k>
  using Traces = SmallMat<n_loc, k, dof_value_t>;
  using error_t = std::remove_cvref_t<decltype(std::declval<const LocalSolverT&>().error_sq(
    std::declval<Edge&>(),
    std::declval<const Operator&>(),
    std::declval<const Traces<1>&>(),
    std::declval<const Args&>()))>;
  static_assert(
    requires(error_t a, error_t b) { a + b; },
    "The error type of the local solver must provide operator+ (e.g. a SmallVec).");

 private:
  HyperGraphT hyper_graph_;
  const LocalSolverT local_solver_;

  std::vector<Operator> store_;              // one operator per distinct key (prepared only)
  std::vector<std::size_t> slot_;            // hyperedge -> index into store_
  std::vector<std::size_t> representative_;  // slot -> a hyperedge with that key
  bool prepared_ = false;

  std::vector<dof_index_t> constraints_;  // sorted global dof indices
  std::vector<dof_value_t> free_mask_;    // 1 on free dofs, 0 on constrained dofs

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

  // ---- sizes, refinement
  // ---------------------------------------------------------------------------

  dof_index_t size_of_system() const { return hyper_graph_.n_global_dofs(); }
  std::size_t n_hyEdges() const { return hyper_graph_.n_hyEdges(); }
  unsigned int get_refinement() { return hyper_graph_.get_refinement(); }
  void set_refinement(unsigned int level)
  {
    hyper_graph_.set_refinement(level);
    release();
    constraints_.clear();
    free_mask_.clear();
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
   * \brief   Build the operators of all hyperedges (one per distinct key) and keep them.
   ************************************************************************************************/
  void prepare(const Args& args = Args())
  {
    release();
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    slot_.resize(n_edges);

    if constexpr (requires(Edge& e) { local_solver_.operator_key(e, Stiffness{}, args); })
    {
      using Key = std::remove_cvref_t<decltype(local_solver_.operator_key(std::declval<Edge&>(),
                                                                          Stiffness{}, args))>;
      std::map<Key, std::size_t> slot_of_key;
      for (std::size_t e = 0; e < n_edges; ++e)
      {
        Edge edge = hyper_graph_[e];
        auto [it, inserted] = slot_of_key.emplace(
          local_solver_.operator_key(edge, Stiffness{}, args), representative_.size());
        if (inserted)
          representative_.push_back(e);
        slot_[e] = it->second;
      }
    }
    else
    {
      representative_.resize(n_edges);
      for (std::size_t e = 0; e < n_edges; ++e)
        slot_[e] = representative_[e] = e;
    }

    store_.resize(representative_.size());
    try
    {
      for (std::size_t s = 0; s < representative_.size(); ++s)
      {
        Edge edge = hyper_graph_[representative_[s]];
        store_[s] = local_solver_.local_operator(edge, Stiffness{}, args);
      }
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

  // ---- sweeps
  // ------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   y = J x with identity rows on constrained dofs. y is overwritten.
   ************************************************************************************************/
  template <class VecInT, class VecOutT>
  void apply(const VecInT& x, VecOutT& y, const Args& args = Args())
  {
    check_sizes(x, y);
    zero(y);
    dispatch(
      [&](auto mode)
      {
        sweep<mode.value>(args,
                          [&](std::size_t, Edge& edge, const Operator& op, const Dofs& dofs)
                          {
                            Traces<1> in, res;
                            gather(dofs, x, in);
                            op.apply(edge, in, args, res);
                            scatter(dofs, res, y);
                          });
      });
    for (dof_index_t i : constraints_)
      y[i] = x[i];
  }
  /*!***********************************************************************************************
   * \brief   r = F(x): flux balance on free dofs, x_i on constrained dofs. r is overwritten.
   ************************************************************************************************/
  template <class VecInT, class VecOutT>
  void residual(const VecInT& x, VecOutT& r, const Args& args = Args())
  {
    check_sizes(x, r);
    zero(r);
    dispatch(
      [&](auto mode)
      {
        sweep<mode.value>(args,
                          [&](std::size_t, Edge& edge, const Operator& op, const Dofs& dofs)
                          {
                            Traces<1> in, res;
                            gather(dofs, x, in);
                            op.residual(edge, in, local_solver_.local_load(edge, args), args, res);
                            scatter(dofs, res, r);
                          });
      });
    for (dof_index_t i : constraints_)
      r[i] = x[i];
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
   * operator by construction). Entries in constrained rows are multiplied by zero (branch-free
   * mask); with \c identity_rows a unit diagonal entry is appended for every constrained dof.
   * Duplicate (row, col) pairs have to be summed by the consumer (as scipy's coo_matrix does).
   ************************************************************************************************/
  template <class IndexVecT, class ValueVecT>
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
          rows[index] = dofs[row];
          cols[index] = dofs[col];
          vals[index] = loc_mat(row, col) * free_mask_[dofs[row]];
        }
    };

    if (prepared_ && store_.size() < n_edges)
    {
      // Shared operators: one dense block per slot, then a copy per hyperedge.
      std::vector<Traces<n_loc> > blocks(store_.size());
      for (std::size_t s = 0; s < store_.size(); ++s)
      {
        Edge edge = hyper_graph_[representative_[s]];
        store_[s].apply(edge, identity, args, blocks[s]);
      }
      for (std::size_t e = 0; e < n_edges; ++e)
      {
        Dofs dofs;
        local_dof_indices(hyper_graph_[e], dofs);
        write_block(e, dofs, blocks[slot_[e]]);
      }
    }
    else
      dispatch(
        [&](auto mode)
        {
          sweep<mode.value>(args,
                            [&](std::size_t e, Edge& edge, const Operator& op, const Dofs& dofs)
                            {
                              Traces<n_loc> loc_mat;
                              op.apply(edge, identity, args, loc_mat);
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
   * \brief   Fold f(edge, operator, local trace) over all hyperedges with operator+.
   ************************************************************************************************/
  template <class VecT, class F, class R>
  R reduce(const VecT& x, const Args& args, F f, R result)
  {
    dispatch(
      [&](auto mode)
      {
        sweep<mode.value>(args,
                          [&](std::size_t, Edge& edge, const Operator& op, const Dofs& dofs)
                          {
                            Traces<1> in;
                            gather(dofs, x, in);
                            result = result + f(edge, op, in);
                          });
      });
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Errors: sum of the local solver's error_sq, square root taken componentwise.
   ************************************************************************************************/
  template <class VecT>
  error_t errors(const VecT& x, const Args& args = Args())
  {
    error_t sum = reduce(
      x, args, [&](Edge& edge, const Operator& op, const Traces<1>& lam)
      { return local_solver_.error_sq(edge, op, lam, args); }, error_t());
    for (unsigned int k = 0; k < sum.size(); ++k)
      sum[k] = std::sqrt(sum[k]);
    return sum;
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
  template <class VecT>
  static void gather(const Dofs& dofs, const VecT& x, Traces<1>& loc)
  {
    for (unsigned int i = 0; i < n_loc; ++i)
      loc[i] = x[dofs[i]];
  }
  template <class VecT>
  static void scatter(const Dofs& dofs, const Traces<1>& loc, VecT& y)
  {
    for (unsigned int i = 0; i < n_loc; ++i)
      y[dofs[i]] += loc[i];
  }
  /*!***********************************************************************************************
   * \brief   Call f with std::true_type / std::false_type telling whether the store is used. The
   *          runtime decision is taken once per public call.
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
   * \brief   Run body(e, edge, operator, dofs) over all hyperedges; the operator comes from the
   *store or is built on the fly, decided at compile time.
   ************************************************************************************************/
  template <bool use_store, class Body>
  void sweep(const Args& args, Body body)
  {
    const std::size_t n_edges = hyper_graph_.n_hyEdges();
    for (std::size_t e = 0; e < n_edges; ++e)
    {
      Edge edge = hyper_graph_[e];
      Dofs dofs;
      local_dof_indices(edge, dofs);
      if constexpr (use_store)
        body(e, edge, store_[slot_[e]], dofs);
      else
      {
        const Operator op = local_solver_.local_operator(edge, Stiffness{}, args);
        body(e, edge, op, dofs);
      }
    }
  }
};  // end of class LocalOperatorLoop

}  // end of namespace GlobalLoop
