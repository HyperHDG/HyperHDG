// The interface of GlobalLoop::Generic on one small instantiation (2D, degree 1, 3 x 3 elements):
// properties a solve does not check by itself, such as that the sweeps without prepare agree with
// the prepared ones, that the assembled matrix is the matrix-free operator in both modes, and how
// constraints, the store and reductions behave. Run for every layout of the local operator
// (local_solver/local_operator.hxx): Condensed (the default), Factored, Streaming (nothing
// stored), and a layout of one's own defined here (SinglePrecision); the layouts must agree with
// each other.
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/parameters/diffusion.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

// A layout of one's own: the condensed operator stored in single precision (half the memory of
// Condensed) and applied in double. A layout is any type with make(solver, edge, kind, args)
// returning an operator with apply and residual; this one has no local_solution, so errors
// factorize on the fly.
struct SinglePrecision
{
  template <class CondensedT>
  struct Operator
  {
    static constexpr unsigned int n_ldofs = CondensedT::n_ldofs, n_loc = CondensedT::n_loc;
    using T = typename CondensedT::T;

    SmallSquareMat<n_loc, float> S;
    SmallMat<n_loc, n_ldofs, float> K;

    template <unsigned int k, class SolverT, class Edge, class Args>
    void apply(const SolverT&,
               Edge&,
               const SmallMat<n_loc, k, T>& in,
               const Args&,
               SmallMat<n_loc, k, T>& out) const
    {
      for (unsigned int col = 0; col < k; ++col)
        for (unsigned int row = 0; row < n_loc; ++row)
        {
          T sum = 0.;
          for (unsigned int j = 0; j < n_loc; ++j)
            sum += (T)S(row, j) * in(j, col);
          out(row, col) = sum;
        }
    }
    template <class SolverT, class Edge, class Args>
    void residual(const SolverT& solver,
                  Edge& edge,
                  const SmallVec<n_loc, T>& lambda,
                  const Args& args,
                  SmallVec<n_loc, T>& r) const
    {
      apply(solver, edge, lambda, args, r);
      const SmallVec<n_ldofs, T> load = solver.local_load(edge, args);
      for (unsigned int row = 0; row < n_loc; ++row)
        for (unsigned int j = 0; j < n_ldofs; ++j)
          r[row] += (T)K(row, j) * load[j];
    }
  };
  template <class SolverT, class Edge, class Kind, class Args>
  static auto make(const SolverT& solver, Edge& edge, Kind kind, const Args& args)
  {
    const auto condensed = LocalOperator::Condensed::make(solver, edge, kind, args);
    Operator<std::remove_cvref_t<decltype(condensed)> > op;
    for (unsigned int i = 0; i < op.n_loc * op.n_loc; ++i)
      op.S[i] = (float)condensed.S[i];
    for (unsigned int i = 0; i < op.n_loc * op.n_ldofs; ++i)
      op.K[i] = (float)condensed.K[i];
    return op;
  }
};

// The same problem with the coefficients declared as varying: the local matrix is factorized
// whole (no fluxes eliminated) and every hyperedge has an operator of its own (no operator_slot).
template <unsigned int space_dimT, typename param_float_t = double>
struct DiffusionSinVarying : tests::DiffusionSin<space_dimT, param_float_t>
{
  static constexpr bool constant_coefficients = false;
};

template <class LayoutT, template <unsigned int, typename> class ParametersT = tests::DiffusionSin>
using LoopT = GlobalLoop::Generic<Topology::Cubic<2, 2>,
                                  Geometry::UnitCube<2, 2, double>,
                                  NodeDescriptor::Cubic<2, 2>,
                                  LocalSolver::Diffusion<2, 1, 2, ParametersT, double, LayoutT> >;
using Vec = std::vector<double>;

int failures = 0;
void check(const bool ok, const char* what)
{
  std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what);
  failures += !ok;
}
double max_difference(const Vec& a, const Vec& b)
{
  double difference = 0.;
  for (std::size_t i = 0; i < a.size(); ++i)
    difference = std::max(difference, std::abs(a[i] - b[i]));
  return difference;
}
template <class Loop>
Vec product(Loop& loop, const Vec& x)
{
  Vec y(x.size());
  loop.apply(x, y);
  return y;
}
template <class Loop>
Vec residual_of(Loop& loop, const Vec& x)
{
  Vec r(x.size());
  loop.residual(x, r);
  return r;
}
// Dense column-major matrix from the triplets of assemble; duplicate entries are summed.
template <class Loop>
Vec assembled(Loop& loop, const bool identity_rows)
{
  const std::size_t n = loop.size_of_system(), m = loop.n_triplets(identity_rows);
  std::vector<unsigned int> rows(m), cols(m);
  Vec vals(m), dense(n * n, 0.);
  loop.assemble(typename Loop::Args(), identity_rows, rows, cols, vals);
  for (std::size_t k = 0; k < m; ++k)
    dense[cols[k] * n + rows[k]] += vals[k];
  return dense;
}
Vec multiply(const Vec& dense, const Vec& x)
{
  Vec y(x.size(), 0.);
  for (std::size_t col = 0; col < x.size(); ++col)
    for (std::size_t row = 0; row < x.size(); ++row)
      y[row] += dense[col * x.size() + row] * x[col];
  return y;
}

// What every layout must agree on: the prepared apply, residual, matrix and error at x.
struct Results
{
  Vec y, r, a;
  double error;
};

template <class LayoutT>
Results test(const char* layout, const double tol)
{
  std::printf("---- %s\n", layout);
  using Loop = LoopT<LayoutT>;
  Loop loop(std::vector<unsigned int>{3, 3});
  const unsigned int n = loop.size_of_system();
  check(loop.n_hyEdges() == 9, "9 hyperedges");
  check(n == 24 * 2, "24 faces with 2 dofs each");

  std::mt19937 generator(1);
  std::uniform_real_distribution<double> uniform(-1., 1.);
  Vec x(n), z(n), zero(n, 0.);
  for (unsigned int i = 0; i < n; ++i)
    x[i] = uniform(generator), z[i] = uniform(generator);

  // ---- constraints
  loop.constrain_dirichlet_nodes();
  const auto constraints = loop.constraints();
  check(constraints.size() == 12 * 2, "Dirichlet constraints: the 12 boundary faces");
  check(std::is_sorted(constraints.begin(), constraints.end()), "constraints are sorted");
  loop.set_constraints(std::vector<unsigned int>(constraints.rbegin(), constraints.rend()));
  check(loop.constraints() == constraints, "set_constraints reproduces them, sorted");
  std::vector<bool> constrained(n, false);
  for (unsigned int i : constraints)
    constrained[i] = true;

  // ---- without prepare (an operator per hyperedge and sweep; named stream below)
  check(!loop.prepared() && loop.n_operators() == 0, "not prepared after construction");
  const Vec y_stream = product(loop, x), r_stream = residual_of(loop, x),
            r0 = residual_of(loop, zero);
  bool identity = true, affine = true;
  for (unsigned int i = 0; i < n; ++i)
    if (constrained[i])
      identity = identity && y_stream[i] == x[i] && r_stream[i] == x[i];
    else
      affine = affine && std::abs(r_stream[i] - r0[i] - y_stream[i]) <= tol;
  check(identity, "apply and residual are the identity on constrained dofs");
  check(affine, "residual(x) - residual(0) = apply(x) on free dofs");
  {
    Vec combination(n);
    for (unsigned int i = 0; i < n; ++i)
      combination[i] = 2. * x[i] - z[i];
    const Vec y_z = product(loop, z), y_combination = product(loop, combination);
    Vec expected(n);
    for (unsigned int i = 0; i < n; ++i)
      expected[i] = 2. * y_stream[i] - y_z[i];
    check(max_difference(y_combination, expected) <= tol, "apply is linear");
  }
  const Vec a_stream = assembled(loop, true);
  check(max_difference(multiply(a_stream, x), y_stream) <= tol,
        "without prepare: assembled matrix (identity rows) times x = apply(x)");
  {
    const Vec a_masked = assembled(loop, false);
    bool zero_rows = true;
    for (unsigned int col = 0; col < n; ++col)
      for (unsigned int i : constraints)
        zero_rows = zero_rows && a_masked[col * n + i] == 0.;
    check(zero_rows, "without identity rows, constrained rows are zero");
    check(loop.n_triplets(true) == loop.n_triplets(false) + constraints.size(),
          "n_triplets(true) adds one entry per constraint");
  }
  bool symmetric = true;
  for (unsigned int row = 0; row < n; ++row)
    for (unsigned int col = 0; col < n; ++col)
      if (!constrained[row] && !constrained[col])
        symmetric = symmetric && std::abs(a_stream[col * n + row] - a_stream[row * n + col]) <= tol;
  check(symmetric, "the free-free block is symmetric");

  // ---- prepared store
  loop.prepare();
  check(loop.prepared(), "prepared after prepare");
  check(loop.n_operators() == 1, "congruent hyperedges share one operator (operator_slot)");
  check(max_difference(product(loop, x), y_stream) <= tol, "prepared apply = apply without");
  check(max_difference(residual_of(loop, x), r_stream) <= tol,
        "prepared residual = residual without");
  check(max_difference(assembled(loop, true), a_stream) <= tol,
        "prepared assemble = assemble without");
  const double error_prepared = GlobalLoop::errors(loop, x)[0];
  const Results results{product(loop, x), residual_of(loop, x), assembled(loop, true),
                        error_prepared};

  // ---- reductions, release
  const double count =
    loop.reduce(x, typename Loop::Args(), [](auto&, const auto&, const auto&) { return 1.; }, 0.);
  check(count == 9., "reduce visits every hyperedge once");
  const double largest = loop.reduce(
    x, typename Loop::Args(),
    [](auto&, const auto&, const auto& lambda)
    {
      double result = 0.;
      for (unsigned int i = 0; i < lambda.size(); ++i)
        result = std::max(result, std::abs(lambda[i]));
      return result;
    },
    0., [](const double a, const double b) { return std::max(a, b); });
  double largest_x = 0.;
  for (const double value : x)
    largest_x = std::max(largest_x, std::abs(value));
  check(largest == largest_x, "reduce with max as the combination");
  loop.release();
  check(!loop.prepared() && loop.n_operators() == 0, "release empties the store");
  check(std::abs(GlobalLoop::errors(loop, x)[0] - error_prepared) <= tol * error_prepared,
        "errors: without prepare = prepared");
  return results;
}

// The Results of test() on the same x, without the checks (for the varying parameters).
template <class Loop>
Results results_only()
{
  Loop loop(std::vector<unsigned int>{3, 3});
  const unsigned int n = loop.size_of_system();
  std::mt19937 generator(1);
  std::uniform_real_distribution<double> uniform(-1., 1.);
  Vec x(n), z(n);
  for (unsigned int i = 0; i < n; ++i)
    x[i] = uniform(generator), z[i] = uniform(generator);
  loop.constrain_dirichlet_nodes();
  loop.prepare();
  check(loop.n_operators() == loop.n_hyEdges(), "varying coefficients: an operator per hyperedge");
  return Results{product(loop, x), residual_of(loop, x), assembled(loop, true),
                 GlobalLoop::errors(loop, x)[0]};
}

int main()
{
  const Results condensed = test<LocalOperator::Condensed>("Condensed (S, K)", 1e-12);
  const Results factored = test<LocalOperator::Factored>("Factored (LU, B, C, D)", 1e-12);
  const Results streaming = test<LocalOperator::Streaming>("Streaming (nothing stored)", 1e-12);
  const Results single = test<SinglePrecision>("SinglePrecision (a layout of one's own)", 1e-5);

  std::printf("---- layouts against each other\n");
  auto agree = [](const Results& a, const Results& b, const double tol)
  {
    auto relative = [](const Vec& u, const Vec& v)
    {
      double scale = 0.;
      for (const double value : v)
        scale = std::max(scale, std::abs(value));
      return max_difference(u, v) / scale;
    };
    return relative(a.y, b.y) <= tol && relative(a.r, b.r) <= tol && relative(a.a, b.a) <= tol &&
           std::abs(a.error - b.error) <= tol * b.error;
  };
  check(agree(condensed, factored, 1e-12), "Condensed = Factored: apply, residual, matrix, error");
  check(agree(streaming, factored, 1e-12), "Streaming = Factored");
  check(sizeof(LoopT<LocalOperator::Streaming>::Operator) == 1, "the Streaming operator is empty");
  // with the fluxes eliminated (constant coefficients) and without (declared varying)
  const Results whole = results_only<LoopT<LocalOperator::Factored, DiffusionSinVarying> >();
  const Results whole_condensed =
    results_only<LoopT<LocalOperator::Condensed, DiffusionSinVarying> >();
  check(agree(factored, whole, 1e-12), "fluxes eliminated = whole local matrix factorized");
  check(agree(whole_condensed, whole, 1e-12), "and Condensed = Factored without elimination");
  check(agree(single, condensed, 1e-6), "SinglePrecision = Condensed up to float rounding");
  return failures;
}
