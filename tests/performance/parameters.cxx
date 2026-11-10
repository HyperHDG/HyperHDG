// What a parameter struct can take at runtime, and what it costs. The local solvers never create a
// parameter struct: they name its static members, and Diffusion passes its functions as template
// arguments to the quadrature (integrate_vol_phiphifunc<..., parameters::inverse_diffusion_coeff>),
// so the functions themselves must be static and fixed at compile time; so must every constant the
// solver branches on with if constexpr / requires (constant_coefficients, the array types). What a
// static function computes may still depend on runtime state, static data members of the struct:
// - Compiled: everything a compile-time constant (the reference),
// - Variable: a runtime coefficient (static inline double), set by the program,
// - Pointer: a runtime function pointer (static inline double (*)(...)),
// - Function: a runtime std::function (may capture state),
// - and dirichlet_nodes as a runtime std::vector (is_dirichlet uses only begin/end).
// Each for a cheap coefficient (1 + x_0^2, where a call dominates) on unshared operators
// (constant_coefficients = false: prepare evaluates the coefficient at the quadrature points of
// every hyperedge) and the residual (local_load evaluates the right-hand side). 2D degree 2.
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/performance/benchmark.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

// The shared part: u = sin(x_0) x_1, Dirichlet data everywhere, variable coefficient.
template <unsigned int space_dimT, typename param_float_t>
struct Common
{
  using P = Point<space_dimT, param_float_t>;
  static constexpr bool constant_coefficients = false;
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t analytic_result(const P& x, const param_float_t = 0.)
  {
    return std::sin(x[0]) * x[1];
  }
  static param_float_t dirichlet_value(const P& x, const param_float_t = 0.)
  {
    return analytic_result(x);
  }
  static param_float_t neumann_value(const P&, const param_float_t = 0.) { return 0.; }
  static param_float_t coefficient(const P& x, const param_float_t = 0.)
  {
    return 1. + x[0] * x[0];
  }
  static param_float_t load(const P& x, const param_float_t = 0.) { return x[0] + x[1]; }
};

template <unsigned int space_dimT, typename param_float_t>
struct Compiled : Common<space_dimT, param_float_t>
{
  using Base = Common<space_dimT, param_float_t>;
  static param_float_t inverse_diffusion_coeff(const typename Base::P& x,
                                               const param_float_t t = 0.)
  {
    return Base::coefficient(x, t);
  }
  static param_float_t right_hand_side(const typename Base::P& x, const param_float_t t = 0.)
  {
    return Base::load(x, t);
  }
};

template <unsigned int space_dimT, typename param_float_t>
struct Variable : Common<space_dimT, param_float_t>
{
  using Base = Common<space_dimT, param_float_t>;
  static inline param_float_t scale = 1.;  // set at runtime
  static param_float_t inverse_diffusion_coeff(const typename Base::P& x, const param_float_t = 0.)
  {
    return 1. + scale * x[0] * x[0];
  }
  static param_float_t right_hand_side(const typename Base::P& x, const param_float_t = 0.)
  {
    return scale * (x[0] + x[1]);
  }
};

template <unsigned int space_dimT, typename param_float_t>
struct Pointer : Common<space_dimT, param_float_t>
{
  using Base = Common<space_dimT, param_float_t>;
  using F = param_float_t (*)(const typename Base::P&, param_float_t);
  static inline F coefficient_pointer = nullptr, load_pointer = nullptr;  // set at runtime
  static param_float_t inverse_diffusion_coeff(const typename Base::P& x,
                                               const param_float_t t = 0.)
  {
    return coefficient_pointer(x, t);
  }
  static param_float_t right_hand_side(const typename Base::P& x, const param_float_t t = 0.)
  {
    return load_pointer(x, t);
  }
};

template <unsigned int space_dimT, typename param_float_t>
struct Function : Common<space_dimT, param_float_t>
{
  using Base = Common<space_dimT, param_float_t>;
  using F = std::function<param_float_t(const typename Base::P&, param_float_t)>;
  static inline F coefficient_function, load_function;  // set at runtime, may capture
  static param_float_t inverse_diffusion_coeff(const typename Base::P& x,
                                               const param_float_t t = 0.)
  {
    return coefficient_function(x, t);
  }
  static param_float_t right_hand_side(const typename Base::P& x, const param_float_t t = 0.)
  {
    return load_function(x, t);
  }
};

// Compiled functions, but the Dirichlet boundary is a runtime list
template <unsigned int space_dimT, typename param_float_t>
struct RuntimeBoundary : Compiled<space_dimT, param_float_t>
{
  static inline std::vector<unsigned int> dirichlet_nodes;
  static inline std::vector<unsigned int> neumann_nodes;
};

template <template <unsigned int, typename> typename ParametersT>
using Loop = GlobalLoop::Generic<Topology::Cubic<2, 2>,
                                 Geometry::UnitCube<2, 2, double>,
                                 NodeDescriptor::Cubic<2, 2>,
                                 LocalSolver::Diffusion<2, 2, 4, ParametersT>>;

int failures = 0;

template <template <unsigned int, typename> typename ParametersT>
std::vector<double> measure(ankerl::nanobench::Bench& bench, const std::string& name)
{
  using ankerl::nanobench::doNotOptimizeAway;
  Loop<ParametersT> loop(std::vector<unsigned int>{32, 32});
  loop.constrain_dirichlet_nodes();
  const unsigned int size = loop.size_of_system();
  std::vector<double> x(size), r(size);
  for (unsigned int i = 0; i < size; ++i)
    x[i] = std::sin(0.1 * i);
  bench.run(name + ": prepare (unshared)", [&] { loop.prepare(); });
  bench.run(name + ": residual",
            [&]
            {
              loop.residual(x, r);
              doNotOptimizeAway(r.data());
            });
  bench.run(name + ": apply, streaming",
            [&]
            {
              loop.release();
              loop.apply(x, r);
              doNotOptimizeAway(r.data());
            });
  loop.prepare();
  loop.residual(x, r);
  r.push_back(loop.constraints().size());
  return r;
}

double deviation(const std::vector<double>& a, const std::vector<double>& b)
{
  double d = a.size() == b.size() ? 0. : 1.;
  for (std::size_t i = 0; i < std::min(a.size(), b.size()); ++i)
    d = std::max(d, std::abs(a[i] - b[i]));
  return d;
}

int main()
{
  using P = Point<2, double>;
  Variable<2, double>::scale = 1.;
  Pointer<2, double>::coefficient_pointer = &Common<2, double>::coefficient;
  Pointer<2, double>::load_pointer = &Common<2, double>::load;
  const double one = 1.;  // captured state
  Function<2, double>::coefficient_function = [one](const P& x, double)
  { return one + x[0] * x[0]; };
  Function<2, double>::load_function = [one](const P& x, double) { return one * (x[0] + x[1]); };
  RuntimeBoundary<2, double>::dirichlet_nodes.assign(Common<2, double>::dirichlet_nodes.begin(),
                                                     Common<2, double>::dirichlet_nodes.end());

  ankerl::nanobench::Bench bench;
  bench.title("2D degree 2, 1024 hyperedges, one operator per hyperedge")
    .unit("hyperedge")
    .batch(1024)
    .epochs(11)
    .warmup(1)
    .relative(true);
  const std::vector<double> reference = measure<Compiled>(bench, "Compiled");
  const std::vector<std::vector<double>> results{
    measure<Variable>(bench, "Variable"), measure<Pointer>(bench, "Pointer"),
    measure<Function>(bench, "Function"), measure<RuntimeBoundary>(bench, "RuntimeBoundary")};
  tests::performance::save(bench, "parameters");

  // every variant computes the same residual (and the same constraints)
  const char* names[] = {"Variable", "Pointer", "Function", "RuntimeBoundary"};
  for (std::size_t k = 0; k < results.size(); ++k)
  {
    const double d = deviation(results[k], reference);
    std::printf("%s  %s: residual as Compiled (%.1e)\n", d <= 1e-12 ? "ok  " : "FAIL", names[k], d);
    failures += d > 1e-12;
  }
  // a runtime boundary changes at runtime: Dirichlet on nodes 1 and 2 only (x_0 = 0 and x_0 = 1)
  RuntimeBoundary<2, double>::dirichlet_nodes = {1, 2};
  Loop<RuntimeBoundary> loop(std::vector<unsigned int>{32, 32});
  loop.constrain_dirichlet_nodes();
  const bool fewer = loop.constraints().size() < reference.back();
  std::printf("%s  RuntimeBoundary: %zu instead of %.0f constrained dofs after the change\n",
              fewer ? "ok  " : "FAIL", loop.constraints().size(), reference.back());
  failures += !fewer;
  return failures;
}
