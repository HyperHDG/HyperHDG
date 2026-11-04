// What a sweep costs per hyperedge, by local solver and layout (local_operator.hxx), on the unit
// cube with DiffusionSin (constant coefficients: congruent elements share one operator):
// - Diffusion, Condensed (keeps S = C A^-1 B + D and K = C A^-1) and Factored (keeps the factors),
// - DiffusionTensor (sum factorization, fast diagonalization),
// - Diffusion with one operator per hyperedge (constant_coefficients = false, as for variable
//   coefficients: no sharing, no explicit elimination of the fluxes), i.e. the store's memory,
// and per variant: prepare, apply and residual on the store, apply streaming (no store: built per
// hyperedge and sweep), assemble (triplets of the condensed matrix). Degrees 1 to 3 in dimension
// DIM (set by CMakeLists.txt). Timings are per hyperedge.
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/local_solver/diffusion_tensor.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/parameters/diffusion.hxx>
#include <tests/performance/benchmark.hxx>

#include <random>
#include <string>
#include <vector>

#ifndef DIM
#define DIM 2  // set by CMakeLists.txt; the default serves editors and linters
#endif

// DiffusionSin, but declared variable: no operator_key (one operator per hyperedge), no
// elimination hint (the factorization is of the whole local matrix)
template <unsigned int space_dimT, typename param_float_t>
struct DiffusionSinUnshared : tests::DiffusionSin<space_dimT, param_float_t>
{
  static constexpr bool constant_coefficients = false;
};

template <class SolverT>
using Loop = GlobalLoop::Generic<Topology::Cubic<DIM, DIM>,
                                 Geometry::UnitCube<DIM, DIM, double>,
                                 NodeDescriptor::Cubic<DIM, DIM>,
                                 SolverT>;

template <class SolverT>
void measure(ankerl::nanobench::Bench& bench, const std::string& name, const unsigned int n)
{
  using ankerl::nanobench::doNotOptimizeAway;
  Loop<SolverT> loop(std::vector<unsigned int>(DIM, n));
  loop.constrain_dirichlet_nodes();
  const unsigned int size = loop.size_of_system();
  std::mt19937 generator(0);
  std::uniform_real_distribution<double> uniform(-1., 1.);
  std::vector<double> x(size), y(size);
  for (double& value : x)
    value = uniform(generator);

  bench.run(name + ": apply, streaming",
            [&]
            {
              loop.apply(x, y);
              doNotOptimizeAway(y.data());
            });
  bench.run(name + ": prepare", [&] { loop.prepare(); });
  bench.run(name + ": apply",
            [&]
            {
              loop.apply(x, y);
              doNotOptimizeAway(y.data());
            });
  bench.run(name + ": residual",
            [&]
            {
              loop.residual(x, y);
              doNotOptimizeAway(y.data());
            });
  std::vector<unsigned int> rows(loop.n_triplets()), cols(rows.size());
  std::vector<double> values(rows.size());
  bench.run(name + ": assemble",
            [&]
            {
              loop.assemble(typename Loop<SolverT>::Args(), false, rows, cols, values);
              doNotOptimizeAway(values.data());
            });
  std::printf("  %s: %zu operator(s) stored, %zu bytes each\n", name.c_str(), loop.n_operators(),
              sizeof(typename Loop<SolverT>::Operator));
}

template <unsigned int degree>
void measure_degree(const unsigned int n)
{
  using LocalOperator::Condensed;
  using LocalOperator::Factored;
  unsigned int n_edges = 1;
  for (unsigned int d = 0; d < DIM; ++d)
    n_edges *= n;
  ankerl::nanobench::Bench bench;
  const std::string title = std::to_string(DIM) + "D degree " + std::to_string(degree) + ", " +
                            std::to_string(n_edges) + " hyperedges";
  bench.title(title).unit("hyperedge").batch(n_edges).epochs(7).warmup(1).relative(false);
  measure<LocalSolver::Diffusion<DIM, degree, 2 * degree, tests::DiffusionSin, double, Condensed>>(
    bench, "Diffusion Condensed", n);
  measure<LocalSolver::Diffusion<DIM, degree, 2 * degree, tests::DiffusionSin, double, Factored>>(
    bench, "Diffusion Factored", n);
  measure<LocalSolver::DiffusionTensor<DIM, degree, 2 * degree, tests::DiffusionSin>>(
    bench, "DiffusionTensor", n);
  measure<LocalSolver::Diffusion<DIM, degree, 2 * degree, DiffusionSinUnshared, double, Condensed>>(
    bench, "Diffusion Condensed, unshared", n);
  tests::performance::save(
    bench, "local_operators_d" + std::to_string(DIM) + "_p" + std::to_string(degree));
}

int main()
{
  // about 1000 hyperedges in 2D, 512 in 3D: the unshared store of 3D degree 3 is 512 * 147 kB
  constexpr unsigned int n = DIM == 2 ? 32 : 8;
  measure_degree<1>(n);
  measure_degree<2>(n);
  measure_degree<3>(n);
}
