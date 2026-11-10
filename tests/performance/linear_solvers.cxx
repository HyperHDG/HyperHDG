// CG against a sparse direct solver (CHOLMOD) on the condensed system of LDG-H diffusion
// (DiffusionSin, Condensed layout), unpreconditioned CG to rtol 1e-10 from zero:
// - matrix-free CG: prepare + CG on Generic::apply,
// - assembled CG: assemble (triplets, compressed by CHOLMOD) + CG on the product with the matrix,
// - CHOLMOD (SparseLA::Cholmod): assemble + analyze (ordering, symbolic) + factorize (numeric) +
//   solve per right-hand side.
// The matrix is that of all dofs, with a diagonal entry on the constrained ones; the direct
// solution is checked against that of CG.
// The direct solver pays once (assemble, analyze, factorize) and then little per right-hand side;
// CG pays per right-hand side, with iterations growing like 1/h. Break-even: the number of
// right-hand sides with the same operator (e.g. time steps of a linear problem with a fixed step)
// after which CHOLMOD is cheaper than matrix-free CG, (setup_direct - prepare) / (cg - solve). The
// fill-in, and hence factorize, grows much faster in 3D than in 2D.
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/parameters/diffusion.hxx>
#include <tests/performance/benchmark.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <optional>
#include <string>
#include <vector>

using Vec = std::vector<double>;
using ankerl::nanobench::doNotOptimizeAway;
using ankerl::nanobench::Result;

struct Row
{
  std::string name;
  std::size_t dofs, nnz;
  double factor_nnz;
  unsigned int iterations;
  double prepare, cg, assemble, csr_cg, analyze, factorize, solve, difference;
};
std::vector<Row> rows;

double last_median(const ankerl::nanobench::Bench& bench)
{
  return bench.results().back().median(Result::Measure::elapsed);
}

template <unsigned int dim, unsigned int degree>
void measure(ankerl::nanobench::Bench& bench, const unsigned int n)
{
  using Loop =
    GlobalLoop::Generic<Topology::Cubic<dim, dim>, Geometry::UnitCube<dim, dim, double>,
                        NodeDescriptor::Cubic<dim, dim>,
                        LocalSolver::Diffusion<dim, degree, 2 * degree, tests::DiffusionSin>>;
  const std::string name =
    std::to_string(dim) + "D p" + std::to_string(degree) + " n" + std::to_string(n);
  Row row{name};
  Loop loop(std::vector<unsigned int>(dim, n));
  loop.constrain_dirichlet_nodes();
  const unsigned int size = loop.size_of_system();
  SparseLA::SolverOptions options;
  options.rtol = 1e-10;
  options.max_iterations = 100000;

  // matrix-free
  bench.run(name + ": prepare", [&] { loop.prepare(); });
  row.prepare = last_median(bench);
  Vec rhs(size), zero(size, 0.), x(size);
  loop.residual(zero, rhs);
  for (double& value : rhs)
    value = -value;
  bench.run(name + ": matrix-free CG",
            [&]
            {
              std::fill(x.begin(), x.end(), 0.);
              row.iterations =
                SparseLA::cg([&](const Vec& in, Vec& out) { loop.apply(in, out); }, rhs, x, options)
                  .iterations;
              doNotOptimizeAway(x.data());
            });
  row.cg = last_median(bench);

  // assembled: the condensed matrix is negative definite, so the constrained dofs get -1 on the
  // diagonal (assemble appends their entries last) and CHOLMOD factorizes -A
  const std::size_t n_triplets = loop.n_triplets(true);
  std::vector<unsigned int> triplet_rows(n_triplets), triplet_cols(n_triplets);
  Vec triplet_vals(n_triplets);
  std::optional<SparseLA::Cholmod> cholmod;
  bench.run(name + ": assemble",
            [&]
            {
              loop.assemble(typename Loop::Args(), true, triplet_rows, triplet_cols, triplet_vals);
              std::fill(triplet_vals.begin() + loop.n_triplets(false), triplet_vals.end(), -1.);
              cholmod.emplace(size, triplet_rows, triplet_cols, triplet_vals, -1.);
            });
  row.assemble = last_median(bench);
  row.dofs = size;
  row.nnz = cholmod->nnz();
  Vec y(size);
  bench.run(name + ": assembled CG",
            [&]
            {
              std::fill(y.begin(), y.end(), 0.);
              SparseLA::cg([&](const Vec& in, Vec& out) { cholmod->multiply(in, out); }, rhs, y,
                           options);
              doNotOptimizeAway(y.data());
            });
  row.csr_cg = last_median(bench);

  bench.run(name + ": CHOLMOD analyze", [&] { cholmod->analyze(); });
  row.analyze = last_median(bench);
  bench.run(name + ": CHOLMOD factorize", [&] { cholmod->factorize(); });
  row.factorize = last_median(bench);
  bench.run(name + ": CHOLMOD solve",
            [&]
            {
              cholmod->solve(rhs, y);
              doNotOptimizeAway(y.data());
            });
  row.solve = last_median(bench);
  row.factor_nnz = cholmod->factor_nnz();
  double difference = 0., scale = 0.;
  for (unsigned int i = 0; i < size; ++i)
  {
    difference = std::max(difference, std::abs(y[i] - x[i]));
    scale = std::max(scale, std::abs(x[i]));
  }
  row.difference = difference / scale;
  rows.push_back(row);
}

int main()
{
  ankerl::nanobench::Bench bench;
  bench.title("CG vs CHOLMOD").unit("solve").epochs(3).epochIterations(1).warmup(0);
  measure<2, 1>(bench, 16);
  measure<2, 1>(bench, 64);
  measure<2, 1>(bench, 256);
  measure<2, 3>(bench, 16);
  measure<2, 3>(bench, 64);
  measure<3, 1>(bench, 4);
  measure<3, 1>(bench, 8);
  measure<3, 1>(bench, 16);
  measure<3, 2>(bench, 4);
  measure<3, 2>(bench, 8);
  measure<3, 2>(bench, 12);
  tests::performance::save(bench, "linear_solvers");

  std::printf("\n%-12s %8s %9s %10s %6s | %9s %9s | %9s %9s | %9s %9s %9s | %8s | %s\n", "case",
              "dofs", "nnz", "factor nnz", "its", "prepare", "MF CG", "assemble", "asm. CG",
              "analyze", "factorize", "solve", "vs CG", "break-even (rhs)");
  for (const Row& r : rows)
  {
    const double setup = r.assemble + r.analyze + r.factorize - r.prepare, saving = r.cg - r.solve;
    std::printf(
      "%-12s %8zu %9zu %10.0f %6u | %9.2e %9.2e | %9.2e %9.2e | %9.2e %9.2e %9.2e | %8.1e | ",
      r.name.c_str(), r.dofs, r.nnz, r.factor_nnz, r.iterations, r.prepare, r.cg, r.assemble,
      r.csr_cg, r.analyze, r.factorize, r.solve, r.difference);
    if (setup <= saving)
      std::printf("1 (direct wins from the first)\n");
    else if (saving > 0.)
      std::printf("%.1f\n", setup / saving);
    else
      std::printf("never\n");
  }
}
