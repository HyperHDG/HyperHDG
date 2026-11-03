// LDG-H diffusion on the unit cube with GlobalLoop::Generic: solve u = sin(pi/2 x_0) on meshes of
// 2, 4, 8 elements per direction and check that the L2 error converges with order degree + 1 on
// the last refinement. On meshes small enough for a dense LU, the assembled matrix is also solved
// directly: it must be non-singular, and the direct solution must agree with the matrix-free CG.
// One executable per (DIM, DEGREE), set by CMake (tests/CMakeLists.txt).
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/wrapper/lapack.hxx>

#include <tests/parameters/diffusion.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

// Set per executable by tests/CMakeLists.txt; the defaults only serve editors and linters that
// compile this file without the CMake flags.
#ifndef DIM
#define DIM 2
#endif
#ifndef DEGREE
#define DEGREE 1
#endif

constexpr unsigned int dim = DIM, degree = DEGREE;
constexpr unsigned int max_direct_dofs = 2000;  // dense LU: 32 MB, well under a second

using Loop = GlobalLoop::Generic<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::Diffusion<dim, degree, 2 * degree, tests::DiffusionSin, double> >;

// Dense LU of the assembled matrix (constrained rows are identity rows, as in apply).
std::vector<double> direct_solve(Loop& loop, std::vector<double> rhs)
{
  const unsigned int size = loop.size_of_system();
  const std::size_t n_triplets = loop.n_triplets(true);
  std::vector<unsigned int> rows(n_triplets), cols(n_triplets);
  std::vector<double> vals(n_triplets), dense((std::size_t)size * size, 0.);
  loop.assemble(Loop::Args(), true, rows, cols, vals);
  for (std::size_t k = 0; k < n_triplets; ++k)
    dense[(std::size_t)cols[k] * size + rows[k]] += vals[k];  // column-major, duplicates summed
  Wrapper::lapack_solve(size, 1, dense.data(), rhs.data());   // throws if the matrix is singular
  return rhs;
}

int main()
{
  int failures = 0;
  double previous = 0.;
  for (unsigned int n = 2; n <= 8; n *= 2)
  {
    Loop loop(std::vector<unsigned int>(dim, n));
    loop.constrain_dirichlet_nodes();
    loop.prepare();

    const unsigned int size = loop.size_of_system();
    std::vector<double> rhs(size), zero(size, 0.);
    loop.residual(zero, rhs);
    for (double& value : rhs)
      value = -value;

    SparseLA::SolverOptions options;
    options.rtol = 1e-12;
    options.max_iterations = 10 * size;
    std::vector<double> x(size, 0.);
    const SparseLA::SolverResult cg =
      SparseLA::cg([&](const std::vector<double>& in, std::vector<double>& out)
                   { loop.apply(in, out); }, rhs, x, options);
    const double error = loop.errors(x)[0];
    std::printf("- {dim: %u, degree: %u, n: %u, dofs: %u, cg: %s, iterations: %u, error: %.6e}\n",
                dim, degree, n, size, SparseLA::to_string(cg.reason), cg.iterations, error);
    if (!cg.converged())
    {
      std::printf("  FAIL: CG did not converge\n");
      ++failures;
    }
    if (!std::isfinite(error))
    {
      std::printf("  FAIL: the error is not finite\n");
      ++failures;
    }

    if (size <= max_direct_dofs)
    {
      const std::vector<double> x_direct = direct_solve(loop, rhs);
      double difference = 0., scale = 0.;
      for (unsigned int i = 0; i < size; ++i)
      {
        difference = std::max(difference, std::abs(x_direct[i] - x[i]));
        scale = std::max(scale, std::abs(x_direct[i]));
      }
      std::printf("  direct vs cg: %.2e (relative, at most 1e-8)\n", difference / scale);
      if (!(difference <= 1e-8 * scale))
        ++failures;
    }

    // Only the last rate is checked: the coarse meshes are preasymptotic.
    if (n == 8)
    {
      const double rate = std::log2(previous / error);
      std::printf("  rate: %.3f (expected at least %.1f)\n", rate, degree + 1 - 0.3);
      if (!(rate >= degree + 1 - 0.3))
        ++failures;
    }
    previous = error;
  }
  return failures;
}
