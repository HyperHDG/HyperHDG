// Timoshenko beams with GlobalLoop::Generic on the cross of domains/cross.geo (run in the source
// tree) with the manufactured solution of tests/parameters/timoshenko.hxx: solve on 2, 4, 8
// subintervals per beam and check that the L2 error of the displacement converges with order
// degree + 1 on the last refinement. The assembled matrix is also solved directly: it must be
// non-singular, and the direct solution must agree with the matrix-free CG. One executable per
// DEGREE, set by CMake (tests/CMakeLists.txt). The records are checked against
// tests/expected/timoshenko_p<DEGREE>.yaml (tests/records.hxx).
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/file.hxx>

#include <tests/direct_solve.hxx>
#include <tests/parameters/timoshenko.hxx>
#include <tests/records.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

// Set per executable by tests/CMakeLists.txt; the default only serves editors and linters that
// compile this file without the CMake flags.
#ifndef DEGREE
#define DEGREE 1
#endif

constexpr unsigned int degree = DEGREE;
constexpr unsigned int max_direct_dofs = 2000;  // dense LU: 32 MB, well under a second

using Loop = GlobalLoop::Generic<
  Topology::File<1, 3>,
  Geometry::File<1, 3>,
  NodeDescriptor::File<1, 3>,
  LocalSolver::TimoshenkoBeam<1, 3, degree, 2 * degree, tests::TimoshenkoCrossParameters>>;

int main()
{
  int failures = 0;
  double previous = 0.;
  for (unsigned int n = 2; n <= 8; n *= 2)
  {
    Loop loop(std::string("domains/cross.geo"));
    loop.set_refinement(n);
    loop.constrain_dirichlet_nodes();
    loop.prepare();

    const unsigned int size = loop.size_of_system();
    std::vector<double> rhs(size), zero(size, 0.);
    loop.residual(zero, rhs);
    for (double& value : rhs)
      value = -value;

    SparseLA::SolverOptions options;
    options.rtol = 1e-12;
    options.max_iterations = 20 * size;
    std::vector<double> x(size, 0.);
    const SparseLA::SolverResult cg =
      SparseLA::cg([&](const std::vector<double>& in, std::vector<double>& out)
                   { loop.apply(in, out); }, rhs, x, options);
    const double error = loop.errors(x)[0];
    tests::record("- {degree: %u, n: %u, dofs: %u, cg: %s, iterations: %u, error: %.6e}\n", degree,
                  n, size, SparseLA::to_string(cg.reason), cg.iterations, error);
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
      const std::vector<double> x_direct = tests::direct_solve(loop, rhs);
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
  return tests::check(failures);
}
