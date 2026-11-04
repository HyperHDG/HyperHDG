// LDG-H diffusion on the unit cube with GlobalLoop::Generic: solve u = sin(pi/2 x_0) on meshes of
// 2, 4, 8 (degree 1: up to 16) elements per direction and check that the L2 error converges with
// order degree + 1 on the last refinement. On meshes small enough for a dense LU, the assembled
// matrix is also solved directly: it must be non-singular, and the direct solution must agree with
// the matrix-free CG. One executable per (DIM, DEGREE, SOLVER), set by CMake
// (tests/CMakeLists.txt): SOLVER is Diffusion or DiffusionTensor (sum factorization); REAL = float
// makes the geometry, the solver and the vectors single precision (one case, 2D degree 1, whose
// errors stay well above the round-off floor of the solve in float).
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/local_solver/diffusion_tensor.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/direct_solve.hxx>
#include <tests/parameters/diffusion.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <type_traits>
#include <vector>

// Set per executable by tests/CMakeLists.txt; the defaults only serve editors and linters that
// compile this file without the CMake flags.
#ifndef DIM
#define DIM 2
#endif
#ifndef DEGREE
#define DEGREE 1
#endif
#ifndef REAL
#define REAL double
#endif
#ifndef SOLVER
#define SOLVER Diffusion
#endif
#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

constexpr unsigned int dim = DIM, degree = DEGREE;
using real = REAL;
constexpr bool single = std::is_same_v<real, float>;
constexpr unsigned int max_direct_dofs = 2000;  // dense LU: 32 MB, well under a second
// Degree 1 is cheap, and in 3D its rate is preasymptotic up to n = 8 (1.67 instead of 1.7).
constexpr unsigned int max_n = degree == 1 ? 16 : 8;

using Loop =
  GlobalLoop::Generic<Topology::Cubic<dim, dim>,
                      Geometry::UnitCube<dim, dim, real>,
                      NodeDescriptor::Cubic<dim, dim>,
                      LocalSolver::SOLVER<dim, degree, 2 * degree, tests::DiffusionSin, real>,
                      real>;

int main()
{
  int failures = 0;
  double previous = 0.;
  for (unsigned int n = 2; n <= max_n; n *= 2)
  {
    Loop loop(std::vector<unsigned int>(dim, n));
    loop.constrain_dirichlet_nodes();
    loop.prepare();

    const unsigned int size = loop.size_of_system();
    std::vector<real> rhs(size), zero(size, 0.);
    loop.residual(zero, rhs);
    for (real& value : rhs)
      value = -value;

    SparseLA::SolverOptions options;
    options.rtol = single ? 1e-6 : 1e-12;
    options.max_iterations = 10 * size;
    std::vector<real> x(size, 0.);
    const SparseLA::SolverResult cg =
      SparseLA::cg([&](const std::vector<real>& in, std::vector<real>& out)
                   { loop.apply(in, out); }, rhs, x, options);
    const double error = loop.errors(x)[0];
    std::printf(
      "- {solver: %s, dim: %u, degree: %u, n: %u, dofs: %u, cg: %s, iterations: %u, "
      "error: %.6e}\n",
      STRINGIFY(SOLVER), dim, degree, n, size, SparseLA::to_string(cg.reason), cg.iterations,
      error);
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
      const std::vector<real> x_direct = tests::direct_solve(loop, rhs);
      double difference = 0., scale = 0.;
      for (unsigned int i = 0; i < size; ++i)
      {
        difference = std::max(difference, (double)std::abs(x_direct[i] - x[i]));
        scale = std::max(scale, (double)std::abs(x_direct[i]));
      }
      const double tolerance = single ? 1e-4 : 1e-8;
      std::printf("  direct vs cg: %.2e (relative, at most %.0e)\n", difference / scale, tolerance);
      if (!(difference <= tolerance * scale))
        ++failures;
    }

    // Only the last rate is checked: the coarse meshes are preasymptotic.
    if (n == max_n)
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
