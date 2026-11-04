// The condensed system solved three ways: CHOLMOD (sparse Cholesky of the assembled free block),
// CG on the assembled CSR matrix, and matrix-free CG (GlobalLoop::Generic::apply). The solutions
// and errors agree, in 2D and 3D. Built only if CMake finds CHOLMOD (SuiteSparse).
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/parameters/diffusion.hxx>
#include <tests/sparse_direct.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <vector>

using Vec = std::vector<double>;
int failures = 0;

template <unsigned int dim, unsigned int degree>
void test(const unsigned int n)
{
  GlobalLoop::Generic<Topology::Cubic<dim, dim>, Geometry::UnitCube<dim, dim, double>,
                      NodeDescriptor::Cubic<dim, dim>,
                      LocalSolver::Diffusion<dim, degree, 2 * degree, tests::DiffusionSin>>
    loop(std::vector<unsigned int>(dim, n));
  loop.constrain_dirichlet_nodes();
  loop.prepare();
  const tests::FreeSystem system(loop);

  tests::Cholmod cholmod(system);
  cholmod.analyze();
  const bool factorized = cholmod.factorize();
  const Vec x_direct = system.expand(cholmod.solve(system.rhs));

  SparseLA::SolverOptions options;
  options.rtol = 1e-13;
  options.max_iterations = 20 * system.size();
  Vec x_free(system.size(), 0.);
  SparseLA::cg([&](const Vec& in, Vec& out) { system.multiply(in, out); }, system.rhs, x_free,
               options);
  const Vec x_csr = system.expand(x_free);

  const unsigned int size = loop.size_of_system();
  Vec rhs(size), zero(size, 0.), x_free_cg(size, 0.);
  loop.residual(zero, rhs);
  for (double& value : rhs)
    value = -value;
  SparseLA::cg([&](const Vec& in, Vec& out) { loop.apply(in, out); }, rhs, x_free_cg, options);

  double d_csr = 0., d_mf = 0., scale = 0.;
  for (unsigned int i = 0; i < size; ++i)
  {
    d_csr = std::max(d_csr, std::abs(x_csr[i] - x_direct[i]));
    d_mf = std::max(d_mf, std::abs(x_free_cg[i] - x_direct[i]));
    scale = std::max(scale, std::abs(x_direct[i]));
  }
  const double e_direct = loop.errors(x_direct)[0], e_cg = loop.errors(x_free_cg)[0];
  const bool ok = factorized && d_csr <= 1e-9 * scale && d_mf <= 1e-9 * scale &&
                  std::abs(e_direct - e_cg) <= 1e-8 * e_direct;
  std::printf(
    "%s  %uD degree %u, n = %u: %zu free dofs, %zu nonzeros, factor %.0f nonzeros; "
    "CSR-CG %.1e, matrix-free CG %.1e from CHOLMOD; errors %.6e / %.6e\n",
    ok ? "ok  " : "FAIL", dim, degree, n, system.size(), system.nnz(), cholmod.factor_nnz(),
    d_csr / scale, d_mf / scale, e_direct, e_cg);
  failures += !ok;
}

int main()
{
  test<2, 1>(8);
  test<2, 3>(8);
  test<3, 2>(4);
  return failures;
}
