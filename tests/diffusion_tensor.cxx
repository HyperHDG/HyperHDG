// LocalSolver::DiffusionTensor (sum-factorized, fast diagonalization) against
// LocalSolver::Diffusion on the same mesh (3 elements per direction): the same discretization, so
// apply, residual and errors agree to round-off, prepared and streaming, and a CG solve gives the
// same error. Degrees 1 to 3 in dimension DIM (set by tests/CMakeLists.txt); DIM = 1 also on
// hyperedges in 2D.
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/local_solver/diffusion_tensor.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <tests/parameters/diffusion.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <vector>

#ifndef DIM
#define DIM 2  // set by tests/CMakeLists.txt; the default serves editors and linters
#endif

using Vec = std::vector<double>;
int failures = 0;
void check(const bool ok,
           const char* what,
           const unsigned int edge_dim,
           const unsigned int space_dim,
           const unsigned int degree,
           const double deviation)
{
  std::printf("%s  %u-%u degree %u: %s (%.1e)\n", ok ? "ok  " : "FAIL", edge_dim, space_dim, degree,
              what, deviation);
  failures += !ok;
}
double relative(const Vec& a, const Vec& b)
{
  double difference = 0., scale = 0.;
  for (std::size_t i = 0; i < a.size(); ++i)
  {
    difference = std::max(difference, std::abs(a[i] - b[i]));
    scale = std::max(scale, std::abs(b[i]));
  }
  return difference / scale;
}

template <unsigned int edge_dim, unsigned int space_dim, unsigned int degree, class SolverT>
using Loop = GlobalLoop::Generic<Topology::Cubic<edge_dim, space_dim>,
                                 Geometry::UnitCube<edge_dim, space_dim, double>,
                                 NodeDescriptor::Cubic<edge_dim, space_dim>,
                                 SolverT>;

template <class LoopT>
Vec product(LoopT& loop, const Vec& x)
{
  Vec y(x.size());
  loop.apply(x, y);
  return y;
}
template <class LoopT>
Vec residual_of(LoopT& loop, const Vec& x)
{
  Vec r(x.size());
  loop.residual(x, r);
  return r;
}
template <class LoopT>
double solve(LoopT& loop)
{
  const unsigned int size = loop.size_of_system();
  Vec rhs(size), zero(size, 0.), x(size, 0.);
  loop.residual(zero, rhs);
  for (double& value : rhs)
    value = -value;
  SparseLA::SolverOptions options;
  options.rtol = 1e-12;
  SparseLA::cg([&](const Vec& in, Vec& out) { loop.apply(in, out); }, rhs, x, options);
  return loop.errors(x)[0];
}

template <unsigned int edge_dim, unsigned int space_dim, unsigned int degree>
void test()
{
  using Reference = Loop<edge_dim, space_dim, degree,
                         LocalSolver::Diffusion<edge_dim, degree, 2 * degree, tests::DiffusionSin>>;
  using Tensor =
    Loop<edge_dim, space_dim, degree,
         LocalSolver::DiffusionTensor<edge_dim, degree, 2 * degree, tests::DiffusionSin>>;
  const std::vector<unsigned int> elements(space_dim, 3);
  Reference reference(elements);
  Tensor tensor(elements);
  reference.constrain_dirichlet_nodes();
  tensor.constrain_dirichlet_nodes();

  const unsigned int n = reference.size_of_system();
  std::mt19937 generator(degree);
  std::uniform_real_distribution<double> uniform(-1., 1.);
  Vec x(n);
  for (double& value : x)
    value = uniform(generator);

  // streaming (no prepare): the tensor operator is built per sweep
  const Vec y_reference = product(reference, x), r_reference = residual_of(reference, x);
  double deviation = relative(product(tensor, x), y_reference);
  check(deviation <= 1e-12, "apply, streaming", edge_dim, space_dim, degree, deviation);
  deviation = relative(residual_of(tensor, x), r_reference);
  check(deviation <= 1e-12, "residual, streaming", edge_dim, space_dim, degree, deviation);
  // prepared (the elements share one operator)
  reference.prepare();
  tensor.prepare();
  deviation = relative(product(tensor, x), y_reference);
  check(deviation <= 1e-12, "apply, prepared", edge_dim, space_dim, degree, deviation);
  deviation = relative(residual_of(tensor, x), r_reference);
  check(deviation <= 1e-12, "residual, prepared", edge_dim, space_dim, degree, deviation);
  const double e_reference = reference.errors(x)[0], e_tensor = tensor.errors(x)[0];
  deviation = std::abs(e_tensor - e_reference) / e_reference;
  check(deviation <= 1e-12, "errors", edge_dim, space_dim, degree, deviation);
  // a solve
  const double s_reference = solve(reference), s_tensor = solve(tensor);
  deviation = std::abs(s_tensor - s_reference) / s_reference;
  check(deviation <= 1e-9, "error of the CG solution", edge_dim, space_dim, degree, deviation);
}

int main()
{
  test<DIM, DIM, 1>();
  test<DIM, DIM, 2>();
  test<DIM, DIM, 3>();
  if constexpr (DIM == 1)
  {
    test<1, 2, 1>();
    test<1, 2, 3>();
  }
  return failures;
}
