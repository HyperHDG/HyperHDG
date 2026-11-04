// The Timoshenko beam solver on GlobalLoop::Generic, on the cross of domains/cross.geo with the
// manufactured solution of tests/parameters/timoshenko.hxx: the local operator (blocks probed from
// the legacy maps) is that of GlobalLoop::Elliptic with the legacy interface (apply and residual
// agree on the free dofs to round-off, for both layouts), and CG solutions converge with the
// expected order on refined beams, degrees 1 to 3.
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/file.hxx>

#include <tests/parameters/timoshenko.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using Vec = std::vector<double>;
int failures = 0;
void check(const bool ok, const std::string& what)
{
  std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what.c_str());
  failures += !ok;
}

template <unsigned int degree, class LayoutT = LocalOperator::Condensed>
using Solver = LocalSolver::
  TimoshenkoBeam<1, 3, degree, 2 * degree, tests::TimoshenkoCrossParameters, double, LayoutT>;
template <class SolverT>
using Loop = GlobalLoop::
  Generic<Topology::File<1, 3>, Geometry::File<1, 3>, NodeDescriptor::File<1, 3>, SolverT>;
const std::string domain = "domains/cross.geo";

// Generic against Elliptic: apply and residual on the free dofs (Elliptic zeroes the Dirichlet
// rows, Generic makes them identity rows).
template <unsigned int degree, class LayoutT>
void compare(const char* layout)
{
  GlobalLoop::Elliptic<Topology::File<1, 3>, Geometry::File<1, 3>, NodeDescriptor::File<1, 3>,
                       Solver<degree, LayoutT>, Vec>
    elliptic(domain, 1.);
  Loop<Solver<degree, LayoutT>> generic(domain);
  elliptic.set_refinement(2);
  generic.set_refinement(2);
  generic.constrain_dirichlet_nodes();
  const unsigned int n = generic.size_of_system();
  std::vector<bool> constrained(n, false);
  for (unsigned int i : generic.constraints())
    constrained[i] = true;
  std::mt19937 generator(degree);
  std::uniform_real_distribution<double> uniform(-1., 1.);
  Vec x(n), y(n), r(n);
  for (double& value : x)
    value = uniform(generator);
  generic.prepare();
  generic.apply(x, y);
  generic.residual(x, r);
  const Vec y_legacy = elliptic.trace_to_flux(x), r_legacy = elliptic.residual_flux(x);
  double dy = 0., dr = 0., sy = 0., sr = 0.;
  for (unsigned int i = 0; i < n; ++i)
    if (!constrained[i])
    {
      dy = std::max(dy, std::abs(y[i] - y_legacy[i]));
      dr = std::max(dr, std::abs(r[i] - r_legacy[i]));
      sy = std::max(sy, std::abs(y_legacy[i]));
      sr = std::max(sr, std::abs(r_legacy[i]));
    }
  char what[160];
  std::snprintf(what, sizeof what, "degree %u, %s: apply = Elliptic (%.1e), residual (%.1e)",
                degree, layout, dy / sy, dr / sr);
  check(dy <= 1e-12 * sy && dr <= 1e-12 * sr, what);
}

// The convergence ladder: CG solutions on 2, 4, 8 subintervals per beam, last rate.
template <unsigned int degree>
void converge()
{
  double previous = 0.;
  for (unsigned int n = 2; n <= 8; n *= 2)
  {
    Loop<Solver<degree>> loop(domain);
    loop.set_refinement(n);
    loop.constrain_dirichlet_nodes();
    loop.prepare();
    const unsigned int size = loop.size_of_system();
    Vec rhs(size), zero(size, 0.), x(size, 0.);
    loop.residual(zero, rhs);
    for (double& value : rhs)
      value = -value;
    SparseLA::SolverOptions options;
    options.rtol = 1e-12;
    options.max_iterations = 20 * size;
    const auto cg =
      SparseLA::cg([&](const Vec& in, Vec& out) { loop.apply(in, out); }, rhs, x, options);
    const double error = loop.errors(x)[0];
    std::printf("- {degree: %u, n: %u, dofs: %u, cg: %s, iterations: %u, error: %.6e}\n", degree, n,
                size, SparseLA::to_string(cg.reason), cg.iterations, error);
    if (!cg.converged())
      check(false, "CG converges");
    if (n == 8)
    {
      const double rate = std::log2(previous / error);
      char what[96];
      std::snprintf(what, sizeof what, "degree %u: rate %.3f (at least %.1f)", degree, rate,
                    degree + 1 - 0.3);
      check(rate >= degree + 1 - 0.3, what);
    }
    previous = error;
  }
}

int main()
{
  compare<1, LocalOperator::Condensed>("Condensed");
  compare<2, LocalOperator::Condensed>("Condensed");
  compare<2, LocalOperator::Factored>("Factored");
  converge<1>();
  converge<2>();
  converge<3>();
  return failures;
}
