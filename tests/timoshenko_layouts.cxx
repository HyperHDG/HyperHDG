// The layouts of the local operator of the Timoshenko beam solver against each other, on beams of
// different lengths and directions in space (domains/skew_beams.geo, run in the source tree):
// BeamStiffness (the default: the frame of a beam, six stiffnesses, the fluxes of loads) must give
// the apply, the residual and the errors of Condensed (S and K from the whole local matrix) up to
// round-off, with and without prepare. One executable per DEGREE, set by CMake
// (tests/CMakeLists.txt).
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include <tests/parameters/timoshenko.hxx>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

// Set per executable by tests/CMakeLists.txt; the default only serves editors and linters that
// compile this file without the CMake flags.
#ifndef DEGREE
#define DEGREE 1
#endif

template <class LayoutT>
using Loop = GlobalLoop::Generic<
  Topology::File<1, 3>,
  Geometry::File<1, 3>,
  NodeDescriptor::File<1, 3>,
  LocalSolver::
    TimoshenkoBeam<1, 3, DEGREE, 2 * DEGREE, tests::TimoshenkoCrossParameters, double, LayoutT>>;
using Vec = std::vector<double>;

int failures = 0;
void check(const bool ok, const char* what)
{
  std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what);
  failures += !ok;
}
bool agree(const Vec& a, const Vec& b, const double tol)
{
  double difference = 0., scale = 0.;
  for (std::size_t i = 0; i < a.size(); ++i)
  {
    difference = std::max(difference, std::abs(a[i] - b[i]));
    scale = std::max(scale, std::abs(b[i]));
  }
  return difference <= tol * scale;
}

int main()
{
  static_assert(sizeof(Loop<LocalSolver::BeamStiffness>::Operator) <=
                  sizeof(Loop<LocalOperator::Condensed>::Operator) / 3,
                "BeamStiffness is much smaller than Condensed.");
  static_assert(offsetof(Loop<LocalSolver::BeamStiffness>::Operator, load_to_flux_single) == 128,
                "apply reads two cache lines: a frame, a length and six stiffnesses.");
  for (const bool prepared : {false, true})
  {
    std::printf("---- %s\n", prepared ? "prepared" : "without prepare");
    Loop<LocalSolver::BeamStiffness> stiffness(std::string("domains/skew_beams.geo"));
    Loop<LocalOperator::Condensed> condensed(std::string("domains/skew_beams.geo"));
    stiffness.set_refinement(3), condensed.set_refinement(3);
    stiffness.constrain_dirichlet_nodes(), condensed.constrain_dirichlet_nodes();
    if (prepared)
      stiffness.prepare(), condensed.prepare();

    const unsigned int n = stiffness.size_of_system();
    std::mt19937 generator(1);
    std::uniform_real_distribution<double> uniform(-1., 1.);
    Vec x(n), y_stiffness(n), y_condensed(n);
    for (double& value : x)
      value = uniform(generator);

    stiffness.apply(x, y_stiffness), condensed.apply(x, y_condensed);
    check(agree(y_stiffness, y_condensed, 1e-12), "apply: BeamStiffness = Condensed");
    stiffness.residual(x, y_stiffness), condensed.residual(x, y_condensed);
    check(agree(y_stiffness, y_condensed, 1e-12), "residual: BeamStiffness = Condensed");
    const double error_stiffness = GlobalLoop::errors(stiffness, x)[0],
                 error_condensed = GlobalLoop::errors(condensed, x)[0];
    check(std::abs(error_stiffness - error_condensed) <= 1e-12 * error_condensed,
          "errors: BeamStiffness = Condensed");
  }
  return failures;
}
