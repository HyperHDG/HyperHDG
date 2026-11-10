// Layouts of the local operator of LocalSolver::TimoshenkoBeam on a beam network: time of
// prepare, apply, residual and assemble, the size of the store, and the agreement of the results
// with those of the first layout. Run by oploop-01-timoshenko-layouts.sh.
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#ifndef DEGREE
#define DEGREE 1
#endif

using Vec = std::vector<double>;
std::string domain;
Vec y_reference, r_reference;

double seconds(const std::chrono::steady_clock::time_point start)
{
  return std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
}
double relative_difference(const Vec& a, const Vec& b)
{
  double difference = 0., scale = 0.;
  for (std::size_t i = 0; i < a.size(); ++i)
  {
    difference = std::max(difference, std::abs(a[i] - b[i]));
    scale = std::max(scale, std::abs(b[i]));
  }
  return difference / scale;
}

template <class LayoutT>
void run(const char* name, const unsigned int n_apply)
{
  using Loop = GlobalLoop::Generic<
    Topology::File<1, 3>, Geometry::File<1, 3>, NodeDescriptor::File<1, 3>,
    LocalSolver::TimoshenkoBeam<1, 3, DEGREE, 2 * DEGREE,
                                LocalSolver::TimoschenkoBeamParametersDefault, double, LayoutT>>;
  Loop loop(domain);
  loop.constrain_dirichlet_nodes();
  const unsigned int size = loop.size_of_system();
  std::mt19937 generator(1);
  std::uniform_real_distribution<double> uniform(-1., 1.);
  Vec x(size), y(size), r(size);
  for (double& value : x)
    value = uniform(generator);

  auto start = std::chrono::steady_clock::now();
  loop.prepare();
  const double t_prepare = seconds(start);

  start = std::chrono::steady_clock::now();
  for (unsigned int i = 0; i < n_apply; ++i)
    loop.apply(x, y);
  const double t_apply = seconds(start) / n_apply;

  start = std::chrono::steady_clock::now();
  loop.residual(x, r);
  const double t_residual = seconds(start);

  const std::size_t n_triplets = loop.n_triplets(true);
  std::vector<unsigned int> rows(n_triplets), cols(n_triplets);
  Vec vals(n_triplets);
  start = std::chrono::steady_clock::now();
  loop.assemble(typename Loop::Args(), true, rows, cols, vals);
  const double t_assemble = seconds(start);

  if (y_reference.empty())
    y_reference = y, r_reference = r;
  std::printf(
    "- {layout: %s, degree: %u, beams: %zu, dofs: %u, prepare_s: %.4f, apply_s: %.5f, "
    "residual_s: %.4f, assemble_s: %.4f, store_MB: %.2f, operator_bytes: %zu, "
    "apply_diff: %.1e, residual_diff: %.1e}\n",
    name, (unsigned int)DEGREE, loop.n_hyEdges(), size, t_prepare, t_apply, t_residual, t_assemble,
    loop.n_operators() * sizeof(typename Loop::Operator) / 1048576.,
    sizeof(typename Loop::Operator), relative_difference(y, y_reference),
    relative_difference(r, r_reference));
  std::fflush(stdout);
}

int main(int argc, char** argv)
{
  domain = argc > 1 ? argv[1] : "domains/cross.geo";
  const unsigned int n_apply = argc > 2 ? std::stoi(argv[2]) : 20;
  run<LocalOperator::Condensed>("Condensed", n_apply);
  run<LocalSolver::BeamStiffness>("BeamStiffness", n_apply);
  run<LocalOperator::Factored>("Factored", n_apply);
  return 0;
}
