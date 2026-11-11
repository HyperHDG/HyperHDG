// What the layouts of the local operator of the Timoshenko beam solver cost per beam on a large
// network: BeamStiffness (the default: the frame of a beam, six stiffnesses, the fluxes of loads)
// against Condensed (S and K) and Factored (the factors of the local matrix). Per layout:
// prepare, apply, residual and assemble on the store, the size of the store, and the agreement of
// apply and residual with the first layout. Degrees 1 to 3.
//
// The network is large on purpose: the store of Condensed and Factored exceeds the caches, so
// their apply pays for memory traffic. By default the cross of domains/cross.geo, refined to
// 100 000 beams (all along the axes, unit stiffnesses); with a domain file as argument, that
// network unrefined, e.g. a mikado network with the stiffnesses of its fibers
// (.geo.h5: configure with HDF5), from https://github.com/HyperHDG/fiber_network:
//   python fiber_network/make_geo.py --mikado 400 --seed 1 --dirichlet xmin xmax ymin ymax \
//     -o mikado.geo.h5        # 109 561 beams
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include <tests/performance/benchmark.hxx>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

using Vec = std::vector<double>;
std::string domain = "domains/cross.geo";
unsigned int refinement = 25000;  // of the default domain: 4 beams each

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

template <unsigned int degree, class LayoutT>
void measure(ankerl::nanobench::Bench& bench, const std::string& name, Vec& y_first, Vec& r_first)
{
  using ankerl::nanobench::doNotOptimizeAway;
  using Loop = GlobalLoop::Generic<
    Topology::File<1, 3>, Geometry::File<1, 3>, NodeDescriptor::File<1, 3>,
    LocalSolver::TimoshenkoBeam<1, 3, degree, 2 * degree,
                                LocalSolver::TimoschenkoBeamParametersDefault, double, LayoutT>>;
  Loop loop(domain);
  if (refinement > 1)
    loop.set_refinement(refinement);
  loop.constrain_dirichlet_nodes();
  const unsigned int size = loop.size_of_system();
  std::mt19937 generator(0);
  std::uniform_real_distribution<double> uniform(-1., 1.);
  Vec x(size), y(size), r(size);
  for (double& value : x)
    value = uniform(generator);

  bench.batch(loop.n_hyEdges());
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
              loop.residual(x, r);
              doNotOptimizeAway(r.data());
            });
  std::vector<unsigned int> rows(loop.n_triplets()), cols(rows.size());
  Vec values(rows.size());
  bench.run(name + ": assemble",
            [&]
            {
              loop.assemble(typename Loop::Args(), false, rows, cols, values);
              doNotOptimizeAway(values.data());
            });
  if (y_first.empty())
    y_first = y, r_first = r;
  std::printf(
    "  %s: %zu beams, %u dofs, store %.1f MB (%zu bytes per beam); against the first "
    "layout: apply %.1e, residual %.1e\n",
    name.c_str(), loop.n_hyEdges(), size,
    loop.n_operators() * sizeof(typename Loop::Operator) / 1048576.,
    sizeof(typename Loop::Operator), relative_difference(y, y_first),
    relative_difference(r, r_first));
}

template <unsigned int degree>
void measure_degree()
{
  ankerl::nanobench::Bench bench;
  bench.title("Timoshenko beams, degree " + std::to_string(degree))
    .unit("beam")
    .epochs(5)
    .epochIterations(1)
    .warmup(1)
    .relative(false);
  Vec y_first, r_first;
  measure<degree, LocalOperator::Condensed>(bench, "Condensed", y_first, r_first);
  measure<degree, LocalSolver::BeamStiffness>(bench, "BeamStiffness", y_first, r_first);
  measure<degree, LocalOperator::Factored>(bench, "Factored", y_first, r_first);
  tests::performance::save(bench, "timoshenko_layouts_p" + std::to_string(degree));
}

int main(int argc, char** argv)
{
  if (argc > 1)
  {
    domain = argv[1];
    refinement = 1;
  }
  measure_degree<1>();
  measure_degree<2>();
  measure_degree<3>();
}
