// Parabolic advection with injection on the network of domains/injection_test.geo:
// advection_hypergraph.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/parabolic.hxx>
#include <HyperHDG/local_solver/advection_parab_ldgh.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include "reproducibles_python/parameters/advection.hxx"

// dim is the hyperedge dimension of the domain file (domains/injection_test.geo)
template <unsigned dim, unsigned degree>
using AdvectionParab = GlobalLoop::Parabolic<
  Topology::File<dim, 2>,
  Geometry::File<dim, 2>,
  NodeDescriptor::File<dim, 2>,
  LocalSolver::AdvectionParab<dim, degree, 2 * degree, injection, double> >;

NB_MODULE(advection_injection_file_parabolic, m)
{
  using namespace reproducibles;
  // constructed with the file name and the vector of solver constants (tau, theta, dt)
  bind_grid<AdvectionParab, std::string, std::vector<double> >(
    m, "AdvectionParab", seq<1>{}, seq<3>{});
}
