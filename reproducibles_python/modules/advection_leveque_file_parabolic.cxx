// Parabolic advection with LeVeque's rotating field on the hypergraph of domains/leVeque_hg.geo:
// advection_convergence_leVeque_hypergraph.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/parabolic.hxx>
#include <HyperHDG/local_solver/advection_parab_ldgh.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include "reproducibles_python/parameters/advection.hxx"

// dim is the hyperedge dimension of the domain file (domains/leVeque_hg.geo)
template <unsigned dim, unsigned degree>
using AdvectionParab = GlobalLoop::Parabolic<
  Topology::File<dim, 3>,
  Geometry::File<dim, 3>,
  NodeDescriptor::File<dim, 3>,
  LocalSolver::AdvectionParab<dim, degree, 2 * degree, LeVequeHG, double> >;

NB_MODULE(advection_leveque_file_parabolic, m)
{
  using namespace reproducibles;
  // constructed with the file name and the vector of solver constants (tau, theta, dt)
  bind_grid<AdvectionParab, std::string, std::vector<double> >(
    m, "AdvectionParab", seq<2>{}, seq<0, 1, 2>{});
}
