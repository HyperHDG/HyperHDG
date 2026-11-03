// Bilaplacian with the mass-approximation eigenvalue loop on the unit cube:
// bilaplacian_convergence_eigenvalue_*.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/mass_approx_eigenvalue.hxx>
#include <HyperHDG/local_solver/bilaplacian_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/bilaplacian.hxx"

template <unsigned dim, unsigned degree>
using Bilaplacian = GlobalLoop::MassApproxEigenvalue<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::Bilaplacian<dim, degree, 2 * degree, TestParametersEigs, double> >;

NB_MODULE(bilaplacian_mass_approx_eigenvalue, m)
{
  using namespace reproducibles;
  bind_grid<Bilaplacian>(m, "Bilaplacian", seq<1, 2>{}, seq<1, 2, 3>{});
}
