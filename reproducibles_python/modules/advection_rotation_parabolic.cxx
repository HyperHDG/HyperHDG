// Parabolic advection on the unit square with LeVeque's rotating velocity field:
// advection_convergence_rotation.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/parabolic.hxx>
#include <HyperHDG/local_solver/advection_parab_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/advection.hxx"

template <unsigned dim, unsigned degree>
using AdvectionParab = GlobalLoop::Parabolic<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::AdvectionParab<dim, degree, 2 * degree, LeVeque, double> >;

NB_MODULE(advection_rotation_parabolic, m)
{
  using namespace reproducibles;
  // constructed with the vector of solver constants (e.g. diffusion coefficient, theta, dt)
  bind_grid<AdvectionParab, std::vector<unsigned int>, std::vector<double> >(
    m, "AdvectionParab", seq<2>{}, seq<0, 1, 2>{});
}
