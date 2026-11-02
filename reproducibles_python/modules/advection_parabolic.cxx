// Parabolic advection on the unit cube, theta scheme: advection_convergence_parabolic.py.
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
  LocalSolver::AdvectionParab<dim, degree, 2 * degree, TestParametersSinParab, double> >;

NB_MODULE(advection_parabolic, m)
{
  using namespace reproducibles;
  // constructed with the vector of solver constants (e.g. diffusion coefficient, theta, dt)
  bind_grid<AdvectionParab, std::vector<unsigned int>, std::vector<double> >(
    m, "AdvectionParab", seq<1, 2>{}, seq<1, 2, 3>{});
}
