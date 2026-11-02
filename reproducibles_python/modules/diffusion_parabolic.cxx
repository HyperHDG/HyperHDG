// Parabolic LDG-H diffusion on the unit cube, theta scheme: diffusion_convergence_parabolic.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/parabolic.hxx>
#include <HyperHDG/local_solver/diffusion_parab_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

template <unsigned dim, unsigned degree>
using DiffusionParab = GlobalLoop::Parabolic<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::DiffusionParab<dim, degree, 2 * degree, TestParametersSinParab, double> >;

NB_MODULE(diffusion_parabolic, m)
{
  using namespace reproducibles;
  // constructed with the vector of solver constants (e.g. diffusion coefficient, theta, dt)
  bind_grid<DiffusionParab, std::vector<unsigned int>, std::vector<double> >(
    m, "DiffusionParab", seq<1, 2>{}, seq<1, 2, 3>{});
}
