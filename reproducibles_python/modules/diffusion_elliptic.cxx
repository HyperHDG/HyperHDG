// LDG-H diffusion on the unit cube, elliptic loop: the classes of diffusion_convergence_elliptic.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

template <unsigned dim, unsigned degree>
using Diffusion = GlobalLoop::Elliptic<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::Diffusion<dim, degree, 2 * degree, TestParametersSinEllipt, double> >;

NB_MODULE(diffusion_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<Diffusion>(m, "Diffusion", seq<1, 2>{}, seq<1, 2, 3>{});
}
