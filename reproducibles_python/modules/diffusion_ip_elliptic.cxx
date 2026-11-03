// Interior penalty diffusion on the unit cube, elliptic loop: the classes of
// diffusion_convergence_elliptic_ip.py. Constructed with the penalty as second argument.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ip.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

template <unsigned dim, unsigned degree>
using DiffusionIP = GlobalLoop::Elliptic<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::DiffusionIP<dim, degree, 2 * degree, TestParametersSinEllipt, double> >;

NB_MODULE(diffusion_ip_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<DiffusionIP>(m, "DiffusionIP", seq<1, 2>{}, seq<1, 2, 3>{});
}
