// Interior penalty diffusion on hypergraphs in 3D, elliptic loop:
// diffusion_hypergraph_convergence_elliptic_ip.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ip.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

template <unsigned dim, unsigned degree>
using DiffusionIP = GlobalLoop::Elliptic<
  Topology::Cubic<dim, 3>,
  Geometry::UnitCube<dim, 3, double>,
  NodeDescriptor::Cubic<dim, 3>,
  LocalSolver::
    DiffusionIP<dim, degree, 2 * degree, HG<dim>::template TestParametersQuadEllipt, double> >;

NB_MODULE(diffusion_ip_hypergraph_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<DiffusionIP>(m, "DiffusionIP", seq<1, 2, 3>{}, seq<1, 2, 3>{});
}
