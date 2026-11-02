// LDG-H diffusion on hypergraphs of dim-dimensional hyperedges in 3D, elliptic loop:
// diffusion_hypergraph_convergence_elliptic.py and diffusion_convergence_hypergraph_elliptic.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

template <unsigned dim, unsigned degree>
using Diffusion = GlobalLoop::Elliptic<
  Topology::Cubic<dim, 3>,
  Geometry::UnitCube<dim, 3, double>,
  NodeDescriptor::Cubic<dim, 3>,
  LocalSolver::Diffusion<dim, degree, 6, HG<dim>::template TestParametersQuadEllipt, double> >;

NB_MODULE(diffusion_hypergraph_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<Diffusion>(m, "Diffusion", seq<1, 2, 3>{}, seq<1, 2, 3>{});
}
