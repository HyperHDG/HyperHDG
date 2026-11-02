// LDG-H diffusion on hypergraphs of 2D hyperedges in 3D with the high-order test solution, elliptic
// loop: diffusion_convergence_hypergraph_elliptic_high_order.py.
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
  LocalSolver::Diffusion<dim, degree, 6, HG<dim>::template TestParametersQuadElliptHO, double> >;

NB_MODULE(diffusion_hypergraph_elliptic_ho, m)
{
  using namespace reproducibles;
  bind_grid<Diffusion>(m, "Diffusion", seq<2>{}, seq<0, 1, 2, 3>{});
}
