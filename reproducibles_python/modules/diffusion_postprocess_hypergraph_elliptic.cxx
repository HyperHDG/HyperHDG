// Postprocessed LDG-H diffusion on hypergraphs in 3D, elliptic loop:
// diffusion_hypergraph_convergence_elliptic_postprocess.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh_pp.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

template <unsigned dim, unsigned degree>
using DiffusionPostprocess =
  GlobalLoop::Elliptic<Topology::Cubic<dim, 3>,
                       Geometry::UnitCube<dim, 3, double>,
                       NodeDescriptor::Cubic<dim, 3>,
                       LocalSolver::DiffusionPostprocess<dim,
                                                         degree,
                                                         2 * degree,
                                                         HG<dim>::template TestParametersQuadEllipt,
                                                         double> >;

NB_MODULE(diffusion_postprocess_hypergraph_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<DiffusionPostprocess>(m, "DiffusionPostprocess", seq<1, 2, 3>{}, seq<1, 2, 3>{});
}
