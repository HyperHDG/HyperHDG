// LDG-H bilaplacian on hypergraphs in 3D, elliptic loop:
// bilaplacian_hypergraph_convergence_elliptic.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/bilaplacian_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/bilaplacian.hxx"

template <unsigned dim, unsigned degree>
using Bilaplacian = GlobalLoop::Elliptic<
  Topology::Cubic<dim, 3>,
  Geometry::UnitCube<dim, 3, double>,
  NodeDescriptor::Cubic<dim, 3>,
  LocalSolver::
    Bilaplacian<dim, degree, 2 * degree, HG<dim>::template TestParametersQuadEllipt, double> >;

NB_MODULE(bilaplacian_hypergraph_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<Bilaplacian>(m, "Bilaplacian", seq<1, 2>{}, seq<1, 2, 3>{});
}
