// LDG-H bilaplacian on the unit cube, elliptic loop: bilaplacian_convergence_elliptic.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/bilaplacian_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/bilaplacian.hxx"

template <unsigned dim, unsigned degree>
using Bilaplacian = GlobalLoop::Elliptic<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::Bilaplacian<dim, degree, 2 * degree, TestParametersSinEllipt, double> >;

NB_MODULE(bilaplacian_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<Bilaplacian>(m, "Bilaplacian", seq<1, 2>{}, seq<1, 2, 3>{});
}
