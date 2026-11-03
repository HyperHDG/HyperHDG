// DiffusionEigs with the nonlinear (Newton) eigenvalue loop on hypergraphs of dim-dimensional
// hyperedges in 3D: diffusion_hypergraph_convergence_eigenvalue_*.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/nonlinear_eigenvalue.hxx>
#include <HyperHDG/local_solver/diffusion_eigs_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

template <unsigned dim, unsigned degree>
using DiffusionEigs = GlobalLoop::NonlinearEigenvalue<
  Topology::Cubic<dim, 3>,
  Geometry::UnitCube<dim, 3, double>,
  NodeDescriptor::Cubic<dim, 3>,
  LocalSolver::DiffusionEigs<dim, degree, 2 * degree, TestParametersEigs, double> >;

NB_MODULE(diffusion_hypergraph_nonlinear_eigenvalue, m)
{
  using namespace reproducibles;
  bind_grid<DiffusionEigs>(m, "DiffusionEigs", seq<1, 2>{}, seq<1, 2, 3>{});
}
