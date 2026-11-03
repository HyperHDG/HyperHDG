// BilaplacianEigs with the shifted-inverse eigenvalue loop on the unit cube:
// bilaplacian_convergence_eigenvalue_*.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/shifted_inverse_eigenvalue.hxx>
#include <HyperHDG/local_solver/bilaplacian_eigs_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/bilaplacian.hxx"

template <unsigned dim, unsigned degree>
using BilaplacianEigs = GlobalLoop::ShiftedInverseEigenvalue<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::BilaplacianEigs<dim, degree, 2 * degree, TestParametersEigs, double> >;

NB_MODULE(bilaplacian_shifted_inverse_eigenvalue, m)
{
  using namespace reproducibles;
  bind_grid<BilaplacianEigs>(m, "BilaplacianEigs", seq<1, 2>{}, seq<1, 2, 3>{});
}
