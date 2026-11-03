// Diffusion with the mass-approximation eigenvalue loop on the unit cube:
// diffusion_convergence_eigenvalue_*.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/mass_approx_eigenvalue.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

template <unsigned dim, unsigned degree>
using Diffusion = GlobalLoop::MassApproxEigenvalue<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::Diffusion<dim, degree, 2 * degree, TestParametersEigs, double> >;

NB_MODULE(diffusion_mass_approx_eigenvalue, m)
{
  using namespace reproducibles;
  bind_grid<Diffusion>(m, "Diffusion", seq<1, 2>{}, seq<1, 2, 3>{});
}
