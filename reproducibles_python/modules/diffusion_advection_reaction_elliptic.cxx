// LDG diffusion-advection-reaction on the unit cube, elliptic loop: advection_diffusion_reaction_convergence_elliptic.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_advection_reaction_ldg.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include "reproducibles_python/parameters/diffusion_advection_reaction.hxx"

template <unsigned dim, unsigned degree>
using DiffusionAdvectionReaction = GlobalLoop::Elliptic<
  Topology::Cubic<dim, dim>,
  Geometry::UnitCube<dim, dim, double>,
  NodeDescriptor::Cubic<dim, dim>,
  LocalSolver::DiffusionAdvectionReaction<dim, degree, 2 * degree, TestParametersSinEllipt, double> >;

NB_MODULE(diffusion_advection_reaction_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<DiffusionAdvectionReaction>(m, "DiffusionAdvectionReaction", seq<1, 2>{}, seq<1, 2, 3>{});
}
