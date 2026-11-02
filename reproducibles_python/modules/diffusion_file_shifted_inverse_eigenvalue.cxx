// DiffusionEigs with the shifted-inverse eigenvalue loop on a network read from file (domains/cross_2d.geo): the *_cross.py eigenvalue scripts.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/shifted_inverse_eigenvalue.hxx>
#include <HyperHDG/local_solver/diffusion_eigs_ldgh.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include "reproducibles_python/parameters/diffusion.hxx"

// dim is the hyperedge dimension; the network file is two-dimensional
template <unsigned dim, unsigned degree>
using DiffusionEigs = GlobalLoop::ShiftedInverseEigenvalue<
  Topology::File<dim, 2>,
  Geometry::File<dim, 2>,
  NodeDescriptor::File<dim, 2>,
  LocalSolver::DiffusionEigs<dim, degree, 2 * degree, TestParametersEigs, double> >;

NB_MODULE(diffusion_file_shifted_inverse_eigenvalue, m)
{
  using namespace reproducibles;
  bind_grid<DiffusionEigs, std::string>(m, "DiffusionEigs", seq<1>{}, seq<1, 2, 3>{});
}
