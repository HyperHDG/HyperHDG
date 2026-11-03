// Timoshenko beams with the default parameters on the network of domains/cross.geo, elliptic loop:
// timoshenko_convergence_cross.py.
#include "bind_grid.hxx"

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

// dim is the hyperedge dimension (1, beams) in three-dimensional space
template <unsigned dim, unsigned degree>
using TimoshenkoBeam =
  GlobalLoop::Elliptic<Topology::File<dim, 3>,
                       Geometry::File<dim, 3>,
                       NodeDescriptor::File<dim, 3>,
                       LocalSolver::TimoshenkoBeam<dim, 3, degree, 2 * degree> >;

NB_MODULE(timoshenko_file_elliptic, m)
{
  using namespace reproducibles;
  bind_grid<TimoshenkoBeam, std::string>(m, "TimoshenkoBeam", seq<1>{}, seq<1, 2, 3, 4, 5, 6>{});
}
