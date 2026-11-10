// LDG-H diffusion on a graph: the edges of the unit square's 2 x 2 grid, degree 1. The condensed
// system of GlobalLoop::Generic is solved with matrix-free CG (SparseLA::cg), the error against
// the analytic solution is printed, and the solution is plotted to
// output/diffusion_elliptic_c++.vtu.
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <examples/parameters/diffusion.hxx>

#include <iostream>
#include <vector>

int main()
{
  const unsigned int poly_degree = 1;
  const unsigned int hyEdge_dim = 1;
  const unsigned int space_dim = 2;
  const unsigned int refinement = 1;

  GlobalLoop::Generic<Topology::Cubic<hyEdge_dim, space_dim>,
                      Geometry::UnitCube<hyEdge_dim, space_dim>,
                      NodeDescriptor::Cubic<hyEdge_dim, space_dim>,
                      LocalSolver::Diffusion<hyEdge_dim, poly_degree, 2 * poly_degree,
                                             HG<hyEdge_dim>::DiffusionElliptic, double> >
    hdg(std::vector<unsigned int>(space_dim, 1 << refinement));
  hdg.constrain_dirichlet_nodes();
  hdg.prepare();  // keep the local operators: every CG iteration applies them

  // The residual at zero is minus the right-hand side.
  const unsigned int size = hdg.size_of_system();
  std::vector<double> rhs(size), solution(size, 0.);
  hdg.residual(solution, rhs);
  for (double& value : rhs)
    value = -value;

  const SparseLA::SolverResult cg = SparseLA::cg(
    [&](const std::vector<double>& x, std::vector<double>& y) { hdg.apply(x, y); }, rhs, solution);
  std::cout << "CG: " << SparseLA::to_string(cg.reason) << " after " << cg.iterations
            << " iterations" << std::endl;
  std::cout << "Error: " << GlobalLoop::errors(hdg, solution)[0] << std::endl;

  hdg.plot_option("fileName", "diffusion_elliptic_c++");
  hdg.plot_option("printFileNumber", "false");
  hdg.plot_option("scale", "0.95");
  hdg.plot_solution(solution);

  return cg.converged() ? 0 : 1;
}
