#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <examples/parameters/diffusion.hxx>

int main()
{
  const unsigned int poly_degree = 1;
  const unsigned int hyEdge_dim = 1;
  const unsigned int space_dim = 2;
  const unsigned int refinement = 1;

  GlobalLoop::Elliptic<Topology::Cubic<hyEdge_dim, space_dim>,
                       Geometry::UnitCube<hyEdge_dim, space_dim>,
                       NodeDescriptor::Cubic<hyEdge_dim, space_dim>,
                       LocalSolver::Diffusion<hyEdge_dim, poly_degree, 2 * poly_degree,
                                              HG<hyEdge_dim>::DiffusionElliptic, double> >
    HDG_wrapper(std::vector<unsigned int>(space_dim, 1 << refinement));

  std::vector<double> vectorRHS = HDG_wrapper.residual_flux(HDG_wrapper.zero_vector());
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)
    vectorRHS[i] *= -1.;

  std::vector<double> vectorSolution = SparseLA::conjugate_gradient(vectorRHS, HDG_wrapper);

  std::cout << "Error: " << HDG_wrapper.errors(vectorSolution)[0] << std::endl;

  HDG_wrapper.plot_option("fileName", "diffusion_elliptic_c++");
  HDG_wrapper.plot_option("printFileNumber", "false");
  HDG_wrapper.plot_option("scale", "0.95");
  HDG_wrapper.plot_solution(vectorSolution);

  return 0;
}
