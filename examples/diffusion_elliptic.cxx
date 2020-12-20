#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <reproducables_python/parameters/diffusion.hxx>

using namespace std;

void diffusion_elliptic()
{
  const unsigned int poly_degree = 1;
  const unsigned int hyEdge_dim = 1;
  const unsigned int space_dim = 2;
  const unsigned int refinement = 1;

  const vector<unsigned int> num_elements(space_dim, 1 << refinement);

  GlobalLoop::Elliptic<Topology::Cubic<hyEdge_dim, space_dim>,
                       Geometry::UnitCube<hyEdge_dim, space_dim>,
                       NodeDescriptor::Cubic<hyEdge_dim, space_dim>,
                       LocalSolver::Diffusion<hyEdge_dim, poly_degree, 2 * poly_degree,
                                              TestParametersSinEllipt, double> >
    HDG_wrapper(num_elements);

  vector<double> vectorRHS = HDG_wrapper.total_flux_vector(HDG_wrapper.return_zero_vector());
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)
    vectorRHS[i] *= -1.;

  vector<double> vectorSolution = SparseLA::conjugate_gradient(vectorRHS, HDG_wrapper);

  cout << "Error: " << HDG_wrapper.calculate_L2_error(vectorSolution) << endl;

  HDG_wrapper.plot_option("fileName", "diffusion_elliptic_c++");
  HDG_wrapper.plot_option("printFileNumber", "false");
  HDG_wrapper.plot_option("scale", "0.95");
  HDG_wrapper.plot_solution(vectorSolution);
}

int main()
{
  diffusion_elliptic();
  return 0;
}
