#include <vector>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/local_solver/timoshenko_network.hxx>

int main(int argc, char** argv)
{
  using Real = double;
  using HDG = GlobalLoop::Elliptic<
    Topology::File<1,3>,
    Geometry::File<1,3>,
    NodeDescriptor::File<1,3>,
    LocalSolver::TimoshenkoBeam<1,3,5,10,LocalSolver::TimoschenkoBeamParametersClamped>
  >;

  if (argc < 2)
    return 1;
  const char* filepath = argv[1];
  HDG hdg(filepath);

  std::vector<Real> rhs = hdg.residual_flux(hdg.zero_vector());
  for (unsigned int i = 0; i < rhs.size(); ++i)
    rhs[i] *= -1.;

  // std::vector<Real> sol = SparseLA::conjugate_gradient(rhs, hdg_wrapper, maxiters, tol);

  std::cout << "Error: " << hdg.errors(rhs)[0] << std::endl;
  // hdg_wrapper.plot_option("fileName", "diffusion_elliptic_c++");
  // hdg_wrapper.plot_option("printFileNumber", "false");
  // hdg_wrapper.plot_option("scale", "0.95");
  // hdg_wrapper.plot_solution(vectorSolution);

  return 0;
}

