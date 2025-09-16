#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/sparse_la.hxx>

int main()
{
  using Real = double;
  using GLoop = GlobalLoop::Elliptic<
    Topology::File<1,3>,
    Geometry::File<1,3>,
    NodeDescriptor::File<1,3>,
    TimoshenkoBeam<1,3,5,10,LocalSolver::TimoshenkoBeamParametersClamped>
  >;

  GLoop hdg_wrapper(filepath);

  std::vector<Real> rhs = hdg_wrapper.residual_flux(hdg_wrapper.zero_vector());
  for (unsigned int i = 0; i < rhs.size(); ++i)
    rhs[i] *= -1.;

  std::vector<Real> sol = SparseLA::conjugate_gradient(rhs, hdg_wrapper, maxiters, tol);

  // std::cout << "Error: " << HDG_wrapper.errors(vectorSolution)[0] << std::endl;
  // hdg_wrapper.plot_option("fileName", "diffusion_elliptic_c++");
  // hdg_wrapper.plot_option("printFileNumber", "false");
  // hdg_wrapper.plot_option("scale", "0.95");
  // hdg_wrapper.plot_solution(vectorSolution);

  return 0;
}

