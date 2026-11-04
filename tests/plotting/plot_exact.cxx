#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>

#include <iostream>
#include <string>
#include <vector>

// VTU plots of a known solution for every (hyperedge dimension, space dimension) up to 3: every
// hypernode is Dirichlet with the value time * x, so the exact LDG-H solution of degree one is
// u = time * x, plotted at time 1. tests/plotting/check_plots.py reads the files through
// ParaView, checks u = x at every point and renders them; plot_vtkhdf does the same for VTKHDF
// (three time steps) and checks those files against HDF5 directly.

template <unsigned int space_dimT, typename param_float_t = double>
struct AllDirichlet
{
  static constexpr std::array<unsigned int, 27U> dirichlet_nodes{
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1.;
  }
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return time * point[0];
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <unsigned int hyEdge_dim, unsigned int space_dim>
void plot()
{
  GlobalLoop::Elliptic<
    Topology::Cubic<hyEdge_dim, space_dim>, Geometry::UnitCube<hyEdge_dim, space_dim, double>,
    NodeDescriptor::Cubic<hyEdge_dim, space_dim>,
    LocalSolver::Diffusion<hyEdge_dim, 1, 2, AllDirichlet, double>, std::vector<double> >
    problem(std::vector<unsigned int>(space_dim, 2U), 1.);
  const std::string name =
    "plot_exact-" + std::to_string(hyEdge_dim) + "-" + std::to_string(space_dim);
  problem.plot_option("fileName", name);
  problem.plot_option("printFileNumber", "false");
  problem.plot_solution(problem.zero_vector(), 1.);
  std::cout << "  output/" << name << ".vtu" << std::endl;
}

int main()
{
  plot<1, 1>();
  plot<1, 2>();
  plot<1, 3>();
  plot<2, 2>();
  plot<2, 3>();
  plot<3, 3>();
  return 0;
}
