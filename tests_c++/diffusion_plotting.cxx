#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>

#include <string>

using namespace std;

/*!*************************************************************************************************
 * \brief   Set the test parameters.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TestParameters
{
  static constexpr std::array<unsigned int, 27U> dirichlet_nodes{
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                               const param_float_t time = 0.)
  {
    return 1.;
  }
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return 0.;
  }
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return point[0];
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& point,
                                     const param_float_t time = 0.)
  {
    return 0.;
  }
};
/*!*************************************************************************************************
 * \brief   Function that tests several aspects of the C++ implementation against a given reference
 *          solution obtained with the Python interface.
 *
 * This function implements an alternative to Executable.py (which usses the Cython interface).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dim, unsigned int space_dim, typename test_float_t>
int do_test()
{
  const std::vector<unsigned int> num_elements(space_dim, 3U);

  GlobalLoop::Elliptic<
    Topology::Cubic<hyEdge_dim, space_dim>, Geometry::UnitCube<hyEdge_dim, space_dim, test_float_t>,
    NodeDescriptor::Cubic<hyEdge_dim, space_dim>,
    LocalSolver::Diffusion<hyEdge_dim, 1, 2, TestParameters, test_float_t>, vector<test_float_t> >
    diffusion_problem(num_elements, (test_float_t)1.);

  std::vector<test_float_t> helper = diffusion_problem.return_zero_vector();

  std::string file_name = "diff_t-" + std::to_string(hyEdge_dim) + "-" + std::to_string(space_dim);
  std::string print_file_number = "false", scale = "0.8";
  std::string plot_edge_boundaries = "false", boundary_scale = "0.9";

  diffusion_problem.plot_option("fileName", file_name);
  diffusion_problem.plot_option("printFileNumber", print_file_number);
  diffusion_problem.plot_option("plotEdgeBoundaries", plot_edge_boundaries);
  diffusion_problem.plot_option("boundaryScale", boundary_scale);
  diffusion_problem.plot_option("scale", scale);

  file_name = "diff_tedge-" + std::to_string(hyEdge_dim) + "-" + std::to_string(space_dim);
  print_file_number = "false", scale = "0.8";
  plot_edge_boundaries = "true", boundary_scale = "0.9";
  std::string plot_edges = "false";
  diffusion_problem.plot_solution(helper);
  diffusion_problem.plot_option("fileName", file_name);
  diffusion_problem.plot_option("printFileNumber", print_file_number);
  diffusion_problem.plot_option("plotEdges", plot_edges);
  diffusion_problem.plot_option("plotEdgeBoundaries", plot_edge_boundaries);
  diffusion_problem.plot_option("boundaryScale", boundary_scale);
  diffusion_problem.plot_option("scale", scale);

  diffusion_problem.plot_solution(helper);

  return 0;
}

int main(int argc, char* argv[])
{
  int result = 0;

  result += do_test<1, 1, float>() + do_test<1, 2, float>() + do_test<1, 3, float>();
  result += do_test<2, 2, float>() + do_test<2, 3, float>() + do_test<3, 3, float>();

  return result;
}
