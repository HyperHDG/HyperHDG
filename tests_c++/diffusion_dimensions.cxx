#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/local_solver/diffusion_uniform_ldgh.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/sparse_la.hxx>

#include <cmath>
#include <string>
#include <type_traits>

using namespace std;
using namespace SparseLA;

/*!*************************************************************************************************
 * \brief   Function that tests several aspects of the C++ implementation against a given reference
 *          solution obtained with the Python interface.
 *
 * \todo    Should we also add naive tests like checking whether return_zero_vector() returns vector
 *          of correct size only containing zeros?
 *
 * This function implements an alternative to Executable.py (which usses the Cython interface).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <unsigned int space_dim, typename float_t>
int do_test_uniform()
{
  float_t solution_tolerance = 0.;
  if constexpr (std::is_same<float_t, double>::value)
    solution_tolerance = 1e-8;
  else if constexpr (std::is_same<float_t, float>::value)
    solution_tolerance = 1e-6;

  bool successful = true;
  const string file_name = "domains/simplex_1_" + to_string(space_dim) + ".geo";

  GlobalLoop::Elliptic<Topology::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       Geometry::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       NodeDescriptor::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       LocalSolver::DiffusionUniform<1, 1, 2, float_t>, std::vector<float_t> >
    diffusion_problem(file_name, (float_t)1.);

  vector<float_t> vectorDirichlet = diffusion_problem.template return_zero_vector();
  vectorDirichlet[0] = (float_t)0.;
  vectorDirichlet[1] = (float_t)1.;

  const vector<unsigned int> index_vector = {0, 1};
  diffusion_problem.read_dirichlet_indices(index_vector);

  vector<float_t> vectorRHS = diffusion_problem.matrix_vector_multiply(vectorDirichlet);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)
    vectorRHS[i] *= (float_t)-1.;

  vector<float_t> solution;
  try
  {
    solution = conjugate_gradient(vectorRHS, diffusion_problem);
  }
  catch (SolveException& exc)
  {
    hy_assert(0 == 1, exc.what());
    successful = false;
  }

  solution = linear_combination((float_t)1., solution, (float_t)1., vectorDirichlet);

  hy_assert(solution.size() == space_dim + 1,
            "Size of solution of C++ program must be size of reference Python solution.");
  if (solution.size() != space_dim + 1)
    successful = false;

  hy_assert(abs(solution[0] - 0.0) < solution_tolerance,
            "Difference between Python's refrence solution ans the solution is too large, i.e. "
              << "it is " << abs(solution[0] - 0.0) << " in the " << 0 << "-th "
              << "component of the solution vector!");
  if (abs(solution[0] - 0.0) >= solution_tolerance)
    successful = false;

  hy_assert(abs(solution[1] - 1.0) < solution_tolerance,
            "Difference between Python's refrence solution ans the solution is too large, i.e. "
              << "it is " << abs(solution[1] - 1.0) << " in the " << 1 << "-th "
              << "component of the solution vector!");
  if (abs(solution[1] - 1.0) >= solution_tolerance)
    successful = false;

  for (unsigned int i = 2; i < solution.size(); ++i)
  {
    hy_assert(abs(solution[i] - 0.5) < solution_tolerance,
              "Difference between Python's refrence solution ans the solution is too large, i.e. "
                << "it is " << abs(solution[i] - 0.5) << " in the " << i << "-th "
                << "component of the solution vector!");
    if (abs(solution[i] - 0.5) >= solution_tolerance)
      successful = false;
  }

  return successful - 1;
}
/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TestParameters
{
  static constexpr std::array<unsigned int, 0U> dirichlet_nodes{};
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
    return 0.;
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
 * \todo    Should we also add naive tests like checking whether return_zero_vector() returns vector
 *          of correct size only containing zeros?
 *
 * This function implements an alternative to Executable.py (which usses the Cython interface).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <unsigned int space_dim, typename float_t>
int do_test_standard()
{
  float_t solution_tolerance = 0.;
  if constexpr (std::is_same<float_t, double>::value)
    solution_tolerance = 1e-8;
  else if constexpr (std::is_same<float_t, float>::value)
    solution_tolerance = 1e-6;

  bool successful = true;
  const string file_name = "domains/simplex_1_" + to_string(space_dim) + ".geo";

  GlobalLoop::Elliptic<Topology::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       Geometry::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       NodeDescriptor::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       LocalSolver::Diffusion<1, 1, 2, TestParameters, float_t>, vector<float_t> >
    diffusion_problem(file_name, (float_t)1.);

  vector<float_t> vectorDirichlet = diffusion_problem.template return_zero_vector();
  vectorDirichlet[0] = (float_t)0.;
  vectorDirichlet[1] = (float_t)1.;

  const vector<unsigned int> index_vector = {0, 1};
  diffusion_problem.read_dirichlet_indices(index_vector);

  vector<float_t> vectorRHS = diffusion_problem.matrix_vector_multiply(vectorDirichlet);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)
    vectorRHS[i] *= (float_t)-1.;

  vector<float_t> solution;
  try
  {
    solution = conjugate_gradient(vectorRHS, diffusion_problem);
  }
  catch (SolveException& exc)
  {
    hy_assert(0 == 1, exc.what());
    successful = false;
  }

  solution = linear_combination((float_t)1., solution, (float_t)1., vectorDirichlet);

  hy_assert(solution.size() == space_dim + 1,
            "Size of solution of C++ program must be size of reference Python solution.");
  if (solution.size() != space_dim + 1)
    successful = false;

  hy_assert(abs(solution[0] - 0.0) < solution_tolerance,
            "Difference between Python's refrence solution ans the solution is too large, i.e. "
              << "it is " << abs(solution[0] - 0.0) << " in the " << 0 << "-th "
              << "component of the solution vector!");
  if (abs(solution[0] - 0.0) >= solution_tolerance)
    successful = false;

  hy_assert(abs(solution[1] - 1.0) < solution_tolerance,
            "Difference between Python's refrence solution ans the solution is too large, i.e. "
              << "it is " << abs(solution[1] - 1.0) << " in the " << 1 << "-th "
              << "component of the solution vector!");
  if (abs(solution[1] - 1.0) >= solution_tolerance)
    successful = false;

  float_t result = 1. / (1. + sqrt(2));

  for (unsigned int i = 2; i < solution.size(); ++i)
  {
    hy_assert(abs(solution[i] - result) < solution_tolerance,
              "Difference between Python's refrence solution ans the solution is too large, i.e. "
                << "it is " << abs(solution[i] - result) << " in the " << i << "-th "
                << "component of the solution vector!");
    if (abs(solution[i] - result) >= solution_tolerance)
      successful = false;
  }

  return successful - 1;
}
/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TestParametersDir
{
  static constexpr std::array<unsigned int, 2U> dirichlet_nodes{0, 1};
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
    return norm_infty(point);
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
 * \todo    Should we also add naive tests like checking whether return_zero_vector() returns vector
 *          of correct size only containing zeros?
 *
 * This function implements an alternative to Executable.py (which usses the Cython interface).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <unsigned int space_dim, typename float_t>
int do_test_standard_dir()
{
  float_t solution_tolerance = 0.;
  if constexpr (std::is_same<float_t, double>::value)
    solution_tolerance = 1e-8;
  else if constexpr (std::is_same<float_t, float>::value)
    solution_tolerance = 1e-6;

  bool successful = true;
  const string file_name = "domains/simplex_1_" + to_string(space_dim) + ".geo";

  GlobalLoop::Elliptic<Topology::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       Geometry::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       NodeDescriptor::File<1, space_dim, vector, Point<space_dim, float_t> >,
                       LocalSolver::Diffusion<1, 1, 2, TestParametersDir, float_t>,
                       vector<float_t> >
    diffusion_problem(file_name, (float_t)1.);

  vector<float_t> vectorRHS = diffusion_problem.template return_zero_vector();
  vectorRHS = diffusion_problem.total_flux_vector(vectorRHS);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)
    vectorRHS[i] *= (float_t)-1.;

  vector<float_t> solution;
  try
  {
    solution = conjugate_gradient(vectorRHS, diffusion_problem);
  }
  catch (SolveException& exc)
  {
    hy_assert(0 == 1, exc.what());
    successful = false;
  }

  hy_assert(solution.size() == space_dim + 1,
            "Size of solution of C++ program must be size of reference Python solution.");
  if (solution.size() != space_dim + 1)
    successful = false;

  hy_assert(abs(solution[0] - 0.0) < solution_tolerance,
            "Difference between Python's refrence solution ans the solution is too large, i.e. "
              << "it is " << abs(solution[0] - 0.0) << " in the " << 0 << "-th "
              << "component of the solution vector!");
  if (abs(solution[0] - 0.0) >= solution_tolerance)
    successful = false;

  hy_assert(abs(solution[1] - 0.0) < solution_tolerance,
            "Difference between Python's refrence solution ans the solution is too large, i.e. "
              << "it is " << abs(solution[1] - 0.0) << " in the " << 1 << "-th "
              << "component of the solution vector!");
  if (abs(solution[1] - 0.0) >= solution_tolerance)
    successful = false;

  float_t result = 1. / (1. + sqrt(2));

  for (unsigned int i = 2; i < solution.size(); ++i)
  {
    hy_assert(abs(solution[i] - result) < solution_tolerance,
              "Difference between Python's refrence solution ans the solution is too large, i.e. "
                << "it is " << abs(solution[i] - result) << " in the " << i << "-th "
                << "component of the solution vector!");
    if (abs(solution[i] - result) >= solution_tolerance)
      successful = false;
  }

  return successful - 1;
}

int main(int argc, char* argv[])
{
  int result = 0;

  result += do_test_uniform<1, float>() + do_test_uniform<1, double>();
  result += do_test_uniform<2, float>() + do_test_uniform<2, double>();
  result += do_test_uniform<3, float>() + do_test_uniform<3, double>();
  result += do_test_uniform<4, float>() + do_test_uniform<4, double>();
  result += do_test_uniform<5, float>() + do_test_uniform<5, double>();

  result += do_test_standard<1, float>() + do_test_standard<1, double>();
  result += do_test_standard<2, float>() + do_test_standard<2, double>();
  result += do_test_standard<3, float>() + do_test_standard<3, double>();
  result += do_test_standard<4, float>() + do_test_standard<4, double>();
  result += do_test_standard<5, float>() + do_test_standard<5, double>();

  result += do_test_standard_dir<1, float>() + do_test_standard_dir<1, double>();
  result += do_test_standard_dir<2, float>() + do_test_standard_dir<2, double>();
  result += do_test_standard_dir<3, float>() + do_test_standard_dir<3, double>();
  result += do_test_standard_dir<4, float>() + do_test_standard_dir<4, double>();
  result += do_test_standard_dir<5, float>() + do_test_standard_dir<5, double>();

  return result;
}
