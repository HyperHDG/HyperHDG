#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_uniform_ldgh.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/sparse_la.hxx>

#include <cmath>
#include <string>
#include <type_traits>

using namespace std;
using namespace SparseLA;

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

  vector<float_t> vectorDirichlet = diffusion_problem.zero_vector();
  vectorDirichlet[0] = (float_t)0.;
  vectorDirichlet[1] = (float_t)1.;

  const vector<unsigned int> index_vector = {0, 1};
  diffusion_problem.read_dirichlet_indices(index_vector);

  vector<float_t> vectorRHS = diffusion_problem.trace_to_flux(vectorDirichlet);
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

int main()
{
  int result = 0;

  result += do_test_uniform<1, float>() + do_test_uniform<1, double>();
  result += do_test_uniform<2, float>() + do_test_uniform<2, double>();
  result += do_test_uniform<3, float>() + do_test_uniform<3, double>();
  result += do_test_uniform<4, float>() + do_test_uniform<4, double>();
  result += do_test_uniform<5, float>() + do_test_uniform<5, double>();

  return result;
}
