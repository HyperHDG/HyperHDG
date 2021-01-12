#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_uniform_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>

#include <type_traits>

using namespace std;
using namespace SparseLA;

/*!*************************************************************************************************
 * \brief   Function that tests several aspects of the C++ implementation against a given reference
 *          solution obtained with the Python interface.
 *
 * This function implements an alternative to Executable.py (which usses the Cython interface).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <typename test_float_t>
int do_test()
{
  test_float_t solution_tolerance = 0.;
  if constexpr (std::is_same<test_float_t, double>::value)
    solution_tolerance = 1e-8;
  else if constexpr (std::is_same<test_float_t, float>::value)
    solution_tolerance = 1e-6;

  bool successful = true;
  const vector<unsigned int> num_elements = {4, 2, 2};

  GlobalLoop::Elliptic<Topology::Cubic<1, 3>, Geometry::UnitCube<1, 3>, NodeDescriptor::Cubic<1, 3>,
                       LocalSolver::DiffusionUniform<1, 1, 2, test_float_t>,
                       std::vector<test_float_t> >
    diffusion_problem(num_elements, num_elements, (test_float_t)1.);

  vector<test_float_t> vectorDirichlet = diffusion_problem.zero_vector();
  vectorDirichlet[0] = (test_float_t)1.;
  vectorDirichlet[vectorDirichlet.size() - 1] = (test_float_t)0.;

  const vector<unsigned int> index_vector = {0, (unsigned int)vectorDirichlet.size() - 1};
  diffusion_problem.read_dirichlet_indices(index_vector);

  vector<test_float_t> vectorRHS = diffusion_problem.trace_to_flux(vectorDirichlet);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)
    vectorRHS[i] *= (test_float_t)-1.;

  vector<test_float_t> solution;
  try
  {
    solution = conjugate_gradient(vectorRHS, diffusion_problem);
  }
  catch (SolveException& exc)
  {
    hy_assert(0 == 1, exc.what());
    successful = false;
  }

  solution = linear_combination((test_float_t)1., solution, (test_float_t)1., vectorDirichlet);

  const std::vector<test_float_t> python_result = {
    1.,         0.6999695,  0.55280737, 0.46359316, 0.41591649, 0.72849089, 0.62353531, 0.52383342,
    0.4428244,  0.39207816, 0.63876017, 0.57986777, 0.5,        0.42013223, 0.36123983, 0.72849089,
    0.62353531, 0.52383342, 0.4428244,  0.39207816, 0.65166809, 0.58551499, 0.5,        0.41448501,
    0.34833191, 0.60792184, 0.5571756,  0.47616658, 0.37646469, 0.27150911, 0.63876017, 0.57986777,
    0.5,        0.42013223, 0.36123983, 0.60792184, 0.5571756,  0.47616658, 0.37646469, 0.27150911,
    0.58408351, 0.53640684, 0.44719263, 0.3000305,  0.};

  hy_assert(solution.size() == python_result.size(),
            "Size of solution of C++ program must be size of reference Python solution.");
  if (solution.size() != python_result.size())
    successful = false;

  for (unsigned int i = 0; i < solution.size(); ++i)
  {
    hy_assert(abs(solution[i] - python_result[i]) < solution_tolerance,
              "Difference between Python's refrence solution ans the solution is too large, i.e. "
                << "it is " << abs(solution[i] - python_result[i]) << " in the " << i << "-th "
                << "component of the solution vector!");
    if (abs(solution[i] - python_result[i]) >= solution_tolerance)
      successful = false;
  }

  return successful - 1;
}

int main()
{
  return do_test<float>() + do_test<double>();
}
