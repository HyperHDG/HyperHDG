#include <HyperHDG/AbstractProblem.hxx>
#include <HyperHDG/SparseLinearAlgebra.hxx>
#include <HyperHDG/Geometry/Cubic.hxx>
#include <HyperHDG/NodeDescriptor/Cubic.hxx>
#include <HyperHDG/LocalSolver/Diffusion.hxx>

#include <cmath>
#include <string>
#include <type_traits>

using namespace std;
using namespace SparseLA;


/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template < unsigned int space_dimT, typename param_float_t = double >
struct TestParameters
{
  static constexpr double pi = acos(-1);
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes
  { 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26 };
  static constexpr std::array<unsigned int, 0U> neumann_nodes {};
  static param_float_t inverse_diffusion_coeff( const Point<space_dimT,param_float_t>& pt )
  { return pi; }
  static param_float_t analytic_result( const Point<space_dimT,param_float_t>& pt )
  { return sin(pi * pt[0]); }
  static param_float_t right_hand_side( const Point<space_dimT,param_float_t>& pt )
  { return pi * sin(pi * pt[0]); }
  static param_float_t dirichlet_value( const Point<space_dimT,param_float_t>& pt )
  { return analytic_result(pt); }
  static param_float_t neumann_value( const Point<space_dimT,param_float_t>& pt )
  { return 0.; }
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
template < unsigned int hyEdge_dim, unsigned int space_dim, typename float_t >
double do_test(const unsigned int iteration)
{
  const std::vector<unsigned int> num_elements(space_dim, (unsigned int) (1 << iteration));

  AbstractProblem
  < 
    Topology::Cubic<hyEdge_dim,space_dim>, Geometry::UnitCube<hyEdge_dim,space_dim,float_t>, 
    NodeDescriptor::Cubic<hyEdge_dim,space_dim>,
    Diffusion<hyEdge_dim,1,2,TestParameters,float_t>
  >  diffusion_problem(num_elements, (float_t) 10.);
  
  vector<float_t> vectorRHS = diffusion_problem.template return_zero_vector<float_t>();
  vectorRHS = diffusion_problem.total_flux_vector(vectorRHS);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)  vectorRHS[i] *= (float_t) -1.;
  
  vector<float_t> solution;
  try { solution = conjugate_gradient( vectorRHS, diffusion_problem ); }
  catch (SparseLASolveException& exc) { hy_assert( 0 == 1 , exc.what() ); }

  std::string file_name = "diff_c-" + std::to_string(hyEdge_dim) + "-" + std::to_string(space_dim);
  std::string print_file_number = "false" , scale = "0.95";

  diffusion_problem.plot_option( "fileName" , file_name );
  diffusion_problem.plot_option( "printFileNumber" , print_file_number );
  diffusion_problem.plot_option( "scale" , scale );

  diffusion_problem.plot_solution(solution);
  
  return diffusion_problem.calculate_L2_error(solution);
}


int main(int argc, char *argv[])
{ 
  for (unsigned int i = 1; i < 8; ++i)
    std::cout << do_test<1,1,double>(i) << std::endl;
  std::cout << std::endl;
  for (unsigned int i = 1; i < 8; ++i)
    std::cout << do_test<1,2,double>(i) << std::endl;
  std::cout << std::endl;
  for (unsigned int i = 1; i < 8; ++i)
    std::cout << do_test<2,2,double>(i) << std::endl;
  return 0;
}