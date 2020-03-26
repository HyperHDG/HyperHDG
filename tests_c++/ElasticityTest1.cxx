#include <HyperHDG/AbstractProblem.hxx>
#include <HyperHDG/Topology/File.hxx>
#include <HyperHDG/Geometry/File.hxx>
#include <LSol_Elasticity.hxx>
#include <HyperHDG/SparseLinearAlgebra.hxx>

#include <string>

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
int main(int argc, char *argv[])
{
  bool successful = true;
  const double solution_tolerance = 1e-8;
  
  std::string filename = "domains/SimpleTriangle.geo";
  AbstractProblem< Topology::File<1,2>, Geometry::File<1,2>, ElasticRods_TensorialUniform<1,2,1,2> >
    diffusion_problem(filename, 1.);

  vector<double> vectorDirichlet = diffusion_problem.return_zero_vector();
  vectorDirichlet[0] = 1.;
  vectorDirichlet[2] = 0.;
  
  const vector<unsigned int> index_vector = { 0 , 1 , 2 , 3 };
  diffusion_problem.read_dirichlet_indices(index_vector);
  
  vector<double> vectorRHS = diffusion_problem.matrix_vector_multiply(vectorDirichlet);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)  vectorRHS[i] *= -1.;
  
  vector<double> solution = conjugate_gradient( vectorRHS, diffusion_problem );
  solution = linear_combination(1., solution, 1., vectorDirichlet);
 
  const std::vector<double> python_result = { 1. , 0. , 0. , 0. , 0.5 , 0.5 };
  
  hy_assert ( solution.size() == python_result.size() ,
              "Size of solution of C++ program must be size of reference Python solution." );
  if ( solution.size() != python_result.size() )  successful = false;
  
  for (unsigned int i = 0; i < solution.size(); ++i)
  {
    hy_assert( abs( solution[i] - python_result[i] ) < solution_tolerance ,
               "Difference between Python's refrence solution ans the solution is too large, i.e. "
               << "it is " << abs( solution[i] - python_result[i] ) << " in the " << i << "-th " <<
               "component of the solution vector!" );
    if ( abs( solution[i] - python_result[i] ) >= solution_tolerance )  successful = false;
  }
  
  return successful - 1;
}
