/*!*************************************************************************************************
 * \file    examples_c++/DiffusionTest.C
 * \brief   File that tests several aspects of the C++ implementation against a given reference
 *          solution obtained with the Python interface.
 *
 * This file implements an alternative to Executable.py (which usses the Cython interface).
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#include <AbstractProblem.hxx>
#include <SparseLinearAlgebra.hxx>

#include <iostream>

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
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
int main(int argc, char *argv[])
{
  bool successful = true;
  const double solution_tolerance = 1e-8;
  const vector<int> num_elements = { 4 , 2 , 2 };
  
  DiffusionProblemRegularNaive<1,3,1> diffusion_problem(num_elements, num_elements, 1.);
  
  vector<double> vectorDirichlet = diffusion_problem.return_zero_vector();
  vectorDirichlet[0] = 1.;
  vectorDirichlet[vectorDirichlet.size()-1] = 0.;
  
  const vector<int> index_vector = { 0 , ((int) vectorDirichlet.size())-1 };
  diffusion_problem.read_dirichlet_indices(index_vector);
  
  vector<double> vectorRHS = diffusion_problem.matrix_vector_multiply(vectorDirichlet);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)  vectorRHS[i] *= -1.;
  
  vector<double> solution;
  try { solution = conjugate_gradient( vectorRHS, diffusion_problem ); }
  catch (SparseLASolveException exc)
  {
    hy_assert( 0 == 1 , exc.what() );
    successful = false;
  }
  
  solution = linear_combination(1., solution, 1., vectorDirichlet);
    
  const std::vector<double> python_result = 
  { 1.,         0.6999695,  0.55280737, 0.46359316, 0.41591649, 0.72849089,
    0.62353531, 0.52383342, 0.4428244,  0.39207816, 0.63876017, 0.57986777,
    0.5,        0.42013223, 0.36123983, 0.72849089, 0.62353531, 0.52383342,
    0.4428244,  0.39207816, 0.65166809, 0.58551499, 0.5,        0.41448501,
    0.34833191, 0.60792184, 0.5571756,  0.47616658, 0.37646469, 0.27150911,
    0.63876017, 0.57986777, 0.5,        0.42013223, 0.36123983, 0.60792184,
    0.5571756,  0.47616658, 0.37646469, 0.27150911, 0.58408351, 0.53640684,
    0.44719263, 0.3000305,  0. };
  
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
  
  if ( successful )
  {
    cout << "Diffusion test 1 was successful!" << endl;
    return 0;
  }
  else  cout << "Diffusion test 1 failed!" << endl;
  
  return -1;
}
