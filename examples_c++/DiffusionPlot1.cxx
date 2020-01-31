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

#include "../AbstractProblem.hxx"
#include "../SparseLinearAlgebra.hxx"
#include "../HyAssert.hxx"
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
  
  DiffusionProblemRegularNaive<2,3,1> diffusion_problem(num_elements, num_elements, 1.);  
  vector<double> vectorDirichlet = diffusion_problem.return_zero_vector();

  diffusion_problem.plot_solution(vectorDirichlet);
}
