/*!*************************************************************************************************
 * \file    examples_c++/PointTest.C
 * \brief   File that tests several aspects of the C++ implementation of Point.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#include <Point.hxx>

#include <array>
#include <random>

using namespace std;

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
  
  const unsigned int space_dim = 5;
  const unsigned int array_len = 15;
  
//  array< Point<space_dim>, array
  
  return -1;
}
