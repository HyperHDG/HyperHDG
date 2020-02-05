/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "HyAssert.hxx"
#include <iostream>
#include <sstream>


void __Hy_Assert(const char* expr_str, bool expr, const char* file, int line, std::stringstream& msg)
{
  if (!expr)
  {
    std::cerr << "Assert failed:  " << msg.str() << std::endl
              << "Expected:       " << expr_str  << std::endl
              << "Source:         " << file << ", line " << line << std::endl;
    abort();
  }
}
