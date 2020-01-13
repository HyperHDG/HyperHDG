/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "HyAssert.h"
#include <iostream>
#include <string>
#include <sstream>


using namespace std;


void __Hy_Assert(const char* expr_str, bool expr, const char* file, int line, stringstream& msg)
{
  if (!expr)
  {
    const string str = msg.str();
    const char* argument = str.c_str();
    cerr << "Assert failed:\t" << argument << "\n"
         << "Expected:\t" << expr_str << "\n"
         << "Source:\t\t" << file << ", line " << line << "\n";
    abort();
  }
}
