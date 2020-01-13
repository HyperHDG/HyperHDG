/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2020
 */


#ifndef HYASSERT_H
#define HYASSERT_H

#include <iostream>

#ifndef NDEBUG
#define hy_assert(Expr, Msg) __Hy_Assert(#Expr, Expr, __FILE__, __LINE__, Msg)
#else
#define hy_assert(Expr, Msg) ;
#endif

void __Hy_Assert(const char* expr_str, bool expr, const char* file, int line, const char* msg)
{
  if (!expr)
  {
    std::cerr << "Assert failed:\t" << msg << "\n"
              << "Expected:\t" << expr_str << "\n"
              << "Source:\t\t" << file << ", line " << line << "\n";
    abort();
  }
}

#endif // end of ifndef HYASSERT_H
