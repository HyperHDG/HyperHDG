#include <iostream>

/*!*************************************************************************************************
 * \file    hy_assert.hxx
 * \brief   This file provides the function \c hy_assert.
 *
 * This is a wrapper file to provide a function that allows to use assertions that are similar to
 * those provided by cassert. That is, we define a macro \c hy_assert that implements assert.
 * If a user wants to use assertions, it is recommended to use \c hy_assert(\c Expr, \c Msg). The
 * use of the function \c __Hy_Assert is \b not recommended.
 *
 * Function \c hy_assert takes two arguments. The first argument is evaluated to a \c boolean and
 * if this returns \c true, nothing is done. If the argument is \c false, the running program is
 * terminated and the second argument is displayed as part of an error message. Here, the second
 * argument is handled as a \c stringstream (without initial \c <<). Thus, for two integers a and b,
 * a function call might look like: hy_assert( a == b , "Integers have not been the same, since a
 * turned out to be " << a << " and b was " << b << "." );
 *
 * Whether this functionality is active or not can be deduced via setting \c NDEBUG, when the code
 * is compiled. Using this functionality makes your program significantly slower. However, usage is
 * highly recommended for testing.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/

#pragma once  // Ensure that file is included only once in a single compilation.

#include <iostream>
#include <sstream>

#if __has_include(<execinfo.h>) && __has_include(<cxxabi.h>)
#include <cxxabi.h>
#include <execinfo.h>
#include <cstdlib>
#include <cstring>
#define HY_HAVE_BACKTRACE
#endif

/*!*************************************************************************************************
 * \brief   Best-effort stack trace to stderr, printed on hy_check / hy_assert failure.
 *
 * glibc backtrace + demangling; function names of the executable's own frames need exported
 * symbols (CMAKE_ENABLE_EXPORTS / -rdynamic, set in the top-level CMakeLists) -- without them
 * (and inside hidden-visibility shared objects like the python modules) frames degrade to raw
 * addresses, which `addr2line -e <binary>` still resolves.
 **************************************************************************************************/
inline void __hy_print_stacktrace()
{
#ifdef HY_HAVE_BACKTRACE
  void* frames[64];
  const int n_frames = backtrace(frames, 64);
  char** symbols = backtrace_symbols(frames, n_frames);
  if (!symbols)
    return;
  std::cerr << "Stack trace (innermost first):" << std::endl;
  for (int i = 1; i < n_frames; ++i)  // frame 0 is this function
  {
    // symbols[i] reads "module(mangled+0xoffset) [addr]"; demangle the middle when present
    char* begin = std::strchr(symbols[i], '(');
    char* plus = begin ? std::strchr(begin, '+') : nullptr;
    char* demangled = nullptr;
    if (begin && plus && plus > begin + 1)
    {
      *plus = '\0';
      int status = 0;
      demangled = abi::__cxa_demangle(begin + 1, nullptr, nullptr, &status);
      *plus = '+';
    }
    std::cerr << "  #" << i - 1 << "  " << (demangled ? demangled : symbols[i]) << std::endl;
    std::free(demangled);
  }
  std::free(symbols);
#endif
}

#define hy_check(Expr, Msg)                                                 \
  do {                                                                       \
    if (!(Expr)) {                                                           \
      std::stringstream __hy_check_text;                                     \
      __hy_check_text << Msg;                                                \
      std::cerr << "Check failed: " << #Expr                                 \
                << "\n  at " << __FILE__ << ":" << __LINE__                  \
                << "\n  " << __hy_check_text.str() << std::endl;             \
      __hy_print_stacktrace();                                               \
      std::abort();                                                          \
    }                                                                        \
  } while (0)

#ifndef NDEBUG

#include <iostream>
#include <sstream>

/*!*************************************************************************************************
 * \brief   The assertion to be used within HyperHDG --- deactivate using -DNDEBUG compile flag.
 *
 * \param   Expr  C++ Expression that can be evaluated to \c true or \c false.
 * \param   Msg   Message that is to be displayed if \c Expr is evaluated to \c false.
 **************************************************************************************************/
// The message stringstream is only constructed on the failure path: building and formatting it
// unconditionally made assert-enabled builds ~400x slower in the per-edge hot loops (ne9-14),
// while the checks themselves are nearly free.
#define hy_assert(Expr, Msg)                                              \
  {                                                                       \
    if (!(Expr)) [[unlikely]]                                             \
    {                                                                     \
      std::stringstream __hy_assertion_text;                              \
      __hy_assertion_text << Msg;                                         \
      __Hy_Assert(#Expr, false, __FILE__, __LINE__, __hy_assertion_text); \
    }                                                                     \
  }                                                                       \
  static_assert(true, "")

// -------------------------------------------------------------------------------------------------
/// \cond EXCLUDE_CODE
// -------------------------------------------------------------------------------------------------

/*!*************************************************************************************************
 * \brief   This function is not (never) to be used.
 *
 * This function is \b not to be used in regular code. It only / solely is defined to allow the use
 * of function \c hy_assert( \c Expr, \c Msg) which is implemented as a macro in file HyAssert.hxx.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
inline void __Hy_Assert(const char* expr_str,
                        bool expr,
                        const char* file,
                        int line,
                        std::stringstream& msg)
{
  if (!expr)
  {
    std::cerr << "Assert failed:  " << msg.str() << std::endl
              << "Expected:       " << expr_str << std::endl
              << "Source:         " << file << ", line " << line << std::endl;
    __hy_print_stacktrace();
    abort();
  }
}

#else  // alternative branch of ifndef NDEBUG
#define hy_assert(Expr, Msg) \
  {                          \
    ;                        \
  }
#endif  // end of ifndef NDEBUG

// -------------------------------------------------------------------------------------------------
/// \endcond
// -------------------------------------------------------------------------------------------------
