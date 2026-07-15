#pragma once  // Ensure that file is included only once in a single compilation.

#include <limits>
#include <type_traits>
#include <utility>

/*!*************************************************************************************************
 * \brief   Unused parametes will neither result in g++, nor in doxygen warnings if wrapped by this.
 **************************************************************************************************/
#define UNUSED(x) /* nothing */

/*!*************************************************************************************************
 * \brief   Dependent false for static_assert in discarded if-constexpr branches.
 *
 * Local-solver dispatch tests candidate call expressions with C++20 requires-expressions directly
 * at the call site (cf. global_loop/prototype.hxx). The fallback branch uses
 * \c static_assert(always_false_v<...>, ...) so that a local solver implementing no usable
 * overload is a compile error instead of a Release-silent \c hy_assert.
 **************************************************************************************************/
template <typename...>
inline constexpr bool always_false_v = false;
/*!*************************************************************************************************
 * \brief   declval-style helper: yields a prvalue copy of its argument (unevaluated contexts only).
 *
 * Inside the dispatch requires-expressions the trailing scalar (time/eigenvalue) must be passed as
 * an rvalue: local solvers commonly default that parameter (e.g. \c time \c = \c 0.) behind a
 * \c hyEdgeT& parameter, and a named lvalue would bind to that reference, making the shorter
 * overload spuriously viable (its body then fails to compile). Not defined — like std::declval.
 **************************************************************************************************/
template <typename T>
std::remove_cvref_t<T> prvalue_of(T&&) noexcept;
