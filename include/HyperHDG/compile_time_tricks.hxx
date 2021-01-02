#pragma once  // Ensure that file is included only once in a single compilation.

#include <limits>
#include <type_traits>
#include <utility>

/*!*************************************************************************************************
 * \brief   Unused parametes will neither result in g++, nor in doxygen warnings if wrapped by this.
 **************************************************************************************************/
#define UNUSED(x) /* nothing */

/*!*************************************************************************************************
 * \brief   Check if some class implements some function with some signature.
 *
 * This macro receives the name of the function that is checked to be implemented with a given
 * signature (not handled to the macro itself) and the name of a struct (output) that can be used
 * to check whether function func is implemented.
 *
 * Having invoked the macro, we are able to use the template struct \c name to check whether
 * function \c fun is a (static or non-static) member function of an element of class \c C, where
 * \c Ret(Args) is the supposed signature.
 *
 * \param[in]   func    The name of the function that is checked to be implemented.
 * \param[out]  name    The resulting struct whose value is true if the function is implemented.
 **************************************************************************************************/
#define HAS_MEMBER_FUNCTION(func, name)                                                            \
  template <typename, typename T>                                                                  \
  struct name                                                                                      \
  {                                                                                                \
    static_assert(std::integral_constant<T, false>::value,                                         \
                  "Second template parameter must be function signature.");                        \
  };                                                                                               \
  template <typename C, typename Ret, typename... Args>                                            \
  struct name<C, Ret(Args...)>                                                                     \
  {                                                                                                \
   private:                                                                                        \
    template <typename T>                                                                          \
    static constexpr auto check(T*) ->                                                             \
      typename std::is_same<decltype(std::declval<T>().func(std::declval<Args>()...)), Ret>::type; \
    template <typename>                                                                            \
    static constexpr std::false_type check(...);                                                   \
    typedef decltype(check<C>(0)) type;                                                            \
                                                                                                   \
   public:                                                                                         \
    static constexpr bool value = type::value;                                                     \
  }

/*!*************************************************************************************************
 * \brief   Calculate the non-negative square-root of a non-negative number at compile time.
 *
 * \tparam      float_t The floating point type with respect to which the square root is evaluated.
 * \param[in]   square  The number whose square-root is evaluated.
 * \retval      root    The square root of the given number.
 **************************************************************************************************/
template <typename float_t>
static constexpr float_t heron_root(const float_t square)
{
  static_assert(square >= 0., "Each square of a number must be non-negative!");

  if constexpr (square == 0.)
    return 0.;

  constexpr float_t corr_fac = 1. - 10. * std::numeric_limits<float_t>::epsilon();
  float_t upper_root = 0.5 * (square + 1.);
  float_t lower_root = square / upper_root;
  while (lower_root < upper_root * corr_fac)
  {
    upper_root = 0.5 * (lower_root + upper_root);
    lower_root = square / upper_root;
  }

  return 0.5 * (lower_root + upper_root);
}
