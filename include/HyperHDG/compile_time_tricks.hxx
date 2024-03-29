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
