#pragma once  // Ensure that file is included only once in a single compilation.

// #include <type_traits>

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
 * \param   func  The name of the function that is checked to be implemented
 * \param   name  The resulting struct whose value is true if the function is implemented
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


template<typename float_t>
static constexpr float_t heron_root(const float_t square)
{
  static_assert( square >= 0., "Each square of a number must be non-negative!");

  float_t root = 0.5 * (square + 1.);
  while( square / root < root - 100. * std::numeric_limits<float_t>::epsilon() )
    root = 0.5 * ( root + square / root );
  
  return root;
}
