#pragma once  // Ensure that file is included only once in a single compilation.

#include <type_traits>

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 *
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template <typename, typename T>
struct not_uses_geometry
{
  static_assert(std::integral_constant<T, false>::value,  // This is my way to write false ;)
                "Second template parameter needs to be of function type.");
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 *
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template <typename C, typename Ret, typename... Args>
struct not_uses_geometry<C, Ret(Args...)>
{
 private:
  template <typename T>
  static constexpr auto check(T*) -> typename std::is_same<
    decltype(std::declval<T>().numerical_flux_from_lambda(std::declval<Args>()...)),
    Ret>::type;

  template <typename>
  static constexpr std::false_type check(...);

  typedef decltype(check<C>(0)) type;

 public:
  static constexpr bool value = type::value;
};
