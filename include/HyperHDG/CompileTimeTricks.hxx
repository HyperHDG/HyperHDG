#pragma once // Ensure that file is included only once in a single compilation.

#include <type_traits>

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct not_uses_geometry
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct not_uses_geometry<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().numerical_flux_from_lambda( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_total_flux
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_total_flux<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().numerical_flux_total( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_initial_flux
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_initial_flux<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().numerical_flux_initial( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_mass_multiply
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_mass_multiply<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().numerical_flux_from_mass( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_total_mass
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_total_mass<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().total_numerical_flux_mass( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_L2_error
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_L2_error<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().calc_L2_error_squared( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_L2_error_temp
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_L2_error_temp<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().calc_L2_error_squared_temp( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_lambda_values
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_lambda_values<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().lambda_values( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_set_data
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_set_data<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().set_data( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};

/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct generates a compile time error if the check is not correctly conducted!
 **************************************************************************************************/
template<typename, typename T>
struct has_flux_der
{
  static_assert( std::integral_constant<T, false>::value, // This is my way to write false ;)
                 "Second template parameter needs to be of function type." );
};
/*!*************************************************************************************************
 * \brief   Check if local solver uses geometry in numerical_flux_from_lambda.
 * 
 * This struct has value true if the function does not use a geometry!
 **************************************************************************************************/
template<typename C, typename Ret, typename... Args>
struct has_flux_der<C, Ret(Args...)>
{
  private:
    template<typename T>
    static constexpr auto check(T*)
    -> typename std::is_same
      < decltype( std::declval<T>().numerical_flux_der( std::declval<Args>()... ) ), Ret >
      ::type;
  
    template<typename>
    static constexpr std::false_type check(...);
    
    typedef decltype(check<C>(0)) type;
  
  public:
    static constexpr bool value = type::value;
};
