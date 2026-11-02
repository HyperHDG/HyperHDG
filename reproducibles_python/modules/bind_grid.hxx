#pragma once
/*!*************************************************************************************************
 * \brief   Bind a (dimension, degree) grid of global-loop instantiations as python classes.
 *
 * The modules of the reproducibles (modules/*.cxx) each define an alias template
 * `template <unsigned dim, unsigned degree> using X = GlobalLoop::...<...>` and call
 * `bind_grid<X>(m, "X", seq<1, 2>{}, seq<1, 2, 3>{})`, which binds the classes
 * `X_D1_P1`, ..., `X_D2_P3`. Extending a module to another dimension or degree is one number in
 * that call; a class the scripts ask for but no module provides is a KeyError at runtime.
 **************************************************************************************************/
#include <HyperHDG/bind_python.hxx>

#include <deque>
#include <string>
#include <utility>
#include <vector>

namespace reproducibles
{
template <unsigned... values>
using seq = std::integer_sequence<unsigned, values...>;

/*!*************************************************************************************************
 * \brief   Name of the class of one instantiation; nanobind keeps the pointer, so the string
 *          lives as long as the module.
 **************************************************************************************************/
inline const char* class_name(const char* prefix, unsigned dim, unsigned degree)
{
  static std::deque<std::string> names;
  names.push_back(std::string(prefix) + "_D" + std::to_string(dim) + "_P" + std::to_string(degree));
  return names.back().c_str();
}

/*!*************************************************************************************************
 * \brief   Bind Loop<dim, degree> as "<prefix>_D<dim>_P<degree>" for every dim and degree.
 *
 * \tparam  Loop        Alias template (dim, degree) -> instantiated global loop.
 * \tparam  TopoArgT    Constructor argument of the topology: std::vector<unsigned int> for the
 *                      cubic classes (default), std::string for the file based ones.
 * \tparam  SolverArgT  Constructor argument of the local solver (the penalty, or the vector of
 *                      constants of the parabolic solvers).
 **************************************************************************************************/
template <template <unsigned, unsigned> class Loop,
          typename TopoArgT = std::vector<unsigned int>,
          typename SolverArgT = double,
          unsigned... dims,
          unsigned... degrees>
void bind_grid(nanobind::module_& m,
               const char* prefix,
               std::integer_sequence<unsigned, dims...>,
               std::integer_sequence<unsigned, degrees...>)
{
  (
    [&]<unsigned dim>()
    {
      (HyperHDG::bind_python<Loop<dim, degrees>, TopoArgT, SolverArgT>(
         m, class_name(prefix, dim, degrees)),
       ...);
    }.template operator()<dims>(),
    ...);
}
}  // end of namespace reproducibles
