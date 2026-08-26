#pragma once

// Explicit nanobind bindings for a global-loop instantiation.
//
// Usage: write a tiny module .cxx that instantiates the global loop you actually need and
// registers it,
//
//   NB_MODULE(my_module, m) { HyperHDG::bind_python<MyGlobalLoop>(m, "MyGlobalLoop"); }
//
// then build it with nanobind_add_module (see python/CMakeLists.txt). This is an alternative
// to the runtime cython compile machinery (import/ + cython/) for workflows that want full
// control over the compile: the set of template instantiations is spelled out in the .cxx,
// CMake produces a plain .so, python imports it.
//
// Every method is guarded with `if constexpr (requires ...)`, so the helper binds whatever
// subset of the global-loop protocol the instantiation provides (an elliptic loop gets no
// set_data/make_initial, an eigenvalue loop no residual_flux, ...) -- one helper for all
// loops instead of a .pyx/.pxd template pair per loop. Vectors cross the boundary as
// std::vector copies; the condensed system crosses as a COO triplet consumed by
// scipy.sparse.

// first include: nanobind pulls in Python.h, which must precede any header defining libc
// feature-test macros, else pyconfig.h redefines _POSIX_C_SOURCE/_XOPEN_SOURCE
#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

namespace HyperHDG
{

namespace detail
{
/*!*************************************************************************************************
 * \brief   Large vector type of a global loop: what zero_vector() returns, else the default.
 *
 * The loops keep \c LargeVecT private, but every loop that owns global vectors hands one out
 * through \c zero_vector(). Eigenvalue loops do not; for those the default matches the
 * default template argument of the loops themselves.
 **************************************************************************************************/
template <typename GlobalLoopT>
struct large_vec
{
  using type = std::vector<double>;
};
template <typename GlobalLoopT>
  requires requires(GlobalLoopT o) { o.zero_vector(); }
struct large_vec<GlobalLoopT>
{
  using type = std::remove_cvref_t<decltype(std::declval<GlobalLoopT&>().zero_vector())>;
};
template <typename GlobalLoopT>
using large_vec_t = typename large_vec<GlobalLoopT>::type;
}  // end of namespace detail

/*!*************************************************************************************************
 * \brief   Bind a global-loop instantiation as a python class.
 *
 * \tparam  GlobalLoopT   The instantiated global loop, e.g. GlobalLoop::Elliptic<...>.
 * \tparam  TopoArgT      Constructor argument of topology (and geometry): std::string for the
 *                        file based classes, std::vector<unsigned int> for the cubic ones.
 * \tparam  SolverArgT    Constructor argument of the local solver, e.g. the penalty tau.
 * \tparam  LargeVecT     Global vector type; deduced from zero_vector() where available.
 * \param   m             The nanobind module the class is added to.
 * \param   name          Name of the class on the python side.
 * \retval  cls           The nanobind class, so callers can add their own methods.
 **************************************************************************************************/
template <typename GlobalLoopT,
          typename TopoArgT = std::string,
          typename SolverArgT = double,
          typename LargeVecT = detail::large_vec_t<GlobalLoopT>>
auto bind_python(nanobind::module_& m, const char* name)
{
  namespace nb = nanobind;
  using Vec = LargeVecT;
  using Scalar = typename Vec::value_type;
  using Idx = unsigned int;

  // COO triplet (rows, cols, values), ready for scipy.sparse.coo_matrix((vals, (rows, cols))):
  // the loop emits one entry per local matrix position, duplicates are summed by scipy.
  using Coo = std::tuple<std::vector<Idx>, std::vector<Idx>, Vec>;
  auto to_coo = [](auto&& mat) -> Coo
  { return {std::move(mat.row_vec), std::move(mat.col_vec), std::move(mat.value_vec)}; };

  auto cls = nb::class_<GlobalLoopT>(m, name);

  // constructors: (topology), (topology, local solver), (topology, geometry, local solver)
  if constexpr (std::is_constructible_v<GlobalLoopT, const TopoArgT&>)
    cls.def(nb::init<const TopoArgT&>(), nb::arg("topology"));
  if constexpr (std::is_constructible_v<GlobalLoopT, const TopoArgT&, const SolverArgT&>)
    cls.def(nb::init<const TopoArgT&, const SolverArgT&>(), nb::arg("topology"),
            nb::arg("local_solver"));
  if constexpr (std::is_constructible_v<GlobalLoopT, const TopoArgT&, const TopoArgT&,
                                        const SolverArgT&>)
    cls.def(nb::init<const TopoArgT&, const TopoArgT&, const SolverArgT&>(), nb::arg("topology"),
            nb::arg("geometry"), nb::arg("local_solver"));

  if constexpr (requires(GlobalLoopT o) { o.size_of_system(); })
    cls.def("size_of_system", [](GlobalLoopT& o) { return o.size_of_system(); });

  if constexpr (requires(GlobalLoopT o) { o.zero_vector(); })
    cls.def("zero_vector", [](GlobalLoopT& o) { return o.zero_vector(); });

  if constexpr (requires(GlobalLoopT o, std::vector<unsigned int> v) {
                  o.read_dirichlet_indices(v);
                })
    cls.def(
      "read_dirichlet_indices",
      [](GlobalLoopT& o, const std::vector<unsigned int>& v) { o.read_dirichlet_indices(v); },
      nb::arg("indices"));
  if constexpr (requires(GlobalLoopT o) { o.dirichlet_indices(); })
    cls.def("dirichlet_indices", [](GlobalLoopT& o) { return o.dirichlet_indices(); });

  if constexpr (requires(GlobalLoopT o, Vec v) { o.trace_to_flux(v, 0.); })
    cls.def(
      "trace_to_flux", [](GlobalLoopT& o, const Vec& v, Scalar t) { return o.trace_to_flux(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.trace_to_mass_flux(v, 0.); })
    cls.def(
      "trace_to_mass_flux",
      [](GlobalLoopT& o, const Vec& v, Scalar t) { return o.trace_to_mass_flux(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.jacobian_of_trace_to_flux(v, 0., v, 0.); })
    cls.def(
      "jacobian_of_trace_to_flux",
      [](GlobalLoopT& o, const Vec& x, Scalar eig, const Vec& x_val, Scalar eig_val)
      { return o.jacobian_of_trace_to_flux(x, eig, x_val, eig_val); },
      nb::arg("x"), nb::arg("eig"), nb::arg("x_val"), nb::arg("eig_val"));

  // condensed system as a COO triplet: assemble once, hand it to scipy, solve there
  if constexpr (requires(GlobalLoopT o) { o.trace_to_flux_mat(0.); })
    cls.def(
      "trace_to_flux_mat", [to_coo](GlobalLoopT& o, Scalar t) { return to_coo(o.trace_to_flux_mat(t)); },
      nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.residual_flux(v, 0.); })
    cls.def(
      "residual_flux", [](GlobalLoopT& o, const Vec& v, Scalar t) { return o.residual_flux(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.set_data(v, 0.); })
    cls.def(
      "set_data", [](GlobalLoopT& o, const Vec& v, Scalar t) { o.set_data(v, t); }, nb::arg("x"),
      nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.make_initial(v, 0.); })
    cls.def(
      "make_initial", [](GlobalLoopT& o, const Vec& v, Scalar t) { return o.make_initial(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.errors(v, 0.); })
    cls.def(
      "errors", [](GlobalLoopT& o, const Vec& v, Scalar t) { return o.errors(v, t); }, nb::arg("x"),
      nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, std::string s) { o.plot_option(s, s); })
    cls.def(
      "plot_option",
      [](GlobalLoopT& o, const std::string& option, const std::string& value)
      { return o.plot_option(option, value); },
      nb::arg("option"), nb::arg("value") = "");

  if constexpr (requires(GlobalLoopT o, Vec v) { o.plot_solution(v, 0.); })
    cls.def(
      "plot_solution", [](GlobalLoopT& o, const Vec& v, Scalar t) { o.plot_solution(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o) { o.get_refinement(); })
    cls.def("get_refinement", [](GlobalLoopT& o) { return o.get_refinement(); });
  if constexpr (requires(GlobalLoopT o) { o.set_refinement(1u); })
    cls.def(
      "set_refinement", [](GlobalLoopT& o, unsigned int level) { o.set_refinement(level); },
      nb::arg("level"));

  return cls;
}

}  // end of namespace HyperHDG
