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
#include <nanobind/ndarray.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <algorithm>
#include <cstddef>
#include <span>
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

// NumPy arrays across the boundary: inputs are read in place (any contiguous 1d array of the
// right dtype), outputs are allocated once here and owned by the returned array.
template <typename T>
using in_array = nanobind::ndarray<const T, nanobind::ndim<1>, nanobind::c_contig>;
template <typename T>
using out_array = nanobind::ndarray<nanobind::numpy, T, nanobind::ndim<1>>;
template <typename T>
using inout_array = nanobind::ndarray<T, nanobind::ndim<1>, nanobind::c_contig>;
template <typename T>
std::span<const T> as_span(const in_array<T>& a)
{
  return {a.data(), a.size()};
}
template <typename T>
std::span<T> as_span(inout_array<T>& a)
{
  return {a.data(), a.size()};
}
template <typename T>
out_array<T> new_array(const std::size_t n, std::span<T>& view)
{
  T* data = new T[n];
  nanobind::capsule owner(data, [](void* p) noexcept { delete[] static_cast<T*>(p); });
  view = std::span<T>(data, n);
  return out_array<T>(data, {n}, owner);
}
template <typename T>
out_array<T> to_array(const std::vector<T>& v)
{
  std::span<T> view;
  auto array = new_array<T>(v.size(), view);
  std::copy(v.begin(), v.end(), view.begin());
  return array;
}
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
      "read_dirichlet_indices", [](GlobalLoopT& o, const std::vector<unsigned int>& v)
      { o.read_dirichlet_indices(v); }, nb::arg("indices"));
  if constexpr (requires(GlobalLoopT o) { o.dirichlet_indices(); })
    cls.def("dirichlet_indices", [](GlobalLoopT& o) { return o.dirichlet_indices(); });

  // NumPy overloads of the vector methods (registered first: a list is not an ndarray and
  // falls through to the std::vector overloads below; an array never reaches the list caster).
  using InArr = detail::in_array<Scalar>;
  auto from_array = [](const InArr& a) { return Vec(a.data(), a.data() + a.size()); };
  if constexpr (requires(GlobalLoopT o, Vec v) { o.trace_to_flux(v, 0.); })
    cls.def(
      "trace_to_flux", [from_array](GlobalLoopT& o, const InArr& x, Scalar t)
      { return detail::to_array(o.trace_to_flux(from_array(x), t)); }, nb::arg("x"),
      nb::arg("time") = 0.);
  if constexpr (requires(GlobalLoopT o, Vec v) { o.residual_flux(v, 0.); })
    cls.def(
      "residual_flux", [from_array](GlobalLoopT& o, const InArr& x, Scalar t)
      { return detail::to_array(o.residual_flux(from_array(x), t)); }, nb::arg("x"),
      nb::arg("time") = 0.);
  if constexpr (requires(GlobalLoopT o, Vec v) { o.errors(v, 0.); })
    cls.def(
      "errors", [from_array](GlobalLoopT& o, const InArr& x, Scalar t)
      { return o.errors(from_array(x), t); }, nb::arg("x"), nb::arg("time") = 0.);
  // condensed system as NumPy triplets (rows, cols, values); factorize_once: uniform meshes only
  if constexpr (requires(GlobalLoopT o) { o.trace_to_flux_mat(0., true); })
    cls.def(
      "trace_to_flux_mat",
      [](GlobalLoopT& o, Scalar t, bool factorize_once)
      {
        auto mat = o.trace_to_flux_mat(t, factorize_once);
        return std::make_tuple(detail::to_array(mat.row_vec), detail::to_array(mat.col_vec),
                               detail::to_array(mat.value_vec));
      },
      nb::arg("time") = 0., nb::arg("factorize_once") = false);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.trace_to_flux(v, 0.); })
    cls.def(
      "trace_to_flux", [](GlobalLoopT& o, const Vec& v, Scalar t) { return o.trace_to_flux(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.trace_to_mass_flux(v, 0.); })
    cls.def(
      "trace_to_mass_flux", [](GlobalLoopT& o, const Vec& v, Scalar t)
      { return o.trace_to_mass_flux(v, t); }, nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.jacobian_of_trace_to_flux(v, 0., v, 0.); })
    cls.def(
      "jacobian_of_trace_to_flux",
      [](GlobalLoopT& o, const Vec& x, Scalar eig, const Vec& x_val, Scalar eig_val)
      { return o.jacobian_of_trace_to_flux(x, eig, x_val, eig_val); }, nb::arg("x"), nb::arg("eig"),
      nb::arg("x_val"), nb::arg("eig_val"));

  // condensed system as a COO triplet: assemble once, hand it to scipy, solve there
  if constexpr (requires(GlobalLoopT o) { o.trace_to_flux_mat(0.); })
    cls.def(
      "trace_to_flux_mat", [to_coo](GlobalLoopT& o, Scalar t)
      { return to_coo(o.trace_to_flux_mat(t)); }, nb::arg("time") = 0.);

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
      "plot_option", [](GlobalLoopT& o, const std::string& option, const std::string& value)
      { return o.plot_option(option, value); }, nb::arg("option"), nb::arg("value") = "");

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

  // ---- local-operator loop protocol (GlobalLoop::LocalOperatorLoop) ----------------------------
  // Like the C++ interface, nothing here allocates global vectors: apply, residual and assemble
  // write into caller-allocated NumPy arrays (n_triplets sizes the triplet arrays). Args is the
  // solver's runtime parameter struct, constructible from the time.
  if constexpr (requires {
                  typename GlobalLoopT::Args;
                  typename GlobalLoopT::index_type;
                })
  {
    using Args = typename GlobalLoopT::Args;
    using Idx = typename GlobalLoopT::index_type;
    using InIdx = detail::in_array<Idx>;
    using OutArr = detail::inout_array<Scalar>;
    using OutIdx = detail::inout_array<Idx>;
    if constexpr (requires(GlobalLoopT o, std::span<const Scalar> x, std::span<Scalar> y, Args a) {
                    o.apply(x, y, a);
                  })
      cls.def(
        "apply",
        [](GlobalLoopT& o, const InArr& x, OutArr& y, Scalar t)
        {
          auto out = detail::as_span(y);
          o.apply(detail::as_span(x), out, Args{t});
        },
        nb::arg("x"), nb::arg("y").noconvert(), nb::arg("time") = 0.);
    if constexpr (requires(GlobalLoopT o, std::span<const Scalar> x, std::span<Scalar> y, Args a) {
                    o.residual(x, y, a);
                  })
      cls.def(
        "residual",
        [](GlobalLoopT& o, const InArr& x, OutArr& r, Scalar t)
        {
          auto out = detail::as_span(r);
          o.residual(detail::as_span(x), out, Args{t});
        },
        nb::arg("x"), nb::arg("r").noconvert(), nb::arg("time") = 0.);
    if constexpr (requires(GlobalLoopT o, std::span<Idx> i, std::span<Scalar> v, Args a) {
                    o.assemble(a, true, i, i, v);
                  })
    {
      cls.def(
        "n_triplets", [](GlobalLoopT& o, bool identity_rows)
        { return o.n_triplets(identity_rows); }, nb::arg("identity_rows") = false);
      cls.def(
        "assemble",
        [](GlobalLoopT& o, OutIdx& rows, OutIdx& cols, OutArr& vals, Scalar t, bool identity_rows)
        {
          auto r = detail::as_span(rows);
          auto c = detail::as_span(cols);
          auto v = detail::as_span(vals);
          o.assemble(Args{t}, identity_rows, r, c, v);
        },
        nb::arg("rows").noconvert(), nb::arg("cols").noconvert(), nb::arg("vals").noconvert(),
        nb::arg("time") = 0., nb::arg("identity_rows") = false);
    }
    if constexpr (requires(GlobalLoopT o, std::span<const Scalar> x, Args a) { o.errors(x, a); })
      cls.def(
        "errors",
        [](GlobalLoopT& o, const InArr& x, Scalar t)
        {
          auto e = o.errors(detail::as_span(x), Args{t});
          return std::vector<Scalar>(e.data().begin(), e.data().end());
        },
        nb::arg("x"), nb::arg("time") = 0.);
    if constexpr (requires(GlobalLoopT o, Args a) { o.prepare(a); })
      cls.def(
        "prepare", [](GlobalLoopT& o, Scalar t) { o.prepare(Args{t}); }, nb::arg("time") = 0.);
    if constexpr (requires(GlobalLoopT o, std::span<const Idx> i) { o.set_constraints(i); })
      cls.def(
        "set_constraints", [](GlobalLoopT& o, const InIdx& i)
        { o.set_constraints(detail::as_span(i)); }, nb::arg("indices"));
    if constexpr (requires(GlobalLoopT o) { o.constraints(); })
      cls.def("constraints", [](GlobalLoopT& o) { return detail::to_array(o.constraints()); });
  }
  if constexpr (requires(GlobalLoopT o) { o.release(); })
    cls.def("release", [](GlobalLoopT& o) { o.release(); });
  if constexpr (requires(GlobalLoopT o) { o.prepared(); })
    cls.def("prepared", [](GlobalLoopT& o) { return o.prepared(); });
  if constexpr (requires(GlobalLoopT o) { o.n_operators(); })
    cls.def("n_operators", [](GlobalLoopT& o) { return o.n_operators(); });
  if constexpr (requires(GlobalLoopT o) { o.n_hyEdges(); })
    cls.def("n_hyEdges", [](GlobalLoopT& o) { return o.n_hyEdges(); });
  if constexpr (requires(GlobalLoopT o) { o.constrain_dirichlet_nodes(); })
    cls.def("constrain_dirichlet_nodes", [](GlobalLoopT& o) { o.constrain_dirichlet_nodes(); });

  return cls;
}

}  // end of namespace HyperHDG
