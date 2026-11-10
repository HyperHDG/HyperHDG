#pragma once

// Explicit nanobind bindings for a global-loop instantiation.
//
// Usage: write a tiny module .cxx that instantiates the global loop you actually need and
// registers it,
//
//   NB_MODULE(my_module, m) { HyperHDG::bind_python<MyGlobalLoop>(m, "MyGlobalLoop"); }
//
// then build it with nanobind_add_module (python/hyperhdg.py generates such a CMake project and
// builds it on demand, cf. python/README.md). This is an alternative
// to the runtime cython compile machinery (import/ + cython/) for workflows that want full
// control over the compile: the set of template instantiations is spelled out in the .cxx,
// CMake produces a plain .so, python imports it.
//
// Every method is guarded with `if constexpr (requires ...)`, so the helper binds whatever
// subset of the global-loop protocol the instantiation provides (an elliptic loop gets no
// set_data/make_initial, an eigenvalue loop no residual_flux, ...) -- one helper for all
// loops instead of a .pyx/.pxd template pair per loop. Vectors cross the boundary as
// std::vector copies; the condensed system crosses as a COO triplet of NumPy arrays consumed by
// scipy.sparse.

// first include: nanobind pulls in Python.h, which must precede any header defining libc
// feature-test macros, else pyconfig.h redefines _POSIX_C_SOURCE/_XOPEN_SOURCE
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h>
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
/*!*************************************************************************************************
 * \brief   Whether the geometry of a global loop is constructed from (and refers to) a topology.
 *
 * The global loops are templates GlobalLoopT<TopologyT, GeometryT, NodeDescriptorT, LocalSolverT,
 * LargeVecT, dof_index_t>. Their constructor (topology, geometry, local solver) passes the second
 * argument to the constructor of GeometryT. The file based geometries are constructed from a
 * topology, a reference to whose domain information they keep, whereas the cubic geometries copy
 * the vector of elements. For other global loops, the answer is true (to be safe).
 **************************************************************************************************/
template <typename GlobalLoopT>
struct geometry_from_topology : std::true_type
{
};
template <template <class, class, class, class, typename, typename> class LoopT,
          class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename LargeVecT,
          typename dof_index_t>
struct geometry_from_topology<
  LoopT<TopologyT, GeometryT, NodeDescriptorT, LocalSolverT, LargeVecT, dof_index_t>>
: std::bool_constant<std::is_same_v<typename GeometryT::constructor_value_type, TopologyT>>
{
};
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

  // COO triplet (rows, cols, values) of NumPy arrays, ready for
  // scipy.sparse.coo_matrix((vals, (rows, cols))): the loop emits one entry per local matrix
  // position, duplicates are summed by scipy. The arrays take over the vectors (no copy).
  auto to_numpy = []<typename T>(std::vector<T>&& vec)
  {
    auto* owned = new std::vector<T>(std::move(vec));
    nb::capsule owner(owned, [](void* ptr) noexcept { delete static_cast<std::vector<T>*>(ptr); });
    return nb::ndarray<nb::numpy, T, nb::ndim<1>>(owned->data(), {owned->size()}, owner);
  };
  auto to_coo = [to_numpy](auto&& mat)
  {
    return std::make_tuple(to_numpy(std::move(mat.row_vec)), to_numpy(std::move(mat.col_vec)),
                           to_numpy(std::move(mat.value_vec)));
  };

  auto cls = nb::class_<GlobalLoopT>(m, name);

  // constructors: (topology), (topology, local solver), (topology, geometry, local solver)
  //
  // The last one is not bound if the geometry is constructed from a topology object, as the file
  // based geometries are, which keep a reference to its domain information: from python, that
  // topology would be a temporary converted from the file name, destroyed after the construction,
  // so that every later use of the geometry would read freed memory. The cubic geometries copy the
  // vector of elements, so their constructor is bound.
  constexpr bool bind_with_topology = std::is_constructible_v<GlobalLoopT, const TopoArgT&>;
  constexpr bool bind_with_solver =
    std::is_constructible_v<GlobalLoopT, const TopoArgT&, const SolverArgT&>;
  constexpr bool bind_with_geometry =
    std::is_constructible_v<GlobalLoopT, const TopoArgT&, const TopoArgT&, const SolverArgT&> &&
    !detail::geometry_from_topology<GlobalLoopT>::value;
  static_assert(bind_with_topology || bind_with_solver || bind_with_geometry,
                "No constructor of GlobalLoopT takes TopoArgT (and SolverArgT): pass the "
                "constructor argument types as template arguments of bind_python, e.g. "
                "std::vector<unsigned int> for the cubic topologies.");
  if constexpr (bind_with_topology)
    cls.def(nb::init<const TopoArgT&>(), nb::arg("topology"));
  if constexpr (bind_with_solver)
    cls.def(nb::init<const TopoArgT&, const SolverArgT&>(), nb::arg("topology"),
            nb::arg("local_solver"));
  if constexpr (bind_with_geometry)
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
  if constexpr (requires(GlobalLoopT o) { o.trace_to_flux_mat(0., true); })
    cls.def(
      "trace_to_flux_mat", [to_coo](GlobalLoopT& o, Scalar t, bool factorize_once)
      { return to_coo(o.trace_to_flux_mat(t, factorize_once)); }, nb::arg("time") = 0.,
      nb::arg("factorize_once") = false);
  else if constexpr (requires(GlobalLoopT o) { o.trace_to_flux_mat(0.); })
    cls.def(
      "trace_to_flux_mat", [to_coo](GlobalLoopT& o, Scalar t)
      { return to_coo(o.trace_to_flux_mat(t)); }, nb::arg("time") = 0.);

  // generalized eigenvalue problems in bulk and skeletal unknowns
  // (GlobalLoop::GeneralizedEigenvalue)
  if constexpr (requires(GlobalLoopT o) { o.stiffness_mat(); })
    cls.def("stiffness_mat", [to_coo](GlobalLoopT& o) { return to_coo(o.stiffness_mat()); });
  if constexpr (requires(GlobalLoopT o) { o.mass_mat(); })
    cls.def("mass_mat", [to_coo](GlobalLoopT& o) { return to_coo(o.mass_mat()); });
  if constexpr (requires(GlobalLoopT o) { o.n_skeleton_dofs(); })
    cls.def("n_skeleton_dofs", [](GlobalLoopT& o) { return o.n_skeleton_dofs(); });
  if constexpr (requires(GlobalLoopT o) { o.n_hyEdges(); })
    cls.def("n_hyEdges", [](GlobalLoopT& o) { return o.n_hyEdges(); });
  if constexpr (requires(GlobalLoopT o) { o.n_loc_bulk_dofs(); })
    cls.def("n_loc_bulk_dofs", [](GlobalLoopT& o) { return o.n_loc_bulk_dofs(); });
  if constexpr (requires(GlobalLoopT o) { o.coefficients(); })
    cls.def("coefficients", [to_numpy](GlobalLoopT& o) { return to_numpy(o.coefficients()); });
  if constexpr (requires(GlobalLoopT o, Vec v, std::vector<double> a) { o.bulk_values(v, a); })
    cls.def(
      "bulk_values", [to_numpy](GlobalLoopT& o, const Vec& v, const std::vector<double>& a)
      { return to_numpy(o.bulk_values(v, a)); }, nb::arg("x"), nb::arg("abscissas"));
  if constexpr (requires(GlobalLoopT o, std::vector<double> a) { o.bulk_points(a); })
    cls.def(
      "bulk_points", [to_numpy](GlobalLoopT& o, const std::vector<double>& a)
      { return to_numpy(o.bulk_points(a)); }, nb::arg("abscissas"));
  if constexpr (requires(GlobalLoopT o) { o.hyEdge_measures(); })
    cls.def("hyEdge_measures",
            [to_numpy](GlobalLoopT& o) { return to_numpy(o.hyEdge_measures()); });
  if constexpr (requires(GlobalLoopT o, Vec v, std::vector<double> a) { o.face_values(v, a); })
    cls.def(
      "face_values", [to_numpy](GlobalLoopT& o, const Vec& v, const std::vector<double>& a)
      { return to_numpy(o.face_values(v, a)); }, nb::arg("x"), nb::arg("abscissas"));
  if constexpr (requires(GlobalLoopT o, std::vector<double> a) { o.face_points(a); })
    cls.def(
      "face_points", [to_numpy](GlobalLoopT& o, const std::vector<double>& a)
      { return to_numpy(o.face_points(a)); }, nb::arg("abscissas"));
  if constexpr (requires(GlobalLoopT o) { o.face_measures(); })
    cls.def("face_measures", [to_numpy](GlobalLoopT& o) { return to_numpy(o.face_measures()); });
  if constexpr (requires(GlobalLoopT o) { o.face_node_types(); })
    cls.def("face_node_types",
            [to_numpy](GlobalLoopT& o) { return to_numpy(o.face_node_types()); });
  if constexpr (requires(GlobalLoopT o) { o.face_node_indices(); })
    cls.def("face_node_indices",
            [to_numpy](GlobalLoopT& o) { return to_numpy(o.face_node_indices()); });

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

  return cls;
}

}  // end of namespace HyperHDG
