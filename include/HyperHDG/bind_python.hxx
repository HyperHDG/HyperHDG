#pragma once

// Explicit nanobind bindings for a global-loop instantiation.
//
// Usage: write a tiny module .cxx that instantiates the global loop you actually need and
// registers it,
//
//   NB_MODULE(my_module, m) { HyperHDG::bind_python<MyGlobalLoop>(m, "MyGlobalLoop"); }
//
// then build it with nanobind_add_module (see python/CMakeLists.txt). This replaces the
// runtime cython compile machinery for workflows that want full control over the compile:
// the set of template instantiations is spelled out in the .cxx, CMake produces a plain .so,
// python imports it.
//
// Every method is guarded with `if constexpr (requires ...)`, so the helper binds whatever
// subset of the global-loop protocol the instantiation provides (elliptic loops get no
// stage/finalize entries, etc.). Vectors cross the boundary as std::vector copies -- simple
// and allocation-honest; use the sparse COO tuple + scipy on the python side for the system.

// nanobind (and thereby Python.h) must come before any header that pulls in libc feature-test
// macros, else pyconfig.h redefines _POSIX_C_SOURCE/_XOPEN_SOURCE
#include <nanobind/nanobind.h>
#include <nanobind/stl/complex.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/tuple.h>
#include <nanobind/stl/vector.h>

#include <HyperHDG/gauss_tableau.hxx>  // Gauss::StageTime

#include <span>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace HyperHDG
{

template <typename GlobalLoopT,
          typename TopoArgT = std::string,
          typename SolverArgT = std::vector<double>>
auto bind_python(nanobind::module_& m, const char* name)
{
  namespace nb = nanobind;
  using Scalar = typename GlobalLoopT::dof_value_t;
  using Vec = std::vector<Scalar>;
  using Idx = unsigned int;

  // COO triplet (rows, cols, values); duplicates are summed by scipy.sparse, matching the
  // PETSc MatSetValuesCOO semantics of the C++ drivers.
  using Coo = std::tuple<std::vector<Idx>, std::vector<Idx>, Vec>;
  auto to_coo = [](auto&& mat) -> Coo {
    mat.eliminate_zeros();
    return {std::move(mat.row_vec), std::move(mat.col_vec), std::move(mat.value_vec)};
  };

  auto cls = nb::class_<GlobalLoopT>(m, name);

  cls.def(nb::init<const TopoArgT&, const SolverArgT&>());

  if constexpr (requires(GlobalLoopT o) { o.size_of_system(); })
    cls.def("size_of_system", [](GlobalLoopT& o) { return o.size_of_system(); });
  if constexpr (requires(GlobalLoopT o) { o.n_local_dofs(); })
    cls.def("n_local_dofs", [](GlobalLoopT& o) { return o.n_local_dofs(); });
  if constexpr (requires(GlobalLoopT o) { o.n_owned_dofs(); })
    cls.def("n_owned_dofs", [](GlobalLoopT& o) { return o.n_owned_dofs(); });
  if constexpr (requires { GlobalLoopT::n_dofs_per_node; })
    cls.def_static("n_dofs_per_node", []() { return GlobalLoopT::n_dofs_per_node; });
  if constexpr (requires { GlobalLoopT::n_gauss_stages(); })
    cls.def_static("n_gauss_stages", []() { return GlobalLoopT::n_gauss_stages(); });
  if constexpr (requires { GlobalLoopT::n_gauss_reps(); })
    cls.def_static("n_gauss_reps", []() { return GlobalLoopT::n_gauss_reps(); });

  if constexpr (requires(GlobalLoopT o) { o.zero_vector(); })
    cls.def("zero_vector", [](GlobalLoopT& o) { return o.zero_vector(); });

  if constexpr (requires(GlobalLoopT o, Vec v) { o.make_initial(v, 0.); })
    cls.def(
      "make_initial", [](GlobalLoopT& o, const Vec& v, double t) { return o.make_initial(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o) { o.template trace_to_flux_mat<Idx, Vec>(0.); })
    cls.def(
      "trace_to_flux_mat",
      [to_coo](GlobalLoopT& o, double t) {
        return to_coo(o.template trace_to_flux_mat<Idx, Vec>(t));
      },
      nb::arg("time") = 0.);

  // Gauss stage entries: the stage index rides in a Gauss::StageTime through the loop's
  // generic entries, exactly as in the PETSc drivers (experiments/timowave.cxx).
  if constexpr (requires(GlobalLoopT o) {
                  o.template trace_to_flux_mat<Idx, Vec>(Gauss::StageTime{0., 0});
                })
    cls.def(
      "trace_to_flux_mat_stage",
      [to_coo](GlobalLoopT& o, Idx rep, double t) {
        return to_coo(
          o.template trace_to_flux_mat<Idx, Vec>(Gauss::StageTime{t, static_cast<int>(rep)}));
      },
      nb::arg("rep"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, std::span<Scalar> s) { o.residual_flux2(s, s, 0.); })
    cls.def(
      "residual_flux",
      [](GlobalLoopT& o, Vec x, double t) {
        Vec ax(x.size(), 0.);
        std::span<Scalar> x_span(x), ax_span(ax);
        o.residual_flux2(x_span, ax_span, t);
        return ax;
      },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, std::span<Scalar> s) {
                  o.residual_flux2(s, s, Gauss::StageTime{0., 0});
                })
    cls.def(
      "residual_flux_stage",
      [](GlobalLoopT& o, Vec x, Idx rep, double t) {
        Vec ax(x.size(), 0.);
        std::span<Scalar> x_span(x), ax_span(ax);
        o.residual_flux2(x_span, ax_span, Gauss::StageTime{t, static_cast<int>(rep)});
        return ax;
      },
      nb::arg("x"), nb::arg("rep"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.set_data(v, 0.); })
    cls.def(
      "set_data", [](GlobalLoopT& o, const Vec& v, double t) { o.set_data(v, t); }, nb::arg("x"),
      nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.set_data(v, Gauss::StageTime{0., 0}); })
    cls.def(
      "set_data_stage",
      [](GlobalLoopT& o, const Vec& v, Idx rep, double t) {
        o.set_data(v, Gauss::StageTime{t, static_cast<int>(rep)});
      },
      nb::arg("x"), nb::arg("rep"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o) { o.finalize_step(); })
    cls.def("finalize_step", [](GlobalLoopT& o) { o.finalize_step(); });

  if constexpr (requires(GlobalLoopT o, Vec v) { o.errors(v, 0.); })
    cls.def(
      "errors", [](GlobalLoopT& o, const Vec& v, double t) { return o.errors(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, Vec v) { o.norms(v, 0.); })
    cls.def(
      "norms", [](GlobalLoopT& o, const Vec& v, double t) { return o.norms(v, t); }, nb::arg("x"),
      nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o, std::string s) { o.plot_option(s, s); })
    cls.def(
      "plot_option",
      [](GlobalLoopT& o, const std::string& opt, const std::string& val) {
        return o.plot_option(opt, val);
      },
      nb::arg("option"), nb::arg("value") = "");

  if constexpr (requires(GlobalLoopT o, Vec v) { o.plot_solution(v, 0.); })
    cls.def(
      "plot_solution", [](GlobalLoopT& o, const Vec& v, double t) { o.plot_solution(v, t); },
      nb::arg("x"), nb::arg("time") = 0.);

  if constexpr (requires(GlobalLoopT o) { o.set_refinement(1u); })
    cls.def("set_refinement", [](GlobalLoopT& o, unsigned int n) { o.set_refinement(n); });
  if constexpr (requires(GlobalLoopT o) { o.get_refinement(); })
    cls.def("get_refinement", [](GlobalLoopT& o) { return o.get_refinement(); });

  return cls;
}

}  // end of namespace HyperHDG
