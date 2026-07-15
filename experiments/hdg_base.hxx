#pragma once

#include <HyperHDG/gauss_tableau.hxx>  // Gauss::StageTime

#include <complex>
#include <vector>
#include <span>
#include <string>

struct HDGBase {
  using Real = double;
  using Idx = unsigned int;
  using Vector = std::vector<Real>;
  using Span = std::span<Real>;

  virtual void plot_solution(const Span& lambda, const Real time = 0.) = 0;
  virtual std::string plot_option(const std::string& option, std::string value = "") = 0;
  virtual Idx size_of_system() = 0;
  virtual Idx n_owned_dofs() = 0;
  virtual Idx n_local_dofs() = 0;
  virtual std::vector<Idx> local_to_global_dofs() = 0;
  virtual Idx n_space_dim() = 0;
  virtual std::vector<Real> owned_point_coords() = 0;
  virtual std::vector<Idx> owned_edges_global() = 0;
  virtual Idx n_dofs_per_node() = 0;
  virtual Vector zero_vector() = 0;
  virtual Vector errors(const Span& x_vec, const Real time = 0.) = 0;
  virtual Vector norms(const Span& x_vec, const Real time = 0.) = 0;
  virtual Vector make_initial(const Vector& x_vec, const Real time = 0.) = 0;
  virtual void make_initial_from_static(const Span& x_vec, const Real time = 0.) = 0;
  virtual sparse_mat<Vector> trace_to_flux_mat(const Real time = 0.) = 0;
  virtual void residual_flux2(Span x_vec, Span vec_Ax, Real time = 0.) = 0;
  virtual void set_data(Span x_vec, const Real time = 0.) = 0;
  virtual void finalize_step() = 0;
  virtual void set_refinement(unsigned int i) = 0;
  // Gauss stage interface (s >= 2, complex stage operators; defaults keep non-Gauss loops valid):
  virtual Idx n_gauss_stages() { return 1; }
  virtual Idx n_gauss_reps() { return 1; }
  virtual void stage_weights(Idx rep, Real& affine, Real& mult, Real& w_re, Real& w_im) {
    hy_check(false, "stage_weights is not available for this global loop");
  }
  virtual sparse_mat<std::vector<std::complex<Real>>> trace_to_flux_mat_stage(Idx stage,
                                                                              Real time = 0.) {
    hy_check(false, "trace_to_flux_mat_stage is not available for this global loop");
    return {};
  }
  virtual void residual_flux_stage(Span re_vec, Span im_vec, Idx stage, Real time = 0.) {
    hy_check(false, "residual_flux_stage is not available for this global loop");
  }
  virtual void set_data_stage(Span re_vec, Span im_vec, Idx stage, Real time = 0.) {
    hy_check(false, "set_data_stage is not available for this global loop");
  }
  virtual ~HDGBase() = default;
};

template<typename HDG>
struct HDGWrapper : HDGBase {
  HDG hdg;

  HDGWrapper(HDG&& h) : hdg(std::move(h)) {}

  void plot_solution(const Span& lambda, const Real time = 0.) {
    hdg.plot_solution(lambda, time);
  }
  std::string plot_option(const std::string& option, std::string value = "") {
    return hdg.plot_option(option, value);
  }
  Idx size_of_system() {
    return hdg.size_of_system();
  }
  Idx n_owned_dofs() {
    return hdg.n_owned_dofs();
  }
  Idx n_local_dofs() {
    return hdg.n_local_dofs();
  }
  std::vector<Idx> local_to_global_dofs() {
    return hdg.local_to_global_dofs();
  }
  Idx n_space_dim() {
    return HDG::space_dim();
  }
  std::vector<Real> owned_point_coords() {
    return hdg.owned_point_coords();
  }
  std::vector<Idx> owned_edges_global() {
    return hdg.owned_edges_global();
  }
  Idx n_dofs_per_node() {
    return HDG::n_dofs_per_node;
  }
  Vector zero_vector() {
    return hdg.zero_vector();
  }
  Vector errors(const Span& x_vec, const Real time = 0.) {
    return hdg.errors(x_vec, time);
  }
  Vector norms(const Span& x_vec, const Real time = 0.) {
    return hdg.norms(x_vec, time);
  }
  Vector make_initial(const Vector& x_vec, const Real time = 0.)
  {
    return hdg.make_initial(x_vec, time);
  }
  void make_initial_from_static(const Span& x_vec, const Real time = 0.)
  {
    return hdg.make_initial_from_static(x_vec, time);
  }
  sparse_mat<Vector> trace_to_flux_mat(const Real time = 0.) {
    return hdg.trace_to_flux_mat(time);
  }
  void residual_flux2(Span x_vec, Span vec_Ax, Real time = 0.) {
    hdg.residual_flux2(x_vec, vec_Ax, time);
  }
  void set_data(Span x_vec, const Real time = 0.) {
    hdg.set_data(x_vec, time);
  }
  void finalize_step() {
    // Gauss step completion; only the hyperbolic loop provides it (optional capability).
    if constexpr (requires { hdg.finalize_step(); })
      hdg.finalize_step();
    else
      hy_check(false, "finalize_step is not available for this global loop");
  }
  void set_refinement(unsigned int i) {
    hdg.set_refinement(i);
  }
  Idx n_gauss_stages() {
    if constexpr (requires { HDG::n_gauss_stages(); })
      return HDG::n_gauss_stages();
    else
      return 1;
  }
  Idx n_gauss_reps() {
    if constexpr (requires { HDG::n_gauss_reps(); })
      return HDG::n_gauss_reps();
    else
      return 1;
  }
  void stage_weights(Idx rep, Real& affine, Real& mult, Real& w_re, Real& w_im) {
    if constexpr (requires { hdg.stage_weights(rep, affine, mult, w_re, w_im); })
      hdg.stage_weights(rep, affine, mult, w_re, w_im);
    else
      hy_check(false, "stage_weights is not available for this global loop");
  }
  // The stage solves reuse the loop's generic entries, instantiated with complex vectors and a
  // Gauss::StageTime as the time argument; the re/im split here exists only because the PETSc
  // vectors of the driver are real. Instantiated only for multi-stage solvers -- the requires
  // guards of the loop cannot help here since the incompatibility (real-only local entries)
  // sits inside the generic loop bodies.
  using CVec = std::vector<std::complex<Real>>;
  using CSpan = std::span<std::complex<Real>>;
  static constexpr bool has_multi_stage() {
    if constexpr (requires { HDG::n_gauss_stages(); })
      return HDG::n_gauss_stages() > 1;
    else
      return false;
  }
  sparse_mat<CVec> trace_to_flux_mat_stage(Idx stage, Real time = 0.) {
    if constexpr (has_multi_stage())
      return hdg.template trace_to_flux_mat<Idx, CVec>(Gauss::StageTime{time, stage});
    else {
      hy_check(false, "trace_to_flux_mat_stage is not available for this global loop");
      return {};
    }
  }
  void residual_flux_stage(Span re_vec, Span im_vec, Idx stage, Real time = 0.) {
    if constexpr (has_multi_stage()) {
      CVec zin(re_vec.size(), 0.), zout(re_vec.size(), 0.);
      CSpan sin(zin), sout(zout);
      hdg.residual_flux2(sin, sout, Gauss::StageTime{time, stage});
      for (size_t k = 0; k < zout.size(); ++k) {
        re_vec[k] += zout[k].real();
        im_vec[k] += zout[k].imag();
      }
    } else
      hy_check(false, "residual_flux_stage is not available for this global loop");
  }
  void set_data_stage(Span re_vec, Span im_vec, Idx stage, Real time = 0.) {
    if constexpr (has_multi_stage()) {
      CVec z(re_vec.size());
      for (size_t k = 0; k < z.size(); ++k)
        z[k] = std::complex<Real>(re_vec[k], im_vec[k]);
      CSpan zs(z);
      hdg.set_data(zs, Gauss::StageTime{time, stage});
    } else
      hy_check(false, "set_data_stage is not available for this global loop");
  }
};

