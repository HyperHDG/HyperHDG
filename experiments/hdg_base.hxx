#pragma once

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
  virtual Idx n_dofs_per_node() = 0;
  virtual Vector zero_vector() = 0;
  virtual Vector errors(const Span& x_vec, const Real time = 0.) = 0;
  virtual Vector norms(const Span& x_vec, const Real time = 0.) = 0;
  virtual Vector make_initial(const Vector& x_vec, const Real time = 0.) = 0;
  virtual sparse_mat<Vector> trace_to_flux_mat(const Real time = 0.) = 0;
  virtual void residual_flux2(Span x_vec, Span vec_Ax, Real time = 0.) = 0;
  virtual void set_data(Span x_vec, const Real time = 0.) = 0;
  virtual void set_refinement(unsigned int i) = 0;
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
  sparse_mat<Vector> trace_to_flux_mat(const Real time = 0.) {
    return hdg.trace_to_flux_mat(time);
  }
  void residual_flux2(Span x_vec, Span vec_Ax, Real time = 0.) {
    hdg.residual_flux2(x_vec, vec_Ax, time);
  }
  void set_data(Span x_vec, const Real time = 0.) {
    hdg.set_data(x_vec, time);
  }
  void set_refinement(unsigned int i) {
    hdg.set_refinement(i);
  }
};

