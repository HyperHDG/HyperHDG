// Quadruple precision (std::float128_t) end to end: Diffusion on the unit square, degree 2, with
// an exact solution in the discrete space (u = x_0^2 + x_0 x_1 + x_1^2), which LDG-H reproduces
// exactly, so the error measures round-off only: about 1e-15 in double, about 1e-30 in
// std::float128_t. That needs the quadrature, the local factorization, CG and the error
// integration in the wide type; LAPACK has no std::float128_t, so the factorization runs on the
// generic LU of dense_la.hxx (namespace LU). Built only with -DHYPERHDG_QUAD_TESTS=ON and a
// compiler with std::float128_t (GCC >= 13).
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/cubic.hxx>

#include <array>
#include <cstdio>
#include <stdfloat>
#include <vector>

// u = sum_d x_d^2 + sum_{d < e} x_d x_e, -Laplace u = -2 dim, Dirichlet data everywhere
template <unsigned int space_dimT, typename param_float_t>
struct Quadratic
{
  static constexpr bool constant_coefficients = true;
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& x,
                                       const param_float_t = 0.)
  {
    param_float_t u = 0.;
    for (unsigned int d = 0; d < space_dimT; ++d)
      for (unsigned int e = d; e < space_dimT; ++e)
        u += x[d] * x[e];
    return u;
  }
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return -2. * space_dimT;
  }
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& x,
                                       const param_float_t = 0.)
  {
    return analytic_result(x);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <typename T>
double error(const double rtol)
{
  using Loop = GlobalLoop::Generic<Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, T>,
                                   NodeDescriptor::Cubic<2, 2>,
                                   LocalSolver::Diffusion<2, 2, 4, Quadratic, T>, T>;
  Loop loop(std::vector<unsigned int>{4, 4});
  loop.constrain_dirichlet_nodes();
  loop.prepare();
  const unsigned int size = loop.size_of_system();
  std::vector<T> rhs(size), zero(size, (T)0.), x(size, (T)0.);
  loop.residual(zero, rhs);
  for (T& value : rhs)
    value = -value;
  SparseLA::SolverOptions options;
  options.rtol = rtol;
  options.max_iterations = 10 * size;
  const SparseLA::SolverResult cg = SparseLA::cg(
    [&](const std::vector<T>& in, std::vector<T>& out) { loop.apply(in, out); }, rhs, x, options);
  const double e = (double)GlobalLoop::errors(loop, x)[0];
  std::printf("  %-12s cg %s after %u iterations, error %.2e\n",
              sizeof(T) == 16 ? "float128_t" : "double", SparseLA::to_string(cg.reason),
              cg.iterations, e);
  return cg.converged() ? e : 1.;
}

int main()
{
  const double e_double = error<double>(1e-14), e_quad = error<std::float128_t>(1e-30);
  const bool ok = e_double < 1e-12 && e_quad < 1e-27 && e_quad < 1e-10 * e_double;
  std::printf("%s  the exact solution is reproduced to round-off: double %.1e, float128_t %.1e\n",
              ok ? "ok  " : "FAIL", e_double, e_quad);
  return !ok;
}
