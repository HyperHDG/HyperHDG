// LDG-H diffusion with GlobalLoop::Generic on one-dimensional hyperedges embedded in a space of
// dimension SPACE_DIM (set per executable by tests/CMakeLists.txt): the simplex of
// domains/simplex_1_<SPACE_DIM>.geo, points 0, e_1, ..., e_{SPACE_DIM}, every pair of points an
// edge. u = 0 at the origin, u = 1 at e_1 (Dirichlet), no source. The solution is linear on every
// edge, so polynomial degrees 1 and 2 must reproduce it up to round-off.
//
// Kirchhoff's law at a free node e_i (i > 1): its neighbours are the origin (edge length 1, u = 0),
// e_1 (length sqrt 2, u = 1) and the other free nodes (length sqrt 2, the same value by symmetry),
// so u(e_i) = c = 1 / (1 + sqrt 2). One function is linear on every edge and takes these values:
// u(x) = x_0 + c (x_1 + ... + x_{SPACE_DIM - 1}).
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <HyperHDG/topology/file.hxx>

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#ifndef SPACE_DIM
#define SPACE_DIM 3  // set by tests/CMakeLists.txt; the default serves editors and linters
#endif
constexpr unsigned int space_dim = SPACE_DIM;

template <unsigned int space_dimT, typename param_float_t = double>
struct Simplex
{
  static constexpr param_float_t c = 1. / (1. + M_SQRT2);
  // The node types of the simplex files are the node indices: the origin and e_1 are Dirichlet.
  static constexpr std::array<unsigned int, 2U> dirichlet_nodes{0, 1};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t = 0.)
  {
    param_float_t value = point[0];
    for (unsigned int d = 1; d < space_dimT; ++d)
      value += c * point[d];
    return value;
  }
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t = 0.)
  {
    return analytic_result(point);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <unsigned int degree>
int solve()
{
  using Loop =
    GlobalLoop::Generic<Topology::File<1, space_dim, std::vector, Point<space_dim, double> >,
                        Geometry::File<1, space_dim, std::vector, Point<space_dim, double> >,
                        NodeDescriptor::File<1, space_dim, std::vector, Point<space_dim, double> >,
                        LocalSolver::Diffusion<1, degree, 2 * degree, Simplex, double> >;
  Loop loop("domains/simplex_1_" + std::to_string(space_dim) + ".geo");
  loop.constrain_dirichlet_nodes();

  const unsigned int size = loop.size_of_system();
  std::vector<double> rhs(size), zero(size, 0.), x(size, 0.);
  loop.residual(zero, rhs);
  for (double& value : rhs)
    value = -value;
  SparseLA::SolverOptions options;
  options.rtol = 1e-14;
  const SparseLA::SolverResult cg =
    SparseLA::cg([&](const std::vector<double>& in, std::vector<double>& out)
                 { loop.apply(in, out); }, rhs, x, options);
  const double error = GlobalLoop::errors(loop, x)[0];
  const bool ok = cg.converged() && error < 1e-12;
  std::printf(
    "%s  space dimension %u, degree %u: %u edges, cg %s after %u iterations, error %.2e\n",
    ok ? "ok  " : "FAIL", space_dim, degree, (unsigned int)loop.n_hyEdges(),
    SparseLA::to_string(cg.reason), cg.iterations, error);
  return !ok;
}

int main()
{
  return solve<1>() + solve<2>();
}
