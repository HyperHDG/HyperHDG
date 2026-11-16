#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/global_loop/kinds.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/local_solver/local_operator.hxx>

#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

#include <array>
#include <cmath>
#include <limits>
#include <type_traits>

namespace LocalSolver
{
namespace detail
{
/*!*************************************************************************************************
 * \brief   Eigendecomposition A = V diag(lambda) V^T of a small symmetric matrix (cyclic Jacobi).
 *
 * Row-major n x n; V orthogonal, its columns the eigenvectors. Free of LAPACK, so that it works
 * for any floating point type.
 **************************************************************************************************/
template <unsigned int n, typename T>
void symmetric_eigen(std::array<T, n * n> a, std::array<T, n * n>& v, std::array<T, n>& lambda)
{
  v.fill((T)0.);
  for (unsigned int i = 0; i < n; ++i)
    v[i * n + i] = (T)1.;
  T scale = (T)0.;
  for (const T& entry : a)
    scale = std::max(scale, (T)std::abs(entry));
  for (unsigned int sweep = 0; sweep < 64; ++sweep)
  {
    T off = (T)0.;
    for (unsigned int p = 0; p < n; ++p)
      for (unsigned int q = p + 1; q < n; ++q)
        off = std::max(off, (T)std::abs(a[p * n + q]));
    if (off <= std::numeric_limits<T>::epsilon() * scale)
      break;
    for (unsigned int p = 0; p < n; ++p)
      for (unsigned int q = p + 1; q < n; ++q)
      {
        if (a[p * n + q] == (T)0.)
          continue;
        // rotation that annihilates a(p, q)
        const T theta = (a[q * n + q] - a[p * n + p]) / ((T)2. * a[p * n + q]);
        const T t = (theta >= (T)0. ? (T)1. : (T)-1.) /
                    ((T)std::abs(theta) + (T)std::sqrt(theta * theta + (T)1.));
        const T c = (T)1. / (T)std::sqrt(t * t + (T)1.), s = t * c;
        for (unsigned int k = 0; k < n; ++k)
        {
          const T akp = a[k * n + p], akq = a[k * n + q];
          a[k * n + p] = c * akp - s * akq;
          a[k * n + q] = s * akp + c * akq;
        }
        for (unsigned int k = 0; k < n; ++k)
        {
          const T apk = a[p * n + k], aqk = a[q * n + k];
          a[p * n + k] = c * apk - s * aqk;
          a[q * n + k] = s * apk + c * aqk;
        }
        for (unsigned int k = 0; k < n; ++k)
        {
          const T vkp = v[k * n + p], vkq = v[k * n + q];
          v[k * n + p] = c * vkp - s * vkq;
          v[k * n + q] = s * vkp + c * vkq;
        }
      }
  }
  for (unsigned int i = 0; i < n; ++i)
    lambda[i] = a[i * n + i];
}
}  // end of namespace detail

/*!*************************************************************************************************
 * \brief   LDG-H diffusion with a sum-factorized, matrix-free local operator (fast
 *diagonalization).
 *
 * The same discretization as \c Diffusion (same shape functions, quadrature and unknowns, so the
 * same discrete solution), for hyperedges that are boxes (mat_r diagonal, axis-aligned local
 * normals, as UnitCube and the one-dimensional hyperedges of File) and a diffusion coefficient that
 * is constant on each hyperedge (the parameters declare constant_coefficients). Then the local
 * problem splits into one-dimensional matrices:
 *   - the flux unknowns q_d have the mass k |E| I (k the inverse diffusion coefficient, the
 *     Legendre basis orthonormal) and are eliminated, q_d = (b_q + G_d u) / (k |E|);
 *   - the Schur complement in u is S_u = sum_d X_d acting in direction d (identity in the others),
 *     X_d = tau |F_d| sum_s b_s b_s^T + H_d G_d / (k |E|) with the 1D derivative matrix G_d, its
 *     partner H_d = -G_d + boundary terms (= G_d^T up to quadrature) and the endpoint values b_s;
 *   - S_u^-1 is applied by fast diagonalization: X_d = V_d Lambda_d V_d^T, transform with V_d^T
 *     direction by direction, divide by sum_d Lambda_d, transform back with V_d;
 *   - B and C (trace to volume and back) are extensions and restrictions along one direction.
 * One application costs about 4 d (p+1)^(d+1) + O(d (p+1)^d) operations and the operator stores
 * O(d (p+1)^2) numbers, against (2 d (p+1)^(d-1))^2 for the condensed matrix of Diffusion; nothing
 * is factorized, and no LAPACK is needed.
 *
 * The load, the errors and the Dirichlet nodes are Diffusion's (composition: a Diffusion member).
 * The local operator implements apply, residual and local_solution for GlobalLoop::Generic; the
 * legacy interface of the other global loops is not provided.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT = DiffusionParametersDefault,
          typename lSol_float_t = double>
class DiffusionTensor
{
  using T = lSol_float_t;
  using Base = Diffusion<hyEdge_dimT, poly_deg, quad_deg, parametersT, T, LocalOperator::Factored>;
  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT> >,
    T>
    integrator;

  static constexpr unsigned int dim = hyEdge_dimT, n1 = poly_deg + 1;
  static constexpr unsigned int n_vol = Hypercube<dim>::pow(n1);       // shape functions
  static constexpr unsigned int n_face = Hypercube<dim - 1>::pow(n1);  // per face

  const Base base_;
  const T tau_;

 public:
  using constructor_value_type = typename Base::constructor_value_type;
  using data_type = typename Base::data_type;
  using Args = typename Base::Args;
  static constexpr unsigned int hyEdge_dim() { return dim; }
  static constexpr unsigned int n_glob_dofs_per_node() { return Base::n_glob_dofs_per_node(); }
  static constexpr unsigned int n_loc = 2 * dim * n_face;
  static constexpr unsigned int n_ldofs = (dim + 1) * n_vol;
  static_assert(n_ldofs == Base::Load::size(), "The local unknowns are those of Diffusion.");

  DiffusionTensor(const constructor_value_type& tau = 1.) : base_(tau), tau_(tau) {}

  template <unsigned int space_dim>
  static constexpr bool constrained_node(const unsigned int node_type)
  {
    return Base::template constrained_node<space_dim>(node_type);
  }

  /*!***********************************************************************************************
   * \brief   The matrix-free local operator of one box-shaped hyperedge.
   ************************************************************************************************/
  struct Operator
  {
    using Mat1 = std::array<T, n1 * n1>;  // row-major: (i_d, j_d)
    T tau, kE_inv;                        // tau and 1 / (k |E|)
    std::array<std::array<T, n1>, 2> b;   // endpoint values b_s(i), s = 0, 1
    std::array<T, dim> face_area;         // |F_d|, both faces in direction d
    std::array<std::array<T, 2>, dim> n;  // normal component n_d of face 2 d + s
    std::array<Mat1, dim> G, H, V;        // 1D derivative, its partner, eigenvectors of X_d
    std::array<std::array<T, n1>, dim> lambda;

    template <unsigned int k, class SolverT, class Edge, class Args>
    void apply(const SolverT&,
               Edge&,
               const SmallMat<n_loc, k, T>& in,
               const Args&,
               SmallMat<n_loc, k, T>& out) const
    {
      for (unsigned int col = 0; col < k; ++col)
      {
        SmallVec<n_loc, T> lambda_col, flux;
        for (unsigned int i = 0; i < n_loc; ++i)
          lambda_col[i] = in(i, col);
        flux_of(lambda_col, SmallVec<n_ldofs, T>(), flux);
        for (unsigned int i = 0; i < n_loc; ++i)
          out(i, col) = flux[i];
      }
    }
    template <class SolverT, class Edge, class Args>
    void residual(const SolverT& solver,
                  Edge& edge,
                  const SmallVec<n_loc, T>& lambda_in,
                  const Args& args,
                  SmallVec<n_loc, T>& r) const
    {
      flux_of(lambda_in, solver.local_load(edge, args), r);
    }
    template <class SolverT, class Edge, class Args>
    SmallVec<n_ldofs, T> local_solution(const SolverT& solver,
                                        Edge& edge,
                                        const SmallVec<n_loc, T>& lambda_in,
                                        const Args& args) const
    {
      return solve(lambda_in, solver.local_load(edge, args));
    }

   private:
    // r = C x + D lambda for the local solution x of the trace lambda and the load
    void flux_of(const SmallVec<n_loc, T>& lambda_in,
                 const SmallVec<n_ldofs, T>& load,
                 SmallVec<n_loc, T>& r) const
    {
      const SmallVec<n_ldofs, T> x = solve(lambda_in, load);
      // flux = C x + D lambda: restriction of -(n_d q_d + tau (u - lambda)) to every face, the
      // sign with which the condensed operator is positive semi-definite (as Diffusion)
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int s = 0; s < 2; ++s)
        {
          const unsigned int face = 2 * d + s;
          for (unsigned int j = 0; j < n_face; ++j)
            r[face * n_face + j] = tau * face_area[d] * lambda_in[face * n_face + j];
          for (unsigned int i = 0; i < n_vol; ++i)
          {
            const unsigned int i_d = digit(i, d);
            r[face * n_face + face_index(i, d)] -=
              face_area[d] * b[s][i_d] * (n[d][s] * x[d * n_vol + i] + tau * x[dim * n_vol + i]);
          }
        }
    }
    // the local solution x = A^-1 (B lambda + l)
    SmallVec<n_ldofs, T> solve(const SmallVec<n_loc, T>& lambda_in,
                               const SmallVec<n_ldofs, T>& load) const
    {
      // b = B lambda + l: extension of the face values along their direction
      SmallVec<n_ldofs, T> rhs = load;
      for (unsigned int d = 0; d < dim; ++d)
        for (unsigned int s = 0; s < 2; ++s)
        {
          const unsigned int face = 2 * d + s;
          for (unsigned int i = 0; i < n_vol; ++i)
          {
            const T value =
              face_area[d] * b[s][digit(i, d)] * lambda_in[face * n_face + face_index(i, d)];
            rhs[dim * n_vol + i] += tau * value;
            rhs[d * n_vol + i] -= n[d][s] * value;
          }
        }
      // u = S_u^-1 (b_u - sum_d H_d b_q_d / (k |E|)), by fast diagonalization
      std::array<T, n_vol> u, work;
      for (unsigned int i = 0; i < n_vol; ++i)
        u[i] = rhs[dim * n_vol + i];
      for (unsigned int d = 0; d < dim; ++d)
      {
        for (unsigned int i = 0; i < n_vol; ++i)
          work[i] = rhs[d * n_vol + i];
        along(H[d], false, d, work);
        for (unsigned int i = 0; i < n_vol; ++i)
          u[i] -= kE_inv * work[i];
      }
      for (unsigned int d = 0; d < dim; ++d)
        along(V[d], true, d, u);
      for (unsigned int i = 0; i < n_vol; ++i)
      {
        T sum = (T)0.;
        for (unsigned int d = 0; d < dim; ++d)
          sum += lambda[d][digit(i, d)];
        u[i] /= sum;
      }
      for (unsigned int d = 0; d < dim; ++d)
        along(V[d], false, d, u);
      // q_d = (b_q_d + G_d u) / (k |E|)
      SmallVec<n_ldofs, T> x;
      for (unsigned int d = 0; d < dim; ++d)
      {
        work = u;
        along(G[d], false, d, work);
        for (unsigned int i = 0; i < n_vol; ++i)
          x[d * n_vol + i] = kE_inv * (rhs[d * n_vol + i] + work[i]);
      }
      for (unsigned int i = 0; i < n_vol; ++i)
        x[dim * n_vol + i] = u[i];
      return x;
    }

   private:
    // n1^d, the step of direction d in the volume index
    static constexpr unsigned int step(const unsigned int d)
    {
      unsigned int result = 1;
      for (unsigned int c = 0; c < d; ++c)
        result *= n1;
      return result;
    }
    // digit d of the volume index i (direction 0 fastest, as Hypercube::index_decompose)
    static constexpr unsigned int digit(const unsigned int i, const unsigned int d)
    {
      return (i / step(d)) % n1;
    }
    // the face index of i on a face in direction d: its digits without digit d
    static constexpr unsigned int face_index(const unsigned int i, const unsigned int d)
    {
      const unsigned int stride = step(d);
      return i % stride + (i / (stride * n1)) * stride;
    }
    // v = M v (or M^T v) along direction d
    static void along(const Mat1& m,
                      const bool transposed_m,
                      const unsigned int d,
                      std::array<T, n_vol>& v)
    {
      const unsigned int stride = step(d);
      const std::array<T, n_vol> in = v;
      for (unsigned int i = 0; i < n_vol; ++i)
      {
        const unsigned int i_d = digit(i, d), base = i - i_d * stride;
        T sum = (T)0.;
        for (unsigned int m_d = 0; m_d < n1; ++m_d)
          sum += (transposed_m ? m[m_d * n1 + i_d] : m[i_d * n1 + m_d]) * in[base + m_d * stride];
        v[i] = sum;
      }
    }
  };

  /*!***********************************************************************************************
   * \brief   Set up the one-dimensional matrices and the eigendecompositions of one hyperedge.
   ************************************************************************************************/
  template <class Edge>
  Operator local_operator(Edge& edge, GlobalLoop::Stiffness, const Args& args) const
  {
    constexpr unsigned int space_dim = std::remove_cvref_t<decltype(Edge::geometry)>::space_dim();
    using parameters = parametersT<space_dim, T>;
    static_assert(
      requires { parameters::constant_coefficients; },
      "DiffusionTensor needs a diffusion coefficient constant on every hyperedge: the "
      "parameters must declare constant_coefficients.");
    static_assert(parameters::constant_coefficients,
                  "DiffusionTensor needs a diffusion coefficient constant on every hyperedge.");
    const auto& geometry = edge.geometry;
    const auto r = geometry.mat_r();
    const T area = geometry.area();
    const T k = parameters::inverse_diffusion_coeff(geometry.map_ref_to_phys(Point<dim, T>((T)0.5)),
                                                    args.time);

    Operator op;
    op.tau = tau_;
    op.kE_inv = (T)1. / (k * area);
    for (unsigned int s = 0; s < 2; ++s)
      for (unsigned int a = 0; a < n1; ++a)
        op.b[s][a] = integrator::shape_fcts_at_bdr[a][s];
    for (unsigned int d = 0; d < dim; ++d)
    {
      if constexpr (hy_debug_mode)
        for (unsigned int c = 0; c < dim; ++c)
        {
          hy_assert(c == d || std::abs(r(c, d)) <= 1e-12 * std::abs(r(d, d)),
                    "DiffusionTensor needs box-shaped hyperedges (mat_r diagonal).");
          hy_assert(c == d || std::abs(geometry.local_normal(2 * d)[c]) <= 1e-12,
                    "DiffusionTensor needs axis-aligned local normals.");
        }
      op.face_area[d] = geometry.face_area(2 * d);
      for (unsigned int s = 0; s < 2; ++s)
        op.n[d][s] = geometry.local_normal(2 * d + s)[d];
      // G_d(i, j) = |E| / r_d int phi_i' phi_j, and its partner H_d = -G_d + boundary terms
      typename Operator::Mat1 x;
      for (unsigned int i = 0; i < n1; ++i)
        for (unsigned int j = 0; j < n1; ++j)
        {
          const T g = area / r(d, d) * integrator::integrate_1D_Dphiphi(i, j);
          T boundary = (T)0., face_mass = (T)0.;
          for (unsigned int s = 0; s < 2; ++s)
          {
            boundary += op.face_area[d] * op.n[d][s] * op.b[s][i] * op.b[s][j];
            face_mass += op.face_area[d] * op.b[s][i] * op.b[s][j];
          }
          op.G[d][i * n1 + j] = g;
          op.H[d][i * n1 + j] = -g + boundary;
          x[i * n1 + j] = tau_ * face_mass;
        }
      // X_d = tau |F_d| sum_s b_s b_s^T + H_d G_d / (k |E|), symmetric up to quadrature rounding
      for (unsigned int i = 0; i < n1; ++i)
        for (unsigned int j = 0; j < n1; ++j)
          for (unsigned int m = 0; m < n1; ++m)
            x[i * n1 + j] += op.kE_inv * op.H[d][i * n1 + m] * op.G[d][m * n1 + j];
      for (unsigned int i = 0; i < n1; ++i)
        for (unsigned int j = i + 1; j < n1; ++j)
          x[i * n1 + j] = x[j * n1 + i] = (T)0.5 * (x[i * n1 + j] + x[j * n1 + i]);
      detail::symmetric_eigen<n1, T>(x, op.V[d], op.lambda[d]);
    }
    return op;
  }
  // The operator is a function of the local problem: shared as those of Diffusion.
  template <class Edge>
    requires requires(const Base& base, Edge& edge, const Args& args) {
      base.operator_slot(edge, GlobalLoop::Stiffness{}, args);
    }
  unsigned int n_operator_slots(Edge& edge, GlobalLoop::Stiffness kind, const Args& args) const
  {
    return base_.n_operator_slots(edge, kind, args);
  }
  template <class Edge>
    requires requires(const Base& base, Edge& edge, const Args& args) {
      base.operator_slot(edge, GlobalLoop::Stiffness{}, args);
    }
  unsigned int operator_slot(Edge& edge, GlobalLoop::Stiffness kind, const Args& args) const
  {
    return base_.operator_slot(edge, kind, args);
  }
  template <class Edge>
  SmallVec<n_ldofs, T> local_load(Edge& edge, const Args& args) const
  {
    return base_.local_load(edge, args);
  }
  template <class Edge, class OperatorT>
  auto error_sq(Edge& edge,
                const OperatorT& op,
                const SmallVec<n_loc, T>& lambda,
                const Args& args) const
  {
    return base_.error_sq(edge, op, lambda, args);
  }
};  // end of class DiffusionTensor

}  // end of namespace LocalSolver
