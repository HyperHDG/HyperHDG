#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/global_loop/kinds.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/wrapper/lapack.hxx>

#include <array>
#include <cmath>

namespace LocalSolver
{
/*!*************************************************************************************************
 * \brief   LDG-H diffusion as a local-operator solver (prototype).
 *
 * Reuses the assembly of \c Diffusion (local matrix, loads, error integrals) and exposes the
 * local-operator protocol of \c GlobalLoop::LocalOperatorLoop:
 * - \c local_operator(edge, Stiffness, args) assembles and LU-factorizes the local matrix once and
 *   tabulates the coupling blocks B (trace -> local rhs), C (local solution -> flux) and the
 *   diagonal D (-tau * face area), i.e. the boundary integrals are evaluated once per operator.
 *   The condensed operator is S = C A^{-1} B + D. The operator is a trivially copyable value.
 * - \c Operator::apply<k> applies S to k traces at once (assembly uses k = n_loc with the
 *identity), \c Operator::residual adds the local load before the solve.
 * - \c local_load(edge, args) is the right-hand side of the local problem: volume source and
 *   Dirichlet data (as in Diffusion, Dirichlet values enter through the local load; the trace
 *   on Dirichlet nodes is constrained to zero by the loop).
 * No boundary masking happens here; constraints are the loop's business.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT = DiffusionParametersDefault,
          typename lSol_float_t = double>
class DiffusionOperator
: public Diffusion<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>
{
  using Base = Diffusion<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>;
  using Base::n_loc_dofs_;
  using Base::n_shape_bdr_;
  using Base::n_shape_fct_;
  using Base::tau_;
  using typename Base::integrator;

 public:
  using Base::Base;  // DiffusionOperator(tau)
  using Base::hyEdge_dim;
  using Base::n_glob_dofs_per_node;
  using typename Base::constructor_value_type;
  using typename Base::data_type;

  /*!***********************************************************************************************
   * \brief   Runtime arguments forwarded by the loop to every call.
   ************************************************************************************************/
  struct Args
  {
    lSol_float_t time = 0.;
  };

  static constexpr unsigned int n_loc = 2 * hyEdge_dimT * n_shape_bdr_;
  using Load = SmallVec<n_loc_dofs_, lSol_float_t>;

  /*!***********************************************************************************************
   * \brief   Default constraint predicate: Dirichlet node types of the parameters.
   ************************************************************************************************/
  template <unsigned int space_dim>
  static constexpr bool constrained_node(const unsigned int node_type)
  {
    return Base::template is_dirichlet<parametersT<space_dim, lSol_float_t> >(node_type);
  }

  /*!***********************************************************************************************
   * \brief   The immutable per-edge operator: LU factors of A plus the coupling blocks.
   ************************************************************************************************/
  struct Operator
  {
    SmallSquareMat<n_loc_dofs_, lSol_float_t> lu;
    std::array<int, n_loc_dofs_> ipiv;
    SmallMat<n_loc_dofs_, n_loc, lSol_float_t> B;
    SmallMat<n_loc, n_loc_dofs_, lSol_float_t> C;
    SmallVec<n_loc, lSol_float_t> D;

    template <unsigned int k, class Edge>
    void apply(const Edge&,
               const SmallMat<n_loc, k, lSol_float_t>& in,
               const Args&,
               SmallMat<n_loc, k, lSol_float_t>& out) const
    {
      SmallMat<n_loc_dofs_, k, lSol_float_t> coeffs = B * in;
      solve(coeffs);
      out = C * coeffs;
      for (unsigned int col = 0; col < k; ++col)
        for (unsigned int row = 0; row < n_loc; ++row)
          out(row, col) += D[row] * in(row, col);
    }
    template <class Edge>
    void residual(const Edge&,
                  const SmallVec<n_loc, lSol_float_t>& lambda,
                  const Load& load,
                  const Args&,
                  SmallVec<n_loc, lSol_float_t>& r) const
    {
      SmallVec<n_loc_dofs_, lSol_float_t> coeffs = B * lambda + load;
      solve(coeffs);
      r = C * coeffs;
      for (unsigned int row = 0; row < n_loc; ++row)
        r[row] += D[row] * lambda[row];
    }
    template <class Edge>
    SmallVec<n_loc_dofs_, lSol_float_t> local_solution(const Edge&,
                                                       const SmallVec<n_loc, lSol_float_t>& lambda,
                                                       const Load& load,
                                                       const Args&) const
    {
      SmallVec<n_loc_dofs_, lSol_float_t> coeffs = B * lambda + load;
      solve(coeffs);
      return coeffs;
    }

   private:
    template <unsigned int k>
    void solve(SmallMat<n_loc_dofs_, k, lSol_float_t>& rhs) const
    {
      // dgetrs reads the factors and pivots only; the wrapper takes non-const pointers.
      Wrapper::lapack_solve_factored(n_loc_dofs_, k, const_cast<lSol_float_t*>(lu.data().data()),
                                     const_cast<int*>(ipiv.data()), rhs.data().data());
    }
  };

  /*!***********************************************************************************************
   * \brief   Build the operator of one hyperedge: factorize A, tabulate B, C, D.
   ************************************************************************************************/
  template <class Edge>
  Operator local_operator(Edge& edge, GlobalLoop::Stiffness, const Args& args) const
  {
    Operator op;
    op.lu = Base::assemble_loc_matrix(tau_, edge, args.time);
    try
    {
      Wrapper::lapack_factorize<n_loc_dofs_, lSol_float_t>(op.lu.data(), op.ipiv);
    }
    catch (Wrapper::LAPACKexception& exc)
    {
      hy_assert(0 == 1, exc.what() << std::endl
                                   << "This can happen if quadrature is too inaccurate!");
      throw exc;
    }

    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      const auto normal = edge.geometry.local_normal(face);
      for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      {
        const unsigned int col = face * n_shape_bdr_ + j;
        op.D[col] = -tau_ * edge.geometry.face_area(face);
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
        {
          const lSol_float_t integral =
            integrator::template integrate_bdr_phipsi<decltype(Edge::geometry)>(i, j, face,
                                                                                edge.geometry);
          op.B(hyEdge_dimT * n_shape_fct_ + i, col) += tau_ * integral;
          op.C(col, hyEdge_dimT * n_shape_fct_ + i) += tau_ * integral;
          for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
          {
            op.B(dim * n_shape_fct_ + i, col) -= normal[dim] * integral;
            op.C(col, dim * n_shape_fct_ + i) += normal[dim] * integral;
          }
        }
      }
    }
    return op;
  }
  /*!***********************************************************************************************
   * \brief   Share one operator among congruent hyperedges (opt-in by the parameters).
   *
   * Only available if the parameters declare \c constant_coefficients: then the local operator
   * depends on the hyperedge's shape only (constraints and loads are not part of it), and the key
   * is the span vectors rounded to 12 digits. Elements of a uniform mesh thus share one operator,
   * which is what factorize_once of Elliptic::trace_to_flux_mat assumes, here by construction.
   ************************************************************************************************/
  template <class Edge>
    requires(
      parametersT<decltype(Edge::geometry)::space_dim(), lSol_float_t>::constant_coefficients)
  std::array<long long, hyEdge_dimT* decltype(Edge::geometry)::space_dim()>
  operator_key(Edge& edge, GlobalLoop::Stiffness, const Args&) const
  {
    constexpr unsigned int space_dim = decltype(Edge::geometry)::space_dim();
    std::array<long long, hyEdge_dimT * space_dim> key;
    for (unsigned int i = 0; i < hyEdge_dimT; ++i)
    {
      const auto span = edge.geometry.span_vec(i);
      for (unsigned int d = 0; d < space_dim; ++d)
        key[i * space_dim + d] = (long long)std::llround(span[d] * 1e12);
    }
    return key;
  }
  /*!***********************************************************************************************
   * \brief   Local right-hand side: volume source and Dirichlet data.
   ************************************************************************************************/
  template <class Edge>
  Load local_load(Edge& edge, const Args& args) const
  {
    return Base::assemble_rhs_from_global_rhs(edge, args.time);
  }
  /*!***********************************************************************************************
   * \brief   Squared L2 error of the primal unknown on one hyperedge (SmallVec: has operator+).
   ************************************************************************************************/
  template <class Edge>
  SmallVec<1, lSol_float_t> error_sq(Edge& edge,
                                     const Operator& op,
                                     const SmallVec<n_loc, lSol_float_t>& lambda,
                                     const Args& args) const
  {
    using parameters = parametersT<decltype(Edge::geometry)::space_dim(), lSol_float_t>;
    const SmallVec<n_loc_dofs_, lSol_float_t> coefficients =
      op.local_solution(edge, lambda, local_load(edge, args), args);
    std::array<lSol_float_t, n_shape_fct_> coeffs;
    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[i + hyEdge_dimT * n_shape_fct_];
    SmallVec<1, lSol_float_t> error;
    error[0] = integrator::template integrate_vol_diffsquare_discana<
      Point<decltype(Edge::geometry)::space_dim(), lSol_float_t>, decltype(Edge::geometry),
      parameters::analytic_result, Point<hyEdge_dimT, lSol_float_t> >(coeffs, edge.geometry,
                                                                      args.time);
    return error;
  }
};  // end of class DiffusionOperator

}  // end of namespace LocalSolver
