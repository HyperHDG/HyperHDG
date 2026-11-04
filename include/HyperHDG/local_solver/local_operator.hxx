#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/wrapper/lapack.hxx>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>

/*!*************************************************************************************************
 * \brief   Layouts of local operators: what a local solver stores of its condensed local problem.
 *
 * A local solver of the hybridized kind (Diffusion, ...) describes the local problem of a
 * hyperedge by its blocks (\c Blocks)
 *
 *   A x = B lambda + l,    flux = C x + D lambda,
 *
 * with the local unknowns x (n_ldofs), the trace lambda (n_loc), and the local load l, which
 * depends on the hyperedge, the time and the data (local_load of the solver). If the leading
 * n_eliminate unknowns form a diagonal block of A (e.g. the flux unknowns of an LDG-H method with
 * an orthonormal basis, affine elements and coefficients constant on the element), the solver says
 * so in the type of its blocks; the factorization then eliminates them explicitly and factorizes
 * the Schur complement of the others only (\c detail::Factors), in every layout.
 *
 * A layout turns the blocks into the operator GlobalLoop::Generic applies: any type with
 *
 *   template <class BlocksT> static auto make(const BlocksT&);
 *
 * returning an operator with apply<k>(edge, in, args, out), residual(edge, lambda, load, args, r)
 * and optionally local_solution(edge, lambda, load, args). It may also provide make_streaming, a
 * representation for an operator that is used once (a sweep without prepare); without it, make
 * serves both. The solver is templated on its layout and has a default; users can write their own
 * (e.g. single-precision or symmetric storage) without touching the solver.
 *
 * Provided: \c Condensed (the default of Diffusion) stores S = C A^-1 B + D and K = C A^-1, so
 * that apply is S lambda and residual is S lambda + K l, without factors; \c Factored stores the
 * factors of A with B, C, D, as the local problem is set up, and also provides local_solution.
 **************************************************************************************************/
namespace LocalOperator
{
/*!*************************************************************************************************
 * \brief   The blocks of a condensed local problem, A x = B lambda + l, flux = C x + D lambda.
 *
 * \tparam  n_eliminateT  The leading n_eliminateT unknowns form a diagonal block of A, uncoupled
 *                        among themselves, and are eliminated explicitly; 0: factorize A whole.
 **************************************************************************************************/
template <unsigned int n_ldofsT, unsigned int n_locT, typename T, unsigned int n_eliminateT = 0>
struct Blocks
{
  static constexpr unsigned int n_ldofs = n_ldofsT, n_loc = n_locT, n_eliminate = n_eliminateT;
  using value_type = T;
  static_assert(n_eliminate < n_ldofs, "At least one unknown is left to the factorization.");

  SmallSquareMat<n_ldofs, T> A;
  SmallMat<n_ldofs, n_loc, T> B;
  SmallMat<n_loc, n_ldofs, T> C;
  SmallSquareMat<n_loc, T> D;
};

namespace detail
{
/*!*************************************************************************************************
 * \brief   rhs = M^-1 rhs (trans = 'N') or M^-T rhs (trans = 'T'), from LU factors of M (getrs;
 *          LAPACK or the generic LU of dense_la.hxx, by the scalar type).
 **************************************************************************************************/
template <unsigned int n, unsigned int k, typename T, std::size_t n_ipiv>
void lu_solve(const SmallSquareMat<n, T>& lu,
              const std::array<int, n_ipiv>& ipiv,
              SmallMat<n, k, T>& rhs,
              const char trans = 'N')
{
  static_assert(n_ipiv == n, "One pivot per row.");
  LU::solve_factored<n, k, T>(lu.data(), ipiv, rhs.data(), trans);
}

/*!*************************************************************************************************
 * \brief   Factors of the local matrix A, with the diagonal block of the leading n_eliminate
 *          unknowns eliminated explicitly.
 *
 * With A = [D_A, A_DR; A_RD, A_RR] (D_A diagonal, of size n_eliminate), A x = b is
 *   x_R = S^-1 (b_R - A_RD D_A^-1 b_D),  x_D = D_A^-1 (b_D - A_DR x_R),
 * with the Schur complement S = A_RR - A_RD D_A^-1 A_DR, of which LU factors are kept. Building
 * costs about 2 n_D n_R^2 + 2/3 n_R^3 instead of 2/3 n^3 flops, a solve about 2 n_R^2 + 4 n_D n_R
 * instead of 2 n^2 per right-hand side, and the factors (2 n_D + n_R) n_R + n_D numbers instead of
 * n^2. Debug builds check that D_A is diagonal (up to 100 epsilon of T, relative).
 **************************************************************************************************/
template <unsigned int n, unsigned int n_eliminate, typename T>
struct Factors
{
  static constexpr unsigned int n_d = n_eliminate, n_r = n - n_eliminate;

  SmallVec<n_d, T> d_inv;
  SmallMat<n_d, n_r, T> a_dr;
  SmallMat<n_r, n_d, T> a_rd;
  SmallSquareMat<n_r, T> lu;
  std::array<int, n_r> ipiv;

  Factors() = default;
  explicit Factors(const SmallSquareMat<n, T>& a)
  {
    for (unsigned int i = 0; i < n_d; ++i)
    {
      if constexpr (hy_debug_mode)
        for (unsigned int j = 0; j < n_d; ++j)
          hy_assert(
            i == j ||
              std::abs(a(i, j)) <= 100 * std::numeric_limits<T>::epsilon() * std::abs(a(i, i)),
            "Eliminated unknowns must form a diagonal block of the local matrix, but entry ("
              << i << ", " << j << ") is " << (double)std::abs(a(i, j)) << " in magnitude.");
      d_inv[i] = (T)1. / a(i, i);
      for (unsigned int j = 0; j < n_r; ++j)
      {
        a_dr(i, j) = a(i, n_d + j);
        a_rd(j, i) = a(n_d + j, i);
      }
    }
    for (unsigned int j = 0; j < n_r; ++j)
      for (unsigned int i = 0; i < n_r; ++i)
      {
        T schur = a(n_d + i, n_d + j);
        for (unsigned int k = 0; k < n_d; ++k)
          schur -= a_rd(i, k) * d_inv[k] * a_dr(k, j);
        lu(i, j) = schur;
      }
    LU::factorize<n_r, T>(lu.data(), ipiv);
  }

  /*!
   * rhs = A^-1 rhs (trans = 'N') or A^-T rhs (trans = 'T'); with A^T the roles of A_DR and A_RD
   * swap and S is solved transposed.
   */
  template <unsigned int k>
  void solve(SmallMat<n, k, T>& rhs, const char trans = 'N') const
  {
    if (trans == 'N')
      solve_impl<false>(rhs);
    else
      solve_impl<true>(rhs);
  }

 private:
  template <bool transposedT, unsigned int k>
  void solve_impl(SmallMat<n, k, T>& rhs) const
  {
    // the blocks of A (or of A^T) below and above the eliminated diagonal block
    auto lower = [this](const unsigned int i, const unsigned int j)
    {
      if constexpr (transposedT)
        return a_dr(j, i);
      else
        return a_rd(i, j);
    };
    auto upper = [this](const unsigned int i, const unsigned int j)
    {
      if constexpr (transposedT)
        return a_rd(j, i);
      else
        return a_dr(i, j);
    };
    SmallMat<n_r, k, T> x_r;
    for (unsigned int col = 0; col < k; ++col)
      for (unsigned int i = 0; i < n_r; ++i)
      {
        T value = rhs(n_d + i, col);
        for (unsigned int j = 0; j < n_d; ++j)
          value -= lower(i, j) * d_inv[j] * rhs(j, col);
        x_r(i, col) = value;
      }
    lu_solve(lu, ipiv, x_r, transposedT ? 'T' : 'N');
    for (unsigned int col = 0; col < k; ++col)
    {
      for (unsigned int i = 0; i < n_d; ++i)
      {
        T value = rhs(i, col);
        for (unsigned int j = 0; j < n_r; ++j)
          value -= upper(i, j) * x_r(j, col);
        rhs(i, col) = d_inv[i] * value;
      }
      for (unsigned int i = 0; i < n_r; ++i)
        rhs(n_d + i, col) = x_r(i, col);
    }
  }
};

/*!*************************************************************************************************
 * \brief   Factors of the whole local matrix A (no unknowns eliminated): its LU factors.
 **************************************************************************************************/
template <unsigned int n, typename T>
struct Factors<n, 0, T>
{
  SmallSquareMat<n, T> lu;
  std::array<int, n> ipiv;

  Factors() = default;
  explicit Factors(const SmallSquareMat<n, T>& a) : lu(a) { LU::factorize<n, T>(lu.data(), ipiv); }
  template <unsigned int k>
  void solve(SmallMat<n, k, T>& rhs, const char trans = 'N') const
  {
    lu_solve(lu, ipiv, rhs, trans);
  }
};

template <class BlocksT>
using FactorsOf = Factors<BlocksT::n_ldofs, BlocksT::n_eliminate, typename BlocksT::value_type>;
}  // end of namespace detail

/*!*************************************************************************************************
 * \brief   The blocks as the local problem is set up: the factors of A, B, C, D.
 *
 * Cheap to build (one factorization), expensive to apply (a solve with A per application), and
 * it provides local_solution. Used for streaming by Condensed.
 **************************************************************************************************/
struct Factored
{
  template <class BlocksT>
  struct Operator
  {
    static constexpr unsigned int n_ldofs = BlocksT::n_ldofs, n_loc = BlocksT::n_loc;
    using T = typename BlocksT::value_type;

    detail::FactorsOf<BlocksT> factors;
    SmallMat<n_ldofs, n_loc, T> B;
    SmallMat<n_loc, n_ldofs, T> C;
    SmallSquareMat<n_loc, T> D;

    template <unsigned int k, class Edge, class Args>
    void apply(const Edge&,
               const SmallMat<n_loc, k, T>& in,
               const Args&,
               SmallMat<n_loc, k, T>& out) const
    {
      SmallMat<n_ldofs, k, T> x = B * in;
      factors.solve(x);
      out = C * x + D * in;
    }
    template <class Edge, class Args>
    void residual(const Edge& edge,
                  const SmallVec<n_loc, T>& lambda,
                  const SmallVec<n_ldofs, T>& load,
                  const Args& args,
                  SmallVec<n_loc, T>& r) const
    {
      r = C * local_solution(edge, lambda, load, args) + D * lambda;
    }
    template <class Edge, class Args>
    SmallVec<n_ldofs, T> local_solution(const Edge&,
                                        const SmallVec<n_loc, T>& lambda,
                                        const SmallVec<n_ldofs, T>& load,
                                        const Args&) const
    {
      SmallVec<n_ldofs, T> x = B * lambda + load;
      factors.solve(x);
      return x;
    }
  };

  template <class BlocksT>
  static Operator<BlocksT> make(const BlocksT& blocks)
  {
    return Operator<BlocksT>{detail::FactorsOf<BlocksT>(blocks.A), blocks.B, blocks.C, blocks.D};
  }
};

/*!*************************************************************************************************
 * \brief   The condensed operator: S = C A^-1 B + D and K = C A^-1, no factors.
 *
 * Building costs a factorization, a solve with n_loc right-hand sides, and one with A^T for K;
 * applying costs one product with S (n_loc^2 instead of about 2 n_ldofs n_loc + n_ldofs^2 for
 * Factored), the residual one with K in addition. Smaller than Factored, since n_loc < n_ldofs.
 * No local_solution: reductions that need it build the Factored operator on the fly. A sweep
 * without prepare uses Factored (make_streaming), which is cheaper to build for one use.
 **************************************************************************************************/
struct Condensed
{
  template <class BlocksT>
  struct Operator
  {
    static constexpr unsigned int n_ldofs = BlocksT::n_ldofs, n_loc = BlocksT::n_loc;
    using T = typename BlocksT::value_type;

    SmallSquareMat<n_loc, T> S;
    SmallMat<n_loc, n_ldofs, T> K;

    template <unsigned int k, class Edge, class Args>
    void apply(const Edge&,
               const SmallMat<n_loc, k, T>& in,
               const Args&,
               SmallMat<n_loc, k, T>& out) const
    {
      out = S * in;
    }
    template <class Edge, class Args>
    void residual(const Edge&,
                  const SmallVec<n_loc, T>& lambda,
                  const SmallVec<n_ldofs, T>& load,
                  const Args&,
                  SmallVec<n_loc, T>& r) const
    {
      r = S * lambda + K * load;
    }
  };

  template <class BlocksT>
  static Operator<BlocksT> make(const BlocksT& blocks)
  {
    const detail::FactorsOf<BlocksT> factors(blocks.A);
    Operator<BlocksT> op;
    auto ainv_b = blocks.B;
    factors.solve(ainv_b);
    op.S = blocks.C * ainv_b + blocks.D;
    // K = C A^-1, i.e. K^T = A^-T C^T: one solve with the transposed factors
    auto kt = transposed(blocks.C);
    factors.solve(kt, 'T');
    op.K = transposed(kt);
    return op;
  }
  template <class BlocksT>
  static auto make_streaming(const BlocksT& blocks)
  {
    return Factored::make(blocks);
  }
};

}  // end of namespace LocalOperator
