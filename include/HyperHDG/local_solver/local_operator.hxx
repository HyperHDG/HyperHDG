#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/wrapper/lapack.hxx>

#include <array>
#include <cstddef>

/*!*************************************************************************************************
 * \brief   Layouts of local operators: what a local solver stores of its condensed local problem.
 *
 * A local solver of the hybridized kind (Diffusion, ...) describes the local problem of a
 * hyperedge by its blocks (\c Blocks)
 *
 *   A x = B lambda + l,    flux = C x + D lambda,
 *
 * with the local unknowns x (n_ldofs), the trace lambda (n_loc), and the local load l, which
 * depends on the hyperedge, the time and the data (local_load of the solver). A layout turns the
 * blocks into the operator GlobalLoop::Generic applies: any type with
 *
 *   template <unsigned n_ldofs, unsigned n_loc, typename T> static auto make(Blocks<...>);
 *
 * returning an operator with apply<k>(edge, in, args, out), residual(edge, lambda, load, args, r)
 * and optionally local_solution(edge, lambda, load, args). It may also provide make_streaming, a
 * representation for an operator that is used once (a sweep without prepare); without it, make
 * serves both. The solver is templated on its layout and has a default; users can write their own
 * (e.g. single-precision or symmetric storage) without touching the solver.
 *
 * Provided: \c Condensed (the default of Diffusion) stores S = C A^-1 B + D and K = C A^-1, so
 * that apply is S lambda and residual is S lambda + K l, without factors; \c Factored stores the
 * LU factors of A with B, C, D, as the local problem is set up, and also provides local_solution.
 **************************************************************************************************/
namespace LocalOperator
{
/*!*************************************************************************************************
 * \brief   The blocks of a condensed local problem, A x = B lambda + l, flux = C x + D lambda.
 **************************************************************************************************/
template <unsigned int n_ldofs, unsigned int n_loc, typename T>
struct Blocks
{
  SmallSquareMat<n_ldofs, T> A;
  SmallMat<n_ldofs, n_loc, T> B;
  SmallMat<n_loc, n_ldofs, T> C;
  SmallSquareMat<n_loc, T> D;
};

namespace detail
{
/*!*************************************************************************************************
 * \brief   LU factors of A (getrf), in place.
 **************************************************************************************************/
template <unsigned int n, typename T, std::size_t n_ipiv>
void factorize(SmallSquareMat<n, T>& lu, std::array<int, n_ipiv>& ipiv)
{
  static_assert(n_ipiv == n, "One pivot per row.");
  Wrapper::lapack_factorize<n, T>(lu.data(), ipiv);
}
/*!*************************************************************************************************
 * \brief   rhs = A^-1 rhs (trans = 'N') or A^-T rhs (trans = 'T'), from the LU factors.
 **************************************************************************************************/
template <unsigned int n, unsigned int k, typename T, std::size_t n_ipiv>
void solve(const SmallSquareMat<n, T>& lu,
           const std::array<int, n_ipiv>& ipiv,
           SmallMat<n, k, T>& rhs,
           const char trans = 'N')
{
  // getrs reads the factors and pivots only; the wrapper takes non-const pointers.
  Wrapper::lapack_solve_factored(n, k, const_cast<T*>(lu.data().data()),
                                 const_cast<int*>(ipiv.data()), rhs.data().data(), trans);
}
}  // end of namespace detail

/*!*************************************************************************************************
 * \brief   The blocks as the local problem is set up: LU factors of A, B, C, D.
 *
 * Cheap to build (one factorization), expensive to apply (a solve with A per application), and
 * it provides local_solution. Used for streaming by Condensed.
 **************************************************************************************************/
struct Factored
{
  template <unsigned int n_ldofs, unsigned int n_loc, typename T>
  struct Operator
  {
    SmallSquareMat<n_ldofs, T> lu;
    std::array<int, n_ldofs> ipiv;
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
      detail::solve(lu, ipiv, x);
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
      detail::solve(lu, ipiv, x);
      return x;
    }
  };

  template <unsigned int n_ldofs, unsigned int n_loc, typename T>
  static Operator<n_ldofs, n_loc, T> make(const Blocks<n_ldofs, n_loc, T>& blocks)
  {
    Operator<n_ldofs, n_loc, T> op;
    op.lu = blocks.A;
    detail::factorize(op.lu, op.ipiv);
    op.B = blocks.B;
    op.C = blocks.C;
    op.D = blocks.D;
    return op;
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
  template <unsigned int n_ldofs, unsigned int n_loc, typename T>
  struct Operator
  {
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

  template <unsigned int n_ldofs, unsigned int n_loc, typename T>
  static Operator<n_ldofs, n_loc, T> make(const Blocks<n_ldofs, n_loc, T>& blocks)
  {
    SmallSquareMat<n_ldofs, T> lu = blocks.A;
    std::array<int, n_ldofs> ipiv;
    detail::factorize(lu, ipiv);
    Operator<n_ldofs, n_loc, T> op;
    SmallMat<n_ldofs, n_loc, T> ainv_b = blocks.B;
    detail::solve(lu, ipiv, ainv_b);
    op.S = blocks.C * ainv_b + blocks.D;
    // K = C A^-1, i.e. K^T = A^-T C^T: one solve with the transposed factors
    SmallMat<n_ldofs, n_loc, T> kt = transposed(blocks.C);
    detail::solve(lu, ipiv, kt, 'T');
    op.K = transposed(kt);
    return op;
  }
  template <unsigned int n_ldofs, unsigned int n_loc, typename T>
  static auto make_streaming(const Blocks<n_ldofs, n_loc, T>& blocks)
  {
    return Factored::make(blocks);
  }
};

}  // end of namespace LocalOperator
