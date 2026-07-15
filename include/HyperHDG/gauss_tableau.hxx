#pragma once  // Ensure that file is included only once in a single compilation.

/*!*************************************************************************************************
 * \brief   Gauss-Legendre collocation tableau in the eigenbasis of the Butcher matrix.
 *
 * Port of the tableau setup of ~/phd/hoRK/horkirk.c (HORKIRKBuildTableau): Gauss nodes/weights by
 * Newton iteration on the Legendre polynomial, Butcher matrix A = W V^{-1} from the collocation
 * conditions, eigendecomposition A = T diag(theta) T^{-1} via LAPACK zgeev (trick 1), conjugate
 * eigenvalue pairing with enforced v_partner = conj(v_rep) and phase-fixed real eigenvector
 * columns (trick 2) so only the ceil(s/2) representatives are ever solved.
 *
 * Difference to horkirk: horkirk works with stage DERIVATIVES (update weights b^T T); the HDG
 * Gauss stepping (hdg_gauss.pdf) works with stage VALUES, so the recombination weights here are
 * w = (b^T T) / theta = (b^T A^{-1} T) and the endpoint update reads
 *   y+ = affine * y^n + sum_reps mult_l * Re(w_l * y_l),   affine = 1 - b^T A^{-1} 1 = (-1)^s.
 **************************************************************************************************/

#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/wrapper/lapack.hxx>

#include <cmath>
#include <complex>
#include <vector>

namespace Gauss
{

/*!*************************************************************************************************
 * \brief   Time argument of a Gauss stage solve: physical end-of-step time plus stage index.
 *
 * Passed through the generic global-loop entries (trace_to_flux_mat / residual_flux2 / set_data)
 * in place of the plain time scalar; stage-aware local solvers unpack it, everything else keeps
 * receiving plain times.
 **************************************************************************************************/
struct StageTime
{
  double time;
  unsigned int stage;
};


/*!*************************************************************************************************
 * \brief   Gauss-Legendre nodes c (mapped to [0,1]) and weights b (sum to 1).
 **************************************************************************************************/
inline void gauss_legendre(const unsigned int s, std::vector<double>& c, std::vector<double>& b)
{
  c.resize(s);
  b.resize(s);
  for (unsigned int i = 0; i < s; i++)
  {
    double x = std::cos(M_PI * (i + 0.75) / (s + 0.5));  // initial guess
    double p0 = 1., p1 = x, dp = 1.;
    for (unsigned int it = 0; it < 100; it++)
    {
      p0 = 1.;
      p1 = x;
      for (unsigned int k = 2; k <= s; k++)  // P_k recurrence; p1 = P_s, p0 = P_{s-1}
      {
        const double p2 = ((2. * k - 1.) * x * p1 - (k - 1.) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      dp = s * (x * p1 - p0) / (x * x - 1.);  // P_s'(x)
      const double dx = -p1 / dp;
      x += dx;
      if (std::abs(dx) < 1e-15)
        break;
    }
    const double weight = 2. / ((1. - x * x) * dp * dp);  // weight on [-1,1]
    c[i] = 0.5 * (x + 1.);                                // map to [0,1]
    b[i] = 0.5 * weight;
  }
}

/*!*************************************************************************************************
 * \brief   Eigenbasis tableau data of the s-stage Gauss collocation method (order 2s).
 **************************************************************************************************/
struct Tableau
{
  unsigned int s;                            // number of stages
  unsigned int n_reps;                       // number of solved representatives, ceil(s/2)
  std::vector<double> c;                     // collocation nodes in (0,1), length s
  std::vector<std::complex<double>> theta;   // representative Butcher eigenvalues, length n_reps
  std::vector<std::complex<double>> omega;   // load weights (T^{-1} 1)_rep
  std::vector<std::complex<double>> w;       // value-form recombination weights (b^T T)_rep/theta
  std::vector<std::complex<double>> tinv;    // T^{-1} rows of the reps, row-major n_reps x s
  std::vector<double> mult;                  // multiplicity: 2 for a conjugate pair, else 1
  double affine;                             // 1 - b^T A^{-1} 1 = R(infinity) = (-1)^s
};

/*!*************************************************************************************************
 * \brief   Build the tableau; see the file header. All dense work is on tiny s x s systems.
 **************************************************************************************************/
inline Tableau build_tableau(const unsigned int s)
{
  using cplx = std::complex<double>;
  constexpr double tol = 1e-9;
  Tableau tab;
  tab.s = s;

  std::vector<double> c, b;
  gauss_legendre(s, c, b);
  tab.c = c;

  // Butcher A = W V^{-1}, V_{jk} = c_j^k, W_{jk} = c_j^{k+1}/(k+1)  (column-major, complexified)
  std::vector<cplx> V(s * s), W(s * s), Vinv(s * s, 0.), A(s * s, 0.);
  for (unsigned int k = 0; k < s; k++)
    for (unsigned int j = 0; j < s; j++)
    {
      V[j + s * k] = std::pow(c[j], (double)k);
      W[j + s * k] = std::pow(c[j], (double)(k + 1)) / (k + 1.);
    }
  std::vector<int> ipiv(s);
  for (unsigned int i = 0; i < s; i++)
    Vinv[i + s * i] = 1.;
  Wrapper::lapack_factorize((int)s, V.data(), ipiv.data());
  Wrapper::lapack_solve_factored((int)s, (int)s, V.data(), ipiv.data(), Vinv.data());
  for (unsigned int j = 0; j < s; j++)
    for (unsigned int i = 0; i < s; i++)
      for (unsigned int k = 0; k < s; k++)
        A[i + s * j] += W[i + s * k] * Vinv[k + s * j];

  // eigendecomposition A = T diag(theta) T^{-1}
  std::vector<cplx> T(s * s), eigval(s);
  Wrapper::lapack_eig((int)s, A.data(), eigval.data(), T.data());

  // pair conjugate eigenvalues; enforce v_partner = conj(v_rep) so the conjugation of the
  // stage solutions is exact and the partners never need solving
  std::vector<bool> taken(s, false), is_rep(s, false);
  std::vector<int> partner(s, -1);
  std::vector<unsigned int> reps;
  for (unsigned int i = 0; i < s; i++)
  {
    if (taken[i])
      continue;
    reps.push_back(i);
    is_rep[i] = true;
    taken[i] = true;
    for (unsigned int j = i + 1; j < s; j++)
      if (!taken[j] && std::abs(eigval[j] - std::conj(eigval[i])) < tol)
      {
        for (unsigned int k = 0; k < s; k++)
          T[k + s * j] = std::conj(T[k + s * i]);
        partner[i] = (int)j;
        taken[j] = true;
        break;
      }
  }

  // real (unpaired) eigenvalues: zgeev returns their eigenvector with an arbitrary complex
  // phase; rotate it real so conj(T) = T P holds exactly (only occurs for odd s)
  for (const unsigned int i : reps)
  {
    if (partner[i] >= 0)
      continue;
    double maxa = -1.;
    unsigned int km = 0;
    for (unsigned int k = 0; k < s; k++)
      if (std::abs(T[k + s * i]) > maxa)
      {
        maxa = std::abs(T[k + s * i]);
        km = k;
      }
    const cplx ph = T[km + s * i] / maxa;  // unit phase of largest entry
    for (unsigned int k = 0; k < s; k++)
      T[k + s * i] /= ph;
  }

  // omega = T^{-1} 1 and the full T^{-1} (both via the factored T; gesv-style copy)
  std::vector<cplx> Tcopy(T), omega_full(s, 1.), Tinv(s * s, 0.);
  for (unsigned int i = 0; i < s; i++)
    Tinv[i + s * i] = 1.;
  Wrapper::lapack_factorize((int)s, Tcopy.data(), ipiv.data());
  Wrapper::lapack_solve_factored((int)s, 1, Tcopy.data(), ipiv.data(), omega_full.data());
  Wrapper::lapack_solve_factored((int)s, (int)s, Tcopy.data(), ipiv.data(), Tinv.data());

  // value-form recombination weights w_j = (b^T T)_j / theta_j and the affine coefficient
  // 1 - b^T A^{-1} 1 = 1 - sum_j w_j omega_j (real by conjugate symmetry; must be (-1)^s)
  std::vector<cplx> w_full(s);
  cplx affine = 1.;
  for (unsigned int j = 0; j < s; j++)
  {
    cplx acc = 0.;
    for (unsigned int i = 0; i < s; i++)
      acc += b[i] * T[i + s * j];
    w_full[j] = acc / eigval[j];
    affine -= w_full[j] * omega_full[j];
  }
  hy_assert(std::abs(affine.imag()) < 1e-10,
            "affine coefficient must be real; got imag = " << affine.imag());
  hy_assert(std::abs(affine.real() - (s % 2 ? -1. : 1.)) < 1e-10,
            "affine coefficient must be R(inf) = (-1)^s; got " << affine.real());

  tab.n_reps = (unsigned int)reps.size();
  hy_assert(tab.n_reps == (s + 1) / 2, "expected ceil(s/2) representatives");
  tab.affine = affine.real();
  for (const unsigned int i : reps)
  {
    tab.theta.push_back(eigval[i]);
    tab.omega.push_back(omega_full[i]);
    tab.w.push_back(w_full[i]);
    tab.mult.push_back(partner[i] >= 0 ? 2. : 1.);
    for (unsigned int j = 0; j < s; j++)
      tab.tinv.push_back(Tinv[i + s * j]);  // row i of T^{-1}
  }
  return tab;
}

}  // end of namespace Gauss
