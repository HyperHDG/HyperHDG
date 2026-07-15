#include <HyperHDG/gauss_tableau.hxx>
#include <HyperHDG/hy_assert.hxx>

#include <cmath>
#include <complex>

// Invariants of the Gauss collocation tableau in the eigenbasis of the Butcher matrix
// (Gauss::build_tableau, ported from ~/phd/hoRK/horkirk.c), for the stage-VALUE form used by
// the HDG Gauss stepping:  y+ = affine * y^n + sum_l mult_l * Re(w_l * y_l).
int main()
{
  using cplx = std::complex<double>;
  constexpr double tol = 1e-12;

  // s = 1 must reproduce the implicit-midpoint (== Crank-Nicolson trajectory) constants exactly.
  {
    const Gauss::Tableau t = Gauss::build_tableau(1);
    hy_check(t.n_reps == 1 && t.mult[0] == 1., "s=1: one real representative");
    hy_check(std::abs(t.c[0] - 0.5) < tol, "s=1: c = 1/2, got " << t.c[0]);
    hy_check(std::abs(t.theta[0] - cplx(0.5)) < tol, "s=1: theta = 1/2");
    hy_check(std::abs(t.omega[0] - cplx(1.)) < tol, "s=1: omega = 1");
    hy_check(std::abs(t.w[0] - cplx(2.)) < tol, "s=1: w = 2");
    hy_check(std::abs(t.affine + 1.) < tol, "s=1: affine = -1");
  }

  // s = 2: the eigenvalue pair is 1/4 +- i*sqrt(3)/12 (normalization-independent).
  {
    const Gauss::Tableau t = Gauss::build_tableau(2);
    hy_check(t.n_reps == 1 && t.mult[0] == 2., "s=2: one conjugate-pair representative");
    hy_check(std::abs(t.theta[0].real() - 0.25) < tol &&
               std::abs(std::abs(t.theta[0].imag()) - std::sqrt(3.) / 12.) < tol,
             "s=2: theta = 1/4 +- i*sqrt(3)/12, got " << t.theta[0].real() << " + "
                                                      << t.theta[0].imag() << "i");
  }

  for (unsigned int s = 1; s <= 4; ++s)
  {
    const Gauss::Tableau t = Gauss::build_tableau(s);
    hy_check(t.n_reps == (s + 1) / 2, "s=" << s << ": expected ceil(s/2) representatives");

    // affine = 1 - b^T A^{-1} 1 is the stability function at infinity, R(inf) = (-1)^s.
    hy_check(std::abs(t.affine - (s % 2 ? -1. : 1.)) < 1e-10,
             "s=" << s << ": affine must be (-1)^s, got " << t.affine);

    double mult_sum = 0.;
    for (unsigned int l = 0; l < t.n_reps; ++l)
    {
      mult_sum += t.mult[l];
      // omega = T^{-1} 1, i.e. the row sum of the T^{-1} rows handed out for the loads.
      cplx row = 0.;
      for (unsigned int j = 0; j < s; ++j)
        row += t.tinv[l * s + j];
      hy_check(std::abs(row - t.omega[l]) < 1e-10,
               "s=" << s << " rep " << l << ": sum of tinv row must equal omega");
    }
    hy_check(mult_sum == (double)s, "s=" << s << ": multiplicities must sum to s");

    // Exactness on y' = q t^{q-1}, y(0) = 0 for q = 1..s: stage values Y_j = c_j^q lie in the
    // collocation space, so the endpoint update must return y(1) = 1 exactly:
    //   affine * 0 + sum_l mult_l * Re(w_l * (T^{-1} Y)_l) = 1.
    for (unsigned int q = 1; q <= s; ++q)
    {
      double update = 0.;
      for (unsigned int l = 0; l < t.n_reps; ++l)
      {
        cplx ty = 0.;
        for (unsigned int j = 0; j < s; ++j)
          ty += t.tinv[l * s + j] * std::pow(t.c[j], (double)q);
        update += t.mult[l] * (t.w[l] * ty).real();
      }
      hy_check(std::abs(update - 1.) < 1e-10,
               "s=" << s << ", q=" << q << ": endpoint update must be exact, got " << update);
    }

    // Exactness on y' = 0 (constant state): affine + sum_l mult_l * Re(w_l * omega_l) = 1.
    double con = t.affine;
    for (unsigned int l = 0; l < t.n_reps; ++l)
      con += t.mult[l] * (t.w[l] * t.omega[l]).real();
    hy_check(std::abs(con - 1.) < 1e-10,
             "s=" << s << ": constant state must be reproduced, got " << con);
  }

  return 0;
}
