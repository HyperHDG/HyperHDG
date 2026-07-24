#pragma once

// TestTimoWave4: the wave4 manufactured solution on the axis-aligned cross (cross2.geo) /
// single arm (single1.geo). Extracted from parameters.hxx so that PETSc-free consumers (the
// nanobind python module) can use it: everything except the PetscOptions-based Init() is
// plain header code. Include after <petsc.h> (as parameters.hxx does) to get Init().

#include <HyperHDG/dense_la.hxx>

#include <array>
#include <cmath>

template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave4
{
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};

  static constexpr param_float_t omega = 2*M_PI;

  // u_k = au[k] cos(w sigma) cos(w t + pu[k]), r_k = ar[k] cos(w sigma) cos(w t + pr[k])
  // on the axis-aligned cross (cross2.geo); derived in experiments/timowave4.m.
  // arms must come in +-pairs so the nodal forces n(0) = -i x r(0) balance.
  // distinct per-component phases: all of u, v, r, s, n, m are nonzero at t = 0
  static constexpr std::array<param_float_t, 3> au{1, 2, 3},  pu{0.3, 0.8, 1.3};
  static constexpr std::array<param_float_t, 3> ar{5, 7, 11}, pr{0.5, 1.1, 1.7};

  // Runtime SPATIAL phase of the standing-wave factor cos(w*(x+y+z) + px). The time phases
  // pu/pr cannot silence the boundary (cos(w t + p) never vanishes identically); px moves the
  // spatial nodes: px = -pi/2 turns the factor into sin(w*s), which is zero at every tip of
  // cross2/single1 (s in {0, +-1}) -- Dirichlet data identically zero for all t. The per-arm
  // forcing formulas below generalize by literally phase-shifting their trig factors (the
  // u-parts and m' still cancel at unit wave speed; the +-arm-pair junction balance is
  // pairwise smoothness through the center, independent of px). Default 0 = classic wave4.
  static inline param_float_t px = 0.;

#ifdef PETSC_VERSION_MAJOR
  static PetscErrorCode Init(const char*)
  {
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-wave4_px", &px, NULL));
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "wave4_px: %g\n", (double)px));
    PetscFunctionReturn(PETSC_SUCCESS);
  }
#endif

  static SmallVec<space_dimT, param_float_t> uvec(const param_float_t time) {
    SmallVec<space_dimT, param_float_t> res(0.);
    for (unsigned int k = 0; k < 3; k++)
      res[k] = au[k]*cos(omega*time + pu[k]);
    return res;
  }
  static SmallVec<space_dimT, param_float_t> dt_uvec(const param_float_t time) {
    SmallVec<space_dimT, param_float_t> res(0.);
    for (unsigned int k = 0; k < 3; k++)
      res[k] = -omega*au[k]*sin(omega*time + pu[k]);
    return res;
  }
  static SmallVec<space_dimT, param_float_t> rvec(const param_float_t time) {
    SmallVec<space_dimT, param_float_t> res(0.);
    for (unsigned int k = 0; k < 3; k++)
      res[k] = ar[k]*cos(omega*time + pr[k]);
    return res;
  }
  static SmallVec<space_dimT, param_float_t> dt_rvec(const param_float_t time) {
    SmallVec<space_dimT, param_float_t> res(0.);
    for (unsigned int k = 0; k < 3; k++)
      res[k] = -omega*ar[k]*sin(omega*time + pr[k]);
    return res;
  }

  // axis of the arm containing point (exactly one coordinate is nonzero off the center)
  static unsigned int arm_axis(const Point<space_dimT, param_float_t>& point) {
    if (point[0] != 0) return 0;
    if (point[1] != 0) return 1;
    return 2;
  }
  // e_ax x vec
  static SmallVec<space_dimT, param_float_t> cross_e(const unsigned int ax,
                                                     const SmallVec<space_dimT, param_float_t>& vec) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[(ax+1)%3] = -vec[(ax+2)%3];
    res[(ax+2)%3] =  vec[(ax+1)%3];
    return res;
  }

  // f = w sin(w x_ax) (e_ax x rvec), per arm (derived in experiments/timowave4.m)
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    const unsigned int ax = arm_axis(point);
    auto res = cross_e(ax, rvec(time));
    res *= omega*sin(omega*point[ax] + px);
    return scalar_product(res, normal);
  }
  // g = w sin(w x_ax) (e_ax x uvec) + cos(w x_ax) (rvec - e_ax (e_ax . rvec)), per arm
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    const unsigned int ax = arm_axis(point);
    auto res = cross_e(ax, uvec(time));
    res *= omega*sin(omega*point[ax] + px);
    auto perp = rvec(time);
    perp[ax] = 0.;
    for (unsigned int k = 0; k < 3; k++)
      res[k] += perp[k]*cos(omega*point[ax] + px);
    return scalar_product(res, normal);
  }
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    auto res = uvec(time);
    res *= cos(omega*(point[0]+point[1]+point[2]) + px);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    auto res = dt_uvec(time);
    res *= cos(omega*(point[0]+point[1]+point[2]) + px);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    auto res = dt_rvec(time);
    res *= cos(omega*(point[0]+point[1]+point[2]) + px);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    auto res = rvec(time);
    res *= cos(omega*(point[0]+point[1]+point[2]) + px);
    return res;
  }

  // Dual fields of the manufactured solution (for the dual-variable error). Both are odd
  // under flipping the edge orientation, so the evaluator passes the edge's axial frame
  // vector d = inner_normal(0), orientation included. With sig = x+y+z (the signed arm
  // coordinate on cross2/single1) and (d.1) = sum_k d_k = +-1 on an axis-aligned arm:
  //   n = -du/ds - d x r = w sin(w sig + px) (d.1) uvec - cos(w sig + px) d x rvec
  //   m = -dr/ds         = w sin(w sig + px) (d.1) rvec
  static SmallVec<space_dimT, param_float_t> analytic_result_n(
      const Point<space_dimT, param_float_t>& point,
      const Point<space_dimT, param_float_t>& axial,
      const param_float_t time = 0.)
  {
    const param_float_t sig = point[0] + point[1] + point[2];
    const param_float_t dsum = axial[0] + axial[1] + axial[2];
    auto res = uvec(time);
    res *= omega * sin(omega*sig + px) * dsum;
    const auto r = rvec(time);
    const param_float_t c = cos(omega*sig + px);
    for (unsigned int k = 0; k < 3; k++)
      res[k] -= c * (axial[(k+1)%3]*r[(k+2)%3] - axial[(k+2)%3]*r[(k+1)%3]);
    return res;
  }
  static SmallVec<space_dimT, param_float_t> analytic_result_m(
      const Point<space_dimT, param_float_t>& point,
      const Point<space_dimT, param_float_t>& axial,
      const param_float_t time = 0.)
  {
    const param_float_t sig = point[0] + point[1] + point[2];
    const param_float_t dsum = axial[0] + axial[1] + axial[2];
    auto res = rvec(time);
    res *= omega * sin(omega*sig + px) * dsum;
    return res;
  }
};
