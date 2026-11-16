#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>

#include <array>
#include <cmath>

/*!*************************************************************************************************
 * \file    timoshenko.hxx
 * \brief   Test problems for the Timoshenko beam solver on networks of straight beams.
 *
 * From the former reproducibles_python/parameters/timoshenko.hxx (notation of Hauck, Malqvist,
 *Rupp, "Arbitrary order approximations at constant cost for Timoshenko beam network models" (2024),
 * eq. (2.2)): on an edge with unit tangent i, the displacement u, rotation r, stress resultant n
 * and moment m satisfy, with unit coefficients,
 *   -(d_x u + i x r) = n,   -d_x r = m,   d_x n = f,   d_x m + i x n = g,
 * with continuity of u, r and balance of n, m at the nodes. The rotation variable phi of the solver
 * is -r, right_hand_side_n is f and right_hand_side_m is -g. The parameter functions receive a
 *point and a direction and return the component of the vector field along it.
 **************************************************************************************************/
namespace tests
{
/*!*************************************************************************************************
 * \brief   Manufactured solution on the cross of beams along the x and the y axis
 *          (domains/cross.geo; also the single beam along y of domains/line.geo).
 *
 * u = cos(pi x) e_z + cos(pi y) e_y and r = -(sin(pi x) e_y + sin(pi y) e_z), zero at the four
 * ends. Beams along x (i = e_x): n = (pi + 1) sin(pi x) e_z, m = pi cos(pi x) e_y, so
 * f = pi (pi + 1) cos(pi x) e_z and g = -(pi^2 + pi + 1) sin(pi x) e_y. Beams along y
 * (i = e_y): n = sin(pi y) (pi e_y + e_x), m = pi cos(pi y) e_z, f = pi cos(pi y) (pi e_y + e_x)
 * and g = -(pi^2 + 1) sin(pi y) e_z. Both beam lines are smooth through the junction, so the
 * junction conditions hold.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TimoshenkoCrossParameters
{
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  using P = Point<space_dimT, param_float_t>;

  static param_float_t right_hand_side_n(const P& point, const P& normal, const param_float_t = 0.)
  {
    return M_PI * (M_PI + 1.) * cos(M_PI * point[0]) * normal[2] *
             (point[1] == 0. && point[2] == 0.) +
           M_PI * cos(M_PI * point[1]) * (M_PI * normal[1] + normal[0]) *
             (point[0] == 0. && point[2] == 0.);
  }
  static param_float_t right_hand_side_m(const P& point, const P& normal, const param_float_t = 0.)
  {
    return (M_PI * M_PI + M_PI + 1.) * sin(M_PI * point[0]) * normal[1] *
             (point[1] == 0. && point[2] == 0.) +
           (M_PI * M_PI + 1.) * sin(M_PI * point[1]) * normal[2] *
             (point[0] == 0. && point[2] == 0.);
  }
  static param_float_t analytic_result_u(const P& point, const P& normal, const param_float_t = 0.)
  {
    return cos(M_PI * point[0]) * normal[2] + cos(M_PI * point[1]) * normal[1];
  }
  static param_float_t analytic_result_phi(const P& point,
                                           const P& normal,
                                           const param_float_t = 0.)
  {
    return sin(M_PI * point[0]) * normal[1] + sin(M_PI * point[1]) * normal[2];
  }
  static param_float_t dirichlet_value_u(const P& p, const P& n, const param_float_t t = 0.)
  {
    return analytic_result_u(p, n, t);
  }
  static param_float_t dirichlet_value_phi(const P& p, const P& n, const param_float_t t = 0.)
  {
    return analytic_result_phi(p, n, t);
  }
};
}  // end of namespace tests
