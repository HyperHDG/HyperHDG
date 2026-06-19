#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <array>
#include <petsc.h>
#include <petscviewerhdf5.h>
#include <HyperHDG/dense_la.hxx>

template <unsigned int space_dimT, typename param_float_t = double>
struct ConstantDiffusionParameters
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 1> dirichlet_nodes{ 1 };
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary. -> we ignore
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1.;
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 1.;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
};

template <unsigned int space_dimT, typename param_float_t = double>
struct TestHeat0
{
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1;
  }

  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return 0.5*point[0]*point[0] + time;
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return 0;
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return analytic_result(point, 0);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <unsigned int space_dimT, typename param_float_t = double>
struct TestHeat
{
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return M_PI*M_PI*space_dimT;
  }

  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    param_float_t p = 1;
    for (unsigned int i = 0; i < space_dimT; i++)
      p *= sin(M_PI*point[i]);
    return p * exp(-time);
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return 0;
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return 0;
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return analytic_result(point, 0);
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <unsigned int space_dimT, typename param_float_t = double>
struct TestWave0
{
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1;
  }

  // u = point[0]*point[0] + time*time;
  // v = 2*time;
  // q = -2*point[0];

  // v
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return 2*time;
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return 0;
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t initial_q(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return -2*point[0];
  }

  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <unsigned int space_dimT, typename param_float_t = double>
struct TestWave1
{
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return space_dimT;
  }

  // c = 1/dim
  // q = - c * dx u
  // u =     sin(pi*(sum(x)+t))
  // v =  pi*cos(pi*(sum(x)+t))
  // q = -pi/dim*cos(pi*(sum(x)+t))

  // v
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    param_float_t p = 0;
    for (unsigned int i = 0; i < space_dimT; i++)
      p += point[i];
    return cos(M_PI*(p+time));
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return 0;
  }

  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t initial(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return analytic_result(point, time);
  }
  static param_float_t initial_q(const Point<space_dimT, param_float_t>& point,
                               const param_float_t time = 0.)
  {
    return -analytic_result(point, time)/inverse_diffusion_coeff(point, time);
  }

  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave1
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    if (point[0] != 0) {
      res[1] =  1;
      res[2] = -1;
    }
    if (point[1] != 0) {
      res[0] =  1;
      res[2] = -1;
    }
    res *= 1+time;
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(1);
    res *= 1+point[0]+point[1]+point[2];
    res *= 1+time;
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(1.);
    res *= 1+point[0]+point[1]+point[2];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }
};

// timowave
template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave2
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    if (point[0] != 0) {
      res[1] = time*point[1];
      res[2] = time*point[2];
    } else if (point[1] != 0) {
      res[0] = time*point[0];
      res[2] = time*point[2];
    }
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res = point;
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] = time*point[0];
    res[1] = time*point[1];
    res[2] = time*point[2];
    return res;
  }
};

// timowave
template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave3
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] =  -6*time*point[0] - 0*15*time*point[0]*point[0]
              -6*time*point[1] + 1*33*time*point[1]*point[1];
    res[1] = -12*time*point[0] + 1*33*time*point[0]*point[0]
             -12*time*point[1];
    res[2] = -18*time*point[0] - 1*21*time*point[0]*point[0]
             -18*time*point[1] - 1*15*time*point[1]*point[1];
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res;
    res[0] = -30*time*point[0]
             -30*time*point[1]+9*time*point[1]*point[1]+ 5*time*point[1]*point[1]*point[1];
    res[1] = -42*time*point[0]+9*time*point[0]*point[0]+ 7*time*point[0]*point[0]*point[0]
             -42*time*point[1];
    res[2] = -66*time*point[0]-6*time*point[0]*point[0]+11*time*point[0]*point[0]*point[0]
             -66*time*point[1]-3*time*point[1]*point[1]+11*time*point[1]*point[1]*point[1];
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  // TODO: construct simple good example
  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] = 1 * time * point[0] * point[0] * point[0]
           + 1 * time * point[1] * point[1] * point[1];
    res[1] = 2 * time * point[0] * point[0] * point[0]
           + 2 * time * point[1] * point[1] * point[1];
    res[2] = 3 * time * point[0] * point[0] * point[0]
           + 3 * time * point[1] * point[1] * point[1];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0);
    res[0] = 1 * point[0] * point[0] * point[0]
           + 1 * point[1] * point[1] * point[1];
    res[1] = 2 * point[0] * point[0] * point[0]
           + 2 * point[1] * point[1] * point[1];
    res[2] = 3 * point[0] * point[0] * point[0]
           + 3 * point[1] * point[1] * point[1];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] =  5 * point[0] * point[0] * point[0]
           +  5 * point[1] * point[1] * point[1];
    res[1] =  7 * point[0] * point[0] * point[0]
           +  7 * point[1] * point[1] * point[1];
    res[2] = 11 * point[0] * point[0] * point[0]
           + 11 * point[1] * point[1] * point[1];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] =  5 * time * point[0] * point[0] * point[0]
           +  5 * time * point[1] * point[1] * point[1];
    res[1] =  7 * time * point[0] * point[0] * point[0]
           +  7 * time * point[1] * point[1] * point[1];
    res[2] = 11 * time * point[0] * point[0] * point[0]
           + 11 * time * point[1] * point[1] * point[1];
    return res;
  }
};

// timowave
template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave4
{
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};

  static constexpr param_float_t omega = 2*M_PI;

  // f
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    return scalar_product(res, normal);
  }
  // g
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] = -omega*(sin(omega*point[1]))*cos(omega*time); // NOTE: again! opposite sign as mathematica
    res[1] = -omega*(sin(omega*point[0]))*cos(omega*time);
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
    SmallVec<space_dimT, param_float_t> res(0.);
    res[2] = cos(omega*(point[0]+point[1]+point[2]))*cos(omega*time);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0);
    res[2] = -omega*cos(omega*(point[0]+point[1]+point[2]))*sin(omega*time);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }
};






// timowave
template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave5
{
  static constexpr param_float_t omega = 2*M_PI;

  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  // f
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(-omega*omega*(point[0]+point[1]+point[2])*cos(omega*time));
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  // g
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    if (point[0] != 0) {
      res[1] =  cos(omega*time);
      res[2] = -cos(omega*time);
    }
    if (point[1] != 0) {
      res[0] =  cos(omega*time);
      res[2] = -cos(omega*time);
    }
    if (point[2] != 0) {
      res[0] = -cos(omega*time);
      res[1] =  cos(omega*time);
    }
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(cos(omega*time)*(point[0]+point[1]+point[2]));
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(-omega*(point[0]+point[1]+point[2])*sin(omega*time));
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  // static SmallVec<space_dimT, param_float_t> initial_n(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
  //   SmallVec<space_dimT, param_float_t> res(-cos(omega*time));
  //   return res;
  // }

  // static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
  //   SmallVec<space_dimT, param_float_t> res(0.);
  //   return res;
  // }
};


// timowave
template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave6
{
  static constexpr param_float_t omega = 2.*M_PI;
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  // f
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(-omega*omega*cos(omega*time));
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  // g
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(cos(omega*time));
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(-sin(omega*time));
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_n(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[2] = -time;
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }
};


// timowave
template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave7
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  // f
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(6*time*(point[0]+point[1]));
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  // g
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    if (point[0] != 0) {
      res[1] = 1+time*time*time;
      res[2] = -1-time*time*time;
    }
    if (point[1] != 0) {
      res[0] = 1+time*time*time;
      res[2] = -1-time*time*time;
    }
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res((1+time*time*time)*(point[0]+point[1]));
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(3*time*time*(point[0]+point[1]));
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_n(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(-1-time*time*time);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }
};

// timowave
template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave8
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  // f
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] = 0;
    res[1] = time;
    res[2] = 0;
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  // g
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] = 0;
    res[1] = 1;
    res[2] = time*point[0];
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[2] = point[0];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(1.);
    res[2] = time*point[0];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_n(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[1] = time*point[0];
    res[2] = -1;
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[2] = -time;
    return res;
  }
};


/*!*************************************************************************************************
 * \brief     Timoschenko Network tensile stiffness experiment.
 *
 *            Applies a prescribed tensile strain to the right Dirichlet boundary
 *            of a clamped beam network. The displacement is computed as
 *            `strain * length` and applied in the x-direction.
 *
 *            Both `length` and `strain` are runtime-configurable static members
 *            and must be set after loading the mesh, before the solve.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 * \authors   Joseph Holten, KIT, 2026--
 **************************************************************************************************/
template <unsigned int dim = 3, typename Scalar = double>
struct TimoshenkoStiffness
{
  using Pt = Point<dim, Scalar>;

  /// Global extent of the domain in x-direction. Must be set at runtime after loading the network.
  static inline Scalar length = 0;

  /// Applied tensile strain (dimensionless)
  static inline Scalar strain = .15;

  /// Strain normal component
  static inline unsigned int comp = 2;

  /// Read runtime parameters: domain extent from the mesh file's "/domain" "size" attribute, and
  /// `strain`/`comp` from PETSc options (each falling back to the static defaults above).
  static PetscErrorCode Init(const char* path)
  {
    PetscViewer viewer;
    PetscReal size[3];
    PetscInt comp_ = comp;
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-strain", &strain, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-comp", &comp_, NULL));
    comp = comp_;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, path, FILE_MODE_READ, &viewer));
    PetscCall(PetscViewerHDF5ReadAttribute(viewer, "/domain", "size", PETSC_DOUBLE, NULL, size));
    PetscCall(PetscViewerDestroy(&viewer));
    length = size[0];
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  static Scalar right_hand_side_n(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0;
  }

  static Scalar right_hand_side_m(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar dirichlet_value_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return analytic_result_u(point, normal);
  }

  static Scalar dirichlet_value_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return analytic_result_phi(point, normal);
  }

  static Scalar analytic_result_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return strain * point[0] * (point[0] > .5 * length) * normal[comp];
  }

  // stub
  static Scalar analytic_result_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Pt initial_u(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_v(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_s(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_r(const Pt& point, const Scalar time = 0.) {
    return {};
  }
};

/*!*************************************************************************************************
 * \brief     Timoschenko Network gaussian stiffness experiment.
 *
 *            Applies a prescribed gaussian strain to the Dirichlet boundary
 *            of a clamped beam network.
 *
 *            Both `length` and `strain` are runtime-configurable static members
 *            and must be set after loading the mesh, before the solve.
 *
 * \authors   Joseph Holten, KIT, 2026--
 **************************************************************************************************/
template <unsigned int dim = 3, typename Scalar = double>
struct TimoshenkoClampedConstant
{
  using Pt = Point<dim, Scalar>;

  /// Applied force
  static inline Scalar force = 1;

  /// Read runtime parameters: domain extent from the mesh file's "/domain" "size" attribute, and
  /// `strain`/`comp` from PETSc options (each falling back to the static defaults above).
  static PetscErrorCode Init()
  {
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-force", &force, NULL));
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  static Scalar right_hand_side_n(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return force;
  }

  static Scalar right_hand_side_m(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar dirichlet_value_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return analytic_result_u(point, normal);
  }

  static Scalar dirichlet_value_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return analytic_result_phi(point, normal);
  }

  static Scalar analytic_result_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  // stub
  static Scalar analytic_result_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Pt initial_u(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_v(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_s(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_r(const Pt& point, const Scalar time = 0.) {
    return {};
  }
};



/*!*************************************************************************************************
 * \brief     Timoschenko Network gaussian stiffness experiment.
 *
 *            Applies a prescribed gaussian strain to the Dirichlet boundary
 *            of a clamped beam network.
 *
 *            Both `length` and `strain` are runtime-configurable static members
 *            and must be set after loading the mesh, before the solve.
 *
 * \authors   Joseph Holten, KIT, 2026--
 **************************************************************************************************/
template <unsigned int dim = 3, typename Scalar = double>
struct TimoshenkoGaussian
{
  using Pt = Point<dim, Scalar>;

  /// Global extent of the domain in x-direction. Must be set at runtime after loading the network.
  static inline Scalar length = 0;

  /// Applied tensile strain (dimensionless)
  static inline Scalar strain = .15;

  /// Spatial variance (sigma^2) of the Gaussian bump relative to length
  static inline Scalar std_x = 1;

  /// Read runtime parameters: domain extent from the mesh file's "/domain" "size" attribute, and
  /// `strain`/`comp` from PETSc options (each falling back to the static defaults above).
  static PetscErrorCode Init(const char* path)
  {
    PetscViewer viewer;
    PetscReal size[3];
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-strain", &strain, NULL));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-std_x", &std_x, NULL));
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, path, FILE_MODE_READ, &viewer));
    PetscCall(PetscViewerHDF5ReadAttribute(viewer, "/domain", "size", PETSC_DOUBLE, NULL, size));
    PetscCall(PetscViewerDestroy(&viewer));
    length = size[0];
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  static Scalar right_hand_side_n(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0;
  }

  static Scalar right_hand_side_m(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar dirichlet_value_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return analytic_result_u(point, normal);
  }

  static Scalar dirichlet_value_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return analytic_result_phi(point, normal);
  }

  static Scalar analytic_result_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    Pt center(length*.5);
    Pt r = point - center;
    r[2] = 0;
    Scalar r2 = scalar_product(r,r);
    Scalar s = std_x * length;
    return strain * length * exp(-r2/(2*s*s)) * normal[2];
  }

  // stub
  static Scalar analytic_result_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Pt initial_u(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_v(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_s(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_r(const Pt& point, const Scalar time = 0.) {
    return {};
  }
};


/*!*************************************************************************************************
 * \brief     Timoschenko Network drumhead tap test.
 *
 *            Models a "drumhead tap": a single localized impulse applied at the center of a
 *            clamped beam network. The distributed load `right_hand_side_n` is the product of a
 *            Gaussian bump in space (centered at the domain center) and a Gaussian bump in time
 *            (centered at `tap_time`):
 *
 *              f(x, t) = amplitude
 *                        * exp(-|x - center|^2 / (2 * var_space))
 *                        * exp(-(t - tap_time)^2 / (2 * var_time)) * e_comp
 *
 *            The force points in the `comp` spatial direction. The center is taken at
 *            `length / 2` in every spatial dimension, where `length` is the domain extent shared
 *            with TimoshenkoStiffness (i.e. a `[0, length]^dim` box is assumed).
 *
 *            `length`, `var_space`, `var_time`, `tap_time`, `amplitude` and `comp` are runtime-
 *            configurable static members and must be set after loading the network, before the
 *            solve.
 *
 * \authors   Joseph Holten, KIT, 2026--
 **************************************************************************************************/
template <unsigned int dim, typename Scalar = double>
struct TimoshenkoDrumhead
{
  using Pt = Point<dim, Scalar>;

  /// Global extent of the domain in x-direction. The tap is centered at `length / 2` in each
  /// spatial dimension. Must be set at runtime after loading the network.
  static inline Scalar length = 0;

  /// Spatial variance (sigma^2) of the Gaussian force bump.
  static inline Scalar std_x = 1;

  /// Temporal variance (sigma^2) of the Gaussian force bump.
  static inline Scalar std_t = 1;

  /// Time at which the tap peaks.
  static inline Scalar tap_time = 0;

  /// Peak force amplitude of the tap.
  static inline Scalar amplitude = 1;

  /// Integral of the Gaussian force bump over time and space.
  static inline Scalar energy = 1;

  /// Spatial component the tap force points in.
  static inline unsigned int comp = 2;

  /// Read runtime parameters: domain extent from the mesh file's "/domain" "size" attribute, and
  /// `strain`/`comp` from PETSc options (each falling back to the static defaults above).
  static PetscErrorCode Init(const char* path)
  {
    PetscViewer viewer;
    PetscReal size[3];
    PetscInt comp_ = comp;
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-std_x", &std_x, NULL));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-std_t", &std_t, NULL));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-tap_time", &tap_time, NULL));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-energy", &energy, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-comp", &comp_, NULL));
    comp = comp_;

    amplitude = energy / (std_x*std_x*std_t);

    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, path, FILE_MODE_READ, &viewer));
    PetscCall(PetscViewerHDF5ReadAttribute(viewer, "/domain", "size", PETSC_DOUBLE, NULL, size));
    PetscCall(PetscViewerDestroy(&viewer));
    length = size[0];
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  static Scalar right_hand_side_n(const Pt& point, const Pt& normal, const Scalar time = 0.)
  {
    // Tap center: in-plane at the domain mid-point. The grid lies in the z=0 plane, so the
    // out-of-plane coordinate must be 0 -- a scalar 0.5*length broadcast would sit 0.5 out of
    // plane and the bump exp(-0.25/(2 std_x^2)) would vanish for any sharp std_x.
    Pt center(0.5 * length);
    center[dim - 1] = 0;
    const Scalar r2 = scalar_product(point - center, point - center);
    const Scalar space_bump = std::exp(-r2 / (2 * std_x*std_x));

    const Scalar dt = time - tap_time;
    const Scalar time_bump = std::exp(-dt * dt / (2 * std_t*std_t));

    Pt force(0.);
    force[comp] = amplitude * space_bump * time_bump;
    return scalar_product(force, normal);
  }

  static Scalar right_hand_side_m(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar dirichlet_value_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar dirichlet_value_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar analytic_result_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar analytic_result_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Pt initial_u(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_v(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_s(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_r(const Pt& point, const Scalar time = 0.) {
    return {};
  }
};


/*!*************************************************************************************************
 * \brief     Timoschenko Network sinusoidal displacement at clamped boundary.
 *
 * \authors   Joseph Holten, KIT, 2026--
 **************************************************************************************************/
template <unsigned int dim, typename Scalar = double>
struct TimoshenkoSinClamp
{
  using Pt = Point<dim, Scalar>;

  /// Global extent of the domain in x-direction.
  static inline Scalar length = 1;

  /// Spatial component the displacement is prescribed in.
  static inline unsigned int comp = 2;

  /// Temporal frequency
  static inline Scalar freq = 1;

  /// Displacement strain as fraction of length
  static inline Scalar strain = .10;

  /// Read runtime parameters: domain extent from the mesh file's "/domain" "size" attribute, and
  /// `strain`/`comp` from PETSc options (each falling back to the static defaults above).
  static PetscErrorCode Init(const char* path)
  {
    PetscViewer viewer;
    PetscReal size[3];
    PetscInt comp_ = comp;
    PetscFunctionBeginUser;
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-strain", &strain, NULL));
    PetscCall(PetscOptionsGetReal(NULL, NULL, "-freq", &freq, NULL));
    PetscCall(PetscOptionsGetInt(NULL, NULL, "-comp", &comp_, NULL));
    comp = comp_;
    PetscCall(PetscViewerHDF5Open(PETSC_COMM_WORLD, path, FILE_MODE_READ, &viewer));
    PetscCall(PetscViewerHDF5ReadAttribute(viewer, "/domain", "size", PETSC_DOUBLE, NULL, size));
    PetscCall(PetscViewerDestroy(&viewer));
    length = size[0];
    PetscFunctionReturn(PETSC_SUCCESS);
  }

  static Scalar right_hand_side_n(const Pt& point, const Pt& normal, const Scalar time = 0.)
  {
    return 0.;
  }

  static Scalar right_hand_side_m(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar dirichlet_value_u(const Pt& point, const Pt& normal, const Scalar time = 0.)
  {
    Pt res(0.);
    res[comp] = strain*length*sin(2*M_PI*freq*time);
    return scalar_product(res, normal);
  }

  static Scalar dirichlet_value_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar analytic_result_u(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Scalar analytic_result_phi(const Pt& point, const Pt& normal, const Scalar = 0.)
  {
    return 0.;
  }

  static Pt initial_u(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_v(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_s(const Pt& point, const Scalar time = 0.) {
    return {};
  }

  static Pt initial_r(const Pt& point, const Scalar time = 0.) {
    return {};
  }
};






#endif // PARAMETERS_H
