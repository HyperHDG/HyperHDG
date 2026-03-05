#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <array>
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
    res[1] = 1;
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
    SmallVec<space_dimT, param_float_t> res(1.);
    res[2] = point[0];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
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
    res[2] = -1;
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
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
    res[1] = 1;
    res[2] = 1;
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
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(1.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_n(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[1] = 1;
    res[2] = -1;
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
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
   * \brief   Inverse diffusionbeam_network_bilaplacian.hxx coefficient in PDE as analytic function.
   ************************************************************************************************/
  // static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
  //                                              const param_float_t = 0.)
  // {
  //   return 1.;
  //   // return 1. / M_PI / M_PI;
  // }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] =  -6*time*point[0] - 0*15*time*point[0]*point[0];
    res[1] = -12*time*point[0] + 1*33*time*point[0]*point[0];
    res[2] = -18*time*point[0] - 1*21*time*point[0]*point[0];
    return scalar_product(res, normal);
    // return M_PI * (M_PI - 1.) * cos(M_PI * point[0]) * normal[2] *
    //          (point[1] == 0. && point[2] == 0.) +
    //        M_PI * cos(M_PI * point[1]) * (M_PI * normal[1] + normal[0]) *
    //          (point[0] == 0. && point[2] == 0.);
    // return -M_PI * cos(M_PI * point[0]) * normal[2] * (point[1] == 0. && point[2] == 0.);
    // return M_PI * M_PI * cos(M_PI * point[0]) * normal[2] * (point[1] == 0. && point[2] == 0.);
    // return M_PI * M_PI * sin(M_PI * point[0]) * normal[0];
    // return M_PI * M_PI * sin(M_PI * point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    // auto res_n = initial_n(point, time);
    SmallVec<space_dimT, param_float_t> res;
    res[0] = -30*time*point[0];
    res[1] = -42*time*point[0]+9*time*point[0]*point[0]+ 7*time*point[0]*point[0]*point[0];
    res[2] = -66*time*point[0]-6*time*point[0]*point[0]+11*time*point[0]*point[0]*point[0];
    return scalar_product(res, normal);
    // return (M_PI * M_PI - M_PI + 1.) * sin(M_PI * point[0]) * normal[1] *
    //          (point[1] == 0. && point[2] == 0.) +
    //        (M_PI * M_PI + 1.) * sin(M_PI * point[1]) * normal[2] *
    //          (point[0] == 0. && point[2] == 0.);
    // return  (M_PI * M_PI + 1.) * sin(M_PI * point[0]) * normal[1] * (point[1] == 0. && point[2]
    // == 0.); return  -M_PI * sin(M_PI * point[0]) * normal[1] * (point[1] == 0. && point[2] ==
    // 0.); return M_PI * M_PI * sin(M_PI * point[0]) * normal[0]; return M_PI * M_PI * sin(M_PI *
    // point[0]) * normal[0];
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
    // return cos(M_PI * point[0]) * normal[2] + cos(M_PI * point[1]) * normal[1];
    // return point[0] * normal[0];
    // return sin(M_PI * point[0]) * normal[0];
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
    // return sin(M_PI * point[0]) * normal[1] + sin(M_PI * point[1]) * normal[2];
    // return point[0] * normal[0];
    // return sin(M_PI * point[0]) * normal[0];
  }

  // TODO: construct simple good example
  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] = 1 * time * point[0] * point[0] * point[0];
    res[1] = 2 * time * point[0] * point[0] * point[0];
    res[2] = 3 * time * point[0] * point[0] * point[0];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(.5 / sqrt(.2) * (3*(2*point[0]-1)*(2*point[0]-1)-1)); // ingore this
    res[0] = 1 * point[0] * point[0] * point[0];
    res[1] = 2 * point[0] * point[0] * point[0];
    res[2] = 3 * point[0] * point[0] * point[0];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] =  5 * point[0] * point[0] * point[0];
    res[1] =  7 * point[0] * point[0] * point[0];
    res[2] = 11 * point[0] * point[0] * point[0];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] =  5 * time * point[0] * point[0] * point[0];
    res[1] =  7 * time * point[0] * point[0] * point[0];
    res[2] = 11 * time * point[0] * point[0] * point[0];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_n(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] = -3*time*point[0]*point[0] - 0 *  5 * time*point[0]*point[0]*point[0];
    res[1] = -6*time*point[0]*point[0] + 1 * 11 * time*point[0]*point[0]*point[0];
    res[2] = -9*time*point[0]*point[0] - 1 *  7 * time*point[0]*point[0]*point[0];
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    res[0] = -15*time*point[0]*point[0];
    res[1] = -21*time*point[0]*point[0];
    res[2] = -33*time*point[0]*point[0];
    return res;
  }
};




#endif // PARAMETERS_H
