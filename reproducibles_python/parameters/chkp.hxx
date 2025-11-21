#include <array>
#include <HyperHDG/dense_la.hxx>

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParameters
{
  static constexpr double scale_t = .01;
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  static constexpr std::array<unsigned int, 1U> right_nodes{2};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                                param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                        param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                      param_float_t t = 0.)
  {
    t *= scale_t;
    return cos(p[0]) * sin(p[1]) * exp(-t);
  }

  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                      param_float_t t = 0.)
  {
    return 0;
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& p,
                                      param_float_t t = 0.)
  {
    t *= scale_t;
    param_float_t r = 0, x = p[0], y = p[1];
    r -= 2 * scale_t * exp(-t) * sin(x) * sin(y);
    r += 6 * exp(-2*t)  * sin(x) * cos(x) * sin(y) * sin(y);
    r -= exp(-t) * sin(y);
    return r;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                        param_float_t t = 0.)
  {
    t *= scale_t;
    return sin(p[0]) * sin(p[1]) * exp(-t);
  }
  
  static constexpr param_float_t kappa=-.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }

  static constexpr param_float_t tau_fr = 4.;
};  

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersZero
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return 0;
  }
  
  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }

  static constexpr param_float_t kappa=-.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }
  static constexpr param_float_t tau_fr = 4.;

};  

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersOne
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  static constexpr std::array<unsigned int, 1U> right_nodes{2};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return 1.;
  }
  
  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }

  static constexpr param_float_t kappa=-.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }
  static constexpr param_float_t tau_fr = 4.;

};  

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersPeakon
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  static constexpr std::array<unsigned int, 1U> right_nodes{2};
  static constexpr param_float_t c = 1.;
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return -exp(-abs(p[0] + p[1] - c * t)) * copysign(1.0, p[0] + p[1] - c * t);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return exp(-abs(p[0] + p[1] - c * t));
  }
  
  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }

  static constexpr param_float_t kappa=-.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }
  static constexpr param_float_t tau_fr = 4.;

};  

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersTime
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  static constexpr std::array<unsigned int, 1U> right_nodes{2};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 1.;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return 1. * t;
  }
  
  static constexpr param_float_t kappa=-.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }
  static constexpr param_float_t tau_fr = 4.;

};  

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersTimeNL
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  static constexpr std::array<unsigned int, 1U> right_nodes{2};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 3*p[0] - 3 * t - 2;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return p[0] - t;
  }
  
  static constexpr param_float_t kappa=-.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }
  static constexpr param_float_t tau_fr = 4.;

};  

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersAcc1
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  static constexpr std::array<unsigned int, 1U> right_nodes{2};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    param_float_t x = p[0], y = p[1], r = 0.;
    r -= exp(-t) * pow(x, 3) * pow(1-x, 3) * pow(y, 3) * pow(1-y, 3);
    r += exp(-t) * (6 * x * pow(1-x, 3) -18 * pow(x, 2) * pow(1-x, 2) + 6 * pow(x, 3) * (1-x)) * pow(y, 3) * pow(1-y, 3);
    r -= exp(-t) * (3 * pow(x, 2) * pow(1-x, 3) -3 * pow(x, 3) * pow(1-x, 2)) * pow(y, 3) * pow(1-y, 3);
    r += 3 * exp(-2*t) * pow(x, 3) * pow(1-x, 3) * (3 * pow(x, 2) * pow(1-x, 3) - 3 * pow(x, 3) * pow(1-x, 2)) * pow(y, 6) * pow(1-y, 6);
    r -= 2 * exp(-2*t) * (3 * pow(x, 2) * pow(1-x, 3) - 3 * pow(x, 3) * pow(1-x, 2)) * (6 * x * pow(1-x, 3) -18 * pow(x, 2) * pow(1-x, 2) + 6 * pow(x, 3) * (1-x)) * pow(y, 6) * pow(1-y, 6);
    r -= exp(-2*t) * pow(x, 3)  *pow(1-x, 3) * (6 * pow(1-x, 3) -36 * x * pow(1-x, 2) + 36 * pow(x, 2) * (1-x) - 6 * pow(x, 3)) * pow(y, 6) * pow(1-y, 6);
    r -= exp(-t) * (1./140 - pow(x, 4) / 4. + 3 * pow(x, 5) / 5. - pow(x, 6) / 2. + pow(x, 7) / 7.) * (6 * y * pow(1-y, 3) -18 * pow(y, 2) * pow(1-y, 2) + 6 * pow(y, 3) * (1-y));
    return r;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    param_float_t x = p[0], y = p[1];
    return exp(-t) * pow(x, 3) * pow(1-x, 3) * pow(y, 3) * pow(1-y, 3);
  }
  
  static constexpr param_float_t kappa = 0.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }
  static constexpr param_float_t tau_fr = 4.;

}; 
template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersLinear
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  static constexpr std::array<unsigned int, 1U> right_nodes{2};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return p[0];
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return p[0] - t;
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 1;
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return p[0] - t;
  }
  
  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }

  static constexpr param_float_t kappa = 0.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }
  static constexpr param_float_t tau_fr = 4.;

};  

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersLinearRHS
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
  static constexpr std::array<unsigned int, 1U> right_nodes{2};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return t;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return t;
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }

  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 1.;
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return t;
  }
  
  static param_float_t reference_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }

  static constexpr param_float_t kappa = 0.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }
  static constexpr param_float_t tau_fr = 4.;

};  
