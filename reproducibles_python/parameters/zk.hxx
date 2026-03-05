#include <array>
#include <HyperHDG/dense_la.hxx>

template <unsigned int space_dimT, typename param_float_t = double>
struct ZKParameters
{
  private:
    static constexpr param_float_t beta = 0.;
    static constexpr param_float_t c = 0.3;
  public:
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 3U> right_nodes{2, 5, 8};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t uh_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return analytic_result(p, t);
  }
  /*!***********************************************************************************************
   * \brief   q values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t qh_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    param_float_t x = p[0], y = p[1];
    param_float_t arg = 0.5 * sqrt(c) * ( (x - c * t) * cos(beta) + y * sin(beta) );
    param_float_t outer = -2. * 3. * c * sinh(arg) / pow(cosh(arg), 3);
    return outer * 0.5 * sqrt(c) * cos(beta);
  }
  /*!***********************************************************************************************
   * \brief   s values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t sh_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    param_float_t x = p[0], y = p[1];
    param_float_t arg = 0.5 * sqrt(c) * ( (x - c * t) * cos(beta) + y * sin(beta) );
    param_float_t outer = -2. * 3. * c * sinh(arg) / pow(cosh(arg), 3);
    return outer * 0.5 * sqrt(c) * sin(beta);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    param_float_t x = p[0], y = p[1];
    param_float_t arg = 0.5 * sqrt(c) * ( (x - c * t) * cos(beta) + y * sin(beta) );
    return 3. * c / pow(cosh(arg), 2);
  }
  
};  


template <unsigned int space_dimT, typename param_float_t = double>
struct ZKParametersOne
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 3U> right_nodes{2, 5, 8};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return 1;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t uh_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return 1;
  }
  /*!***********************************************************************************************
   * \brief   Neumann values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t qh_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   s values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t sh_value(const Point<space_dimT, param_float_t>& p,
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
    return 1;
  }
  
};  