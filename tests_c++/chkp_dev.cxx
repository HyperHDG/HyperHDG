#include <vector>
#include <HyperHDG/global_loop/ch-kp.hxx>
#include <HyperHDG/local_solver/ch-kp.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/dense_la.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>

template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParameters
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 8U> dirichlet_nodes{
  1, 2, 3, 4, 5, 6, 7, 8};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>& p,
                               const param_float_t t = 0.)
  {
    return sin(p[0]) * sin(p[1]) * exp(-t);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return sin(p[0]) * sin(p[1]) * exp(-t);
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
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return sin(p[0]) * sin(p[1]) * exp(-t);
  }
  
  static constexpr param_float_t kappa=-1.;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }

};  


int main()
{
  SmallVec<2, unsigned int> top_con;
  top_con[0] = 1;
  top_con[1] = 1;
  std::vector<double> def_con={1., 3., 4., 1., 1., 1., 1., 1.,-3., -2.} ;
  GlobalLoop::Nonlinear<Topology::Cubic<2, 2>,
                       Geometry::UnitCube<2, 2, double>,
                       NodeDescriptor::Cubic<2, 2>,
                       LocalSolver::Chkp<2, 1, 3, ChkpParameters> >
    problem(top_con, def_con);
  std::vector<double> v = problem.make_initial(problem.zero_vector(), 0.);
  std::vector<double> w = problem.residual_flux(v, 0.);
  std::for_each(w.begin(), w.end(), [](double u) {std::cout << u << "\n";});
  std::cout << "\n";
                       
  std::for_each(v.begin(), v.end(), [](double u) {std::cout << u << "\n";});
  return 0;
}
