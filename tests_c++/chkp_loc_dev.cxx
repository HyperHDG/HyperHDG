#include <HyperHDG/local_solver/ch-kp.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/dense_la.hxx>
#include <algorithm>
#include <array>
#include <vector>
#include <cmath>
#include <iostream>

#include <array>
#include <HyperHDG/dense_la.hxx>


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
  static constexpr std::array<unsigned int, 2U> neumann_nodes{1, 2};
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
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>& p,
                                     const param_float_t t = 0.)
  {
    return cos(p[0]) * sin(p[1]) * exp(-t);
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
  
  static constexpr param_float_t kappa=-.5;
  
  static param_float_t tau_f(param_float_t arg)
  {
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }

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
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t t = 0.)
  {
    return 1;
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
    return 1;
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

enum Bdr_type {UP, DOWN, LEFT, RIGHT, UNDEFINED};

template <typename hyEdgeT>
Bdr_type get_bdr(hyEdgeT& he, const unsigned int bdr)
{
  SmallVec<2> loc_normal = he.geometry.local_normal(bdr);
  const double eps = ldexp(1, -40);
  if (loc_normal[1] * loc_normal[1] < eps)
  {
    if (loc_normal[0] > 0)
      return RIGHT;
    else
      return LEFT;
  }
  else if (loc_normal[0] * loc_normal[0] < eps)
  {
    if (loc_normal[1] > 0)
      return UP;
    else
      return DOWN;
  }
  else
  {
    return UNDEFINED;
  }
}

template <typename CVecT, typename LVecT>
void print_bdr_values(const LVecT& lambda, const CVecT& coeff, const Bdr_type bdr)
{
  unsigned int sa, si;
  int f;
  const double r3 = sqrt(3);
  switch (bdr)
  {
    case UP:
      std::cout << "Upper boundary\n";
      sa = 1;
      si = 2;
      f = 1;
      break;
    case DOWN:
      std::cout << "Lower boundary\n";
      sa = 1;
      si = 2;
      f = -1;
      break;
    case LEFT:
      std::cout << "Left boundary\n";
      sa = 2;
      si = 1;
      f = -1;
      break;
    case RIGHT:
      std::cout << "Right boundary\n";
      sa = 2;
      si = 1;
      f = 1;
      break;
    default:
      std::cout << "Undefined boundary!\n";
      return;
  }
  std::cout << "u:\t" << coeff[0] + f * r3 * coeff[si] << "\t" << coeff[sa] + f * r3 * coeff[sa + si] << "\n";
  std::cout << "u^:\t" << lambda[0] << "\t" << lambda[1] << "\n";
  if (bdr == UP or bdr == DOWN)
    return;
  std::cout << "q^:\t" << lambda[2] << "\t" << lambda[3] << "\n";
  if (bdr == RIGHT)
    std::cout << "v^:\t" << lambda[4] << "\t" << lambda[5] << "\n";
  std::cout << "q:\t" << coeff[4] + f * r3 * coeff[4 + si] << "\t" << coeff[4 + sa] + f * r3 * coeff[4 + sa + si] << "\n";
  std::cout << "p:\t" << coeff[2 * 4] + f * r3 * coeff[2 * 4 + si] << "\t" << coeff[2 * 4 + sa] + f * r3 * coeff[2 * 4 + sa + si] << "\n";
  std::cout << "v:\t" << coeff[4 * 4] + f * r3 * coeff[4 * 4 + si] << "\t" << coeff[4 * 4 + sa] + f * r3 * coeff[4 * 4 + sa + si] << "\n";
  std::cout << "z:\t" << coeff[5 * 4] + f * r3 * coeff[5 * 4 + si] << "\t" << coeff[5 * 4 + sa] + f * r3 * coeff[5 * 4 + sa + si] << "\n";
  return;
}

template <typename VecT>
void print_coeff(const VecT& coeff)
{
  std::cout << "u:\t" << coeff[0 * 4 + 0] << "\t" << coeff[0 * 4 + 1] << "\t" << coeff[0 * 4 + 2] << "\t" << coeff[0 * 4 + 3] << "\n";
  std::cout << "q:\t" << coeff[1 * 4 + 0] << "\t" << coeff[1 * 4 + 1] << "\t" << coeff[1 * 4 + 2] << "\t" << coeff[1 * 4 + 3] << "\n";
  std::cout << "p:\t" << coeff[2 * 4 + 0] << "\t" << coeff[2 * 4 + 1] << "\t" << coeff[2 * 4 + 2] << "\t" << coeff[2 * 4 + 3] << "\n";
  std::cout << "s:\t" << coeff[3 * 4 + 0] << "\t" << coeff[3 * 4 + 1] << "\t" << coeff[3 * 4 + 2] << "\t" << coeff[3 * 4 + 3] << "\n";
  std::cout << "v:\t" << coeff[4 * 4 + 0] << "\t" << coeff[4 * 4 + 1] << "\t" << coeff[4 * 4 + 2] << "\t" << coeff[4 * 4 + 3] << "\n";
  std::cout << "z:\t" << coeff[5 * 4 + 0] << "\t" << coeff[5 * 4 + 1] << "\t" << coeff[5 * 4 + 2] << "\t" << coeff[5 * 4 + 3] << "\n";
  std::cout << "r:\t" << coeff[6 * 4 + 0] << "\t" << coeff[6 * 4 + 1] << "\t" << coeff[6 * 4 + 2] << "\t" << coeff[6 * 4 + 3] << "\n";
  return;
}

int main() 
{
  typedef LocalSolver::Chkp<2, 1, 3, ChkpParametersOne> lst;
  HDGHyperGraph<lst::n_glob_dofs_per_node(),
                Topology::File<2, 2>,
                Geometry::File<2, 2>,
                NodeDescriptor::File<2, 2>,
                lst::data_type>
    hg("domains/square.geo");
  hg.set_refinement(1);
  lst ls ;
  std::array< std::array< double, 6 >, 4> lambda_n, res_flux;
  std::vector<double> xv;
  for(unsigned int i = 0; i < hg.n_global_dofs(); i++)
    xv.push_back( (double) i);
  SmallVec<28> coeff(0.);
  SmallVec<4, unsigned int> hyEdge_hyNodes;
  std::for_each(hg.begin(), hg.end(), [&](auto he)
      {
        hyEdge_hyNodes = he.topology.get_hyNode_indices();
        for (unsigned int n = 0; n < 4; ++n)
          hg.hyNode_factory().get_dof_values(hyEdge_hyNodes[n], xv, lambda_n[n]);
        ls.make_initial(lambda_n, he);
        ls.make_skeleton(lambda_n, he, 1.);
        std::cout << "local lambda\n";
        for (unsigned int n = 0; n < 4; ++n)
        {
          std::for_each(lambda_n[n].begin(), lambda_n[n].end(), [](auto i){std::cout << i <<"\t";});
          std::cout << "\n";
          std::cout << he.data.uh_old[n];
        }
        std::cout << "u_old:\t";
        std::cout << he.data.u_old;
/*          
        std::cout << ls.newton(lambda_n, coeff, he, 0.) << std::endl;
        ls.newton(lambda_n, coeff, he, 1.);
        coeff[0] = 1.;
        SmallVec<28> res = ls.get_residual(lambda_n, coeff, he, 1.);
        std::cout << "Residuen:\n" << res;
        print_coeff(coeff);
        for (unsigned int n = 0; n < 4; ++n)
          print_bdr_values(lambda_n[n], coeff, get_bdr(he, n));
        ls.residual_flux(lambda_n, res_flux, he, 1.);
        std::cout << "Kopplungsbeitrag:\n";
        for (unsigned int n = 0; n < 4; ++n)
        {
          std::for_each(res_flux[n].begin(), res_flux[n].end(), [](auto i){std::cout << i <<"\t";});
          std::cout << "\n";
        }
        //int i = 0;
        //std::for_each(res.begin(), res.end(), [&i](double e) {std::cout << i++ << "\t" << e << "\n";});
  */      
        std::cout << "Jacobi analytisch \n" << ls.jacobi(lambda_n, coeff, he, 0.);
        std::cout << "Jacobi numerisch \n" << ls.jacobi(lambda_n, coeff, he, 0.) - ls.finite(lambda_n, coeff, he, .0001);
        std::cout << "\n";
      });
  return 0;
}
