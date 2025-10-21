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
  static constexpr param_float_t tau_fr=-.5;
  
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

int main() 
{
  typedef LocalSolver::Chkp<2, 1, 3, ChkpParametersOne> lst;
  HDGHyperGraph<lst::n_glob_dofs_per_node(),
                Topology::File<2, 2>,
                Geometry::File<2, 2>,
                NodeDescriptor::File<2, 2>,
                lst::data_type>
    hg("domains/square.geo");
  hg.set_refinement(2);
  lst ls ;
  std::array< std::array< double, 6 >, 4> lambda_n, res_flux, dir, out;
  std::vector<double> xv;
  for(unsigned int i = 0; i < hg.n_global_dofs(); i++)
//    xv.push_back( (double) i);
    xv.push_back( (double) 0);
  SmallVec<28> coeff(1.);
/*  for (unsigned int i = 0; i < 4; ++i)
    coeff(4 + i, 0) = 0.;*/
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
          lambda_n[n][2] = 1.;
          std::for_each(lambda_n[n].begin(), lambda_n[n].end(), [](auto i){std::cout << i <<"\t";});
          std::cout << "\n";
          std::cout << he.data.uh_old[n];
        }
        std::cout << "u_old:\t";
        std::cout << he.data.u_old;
          
        std::cout << "Jacobi analytisch \n" << ls.jacobi(lambda_n, coeff, he, 0.);
        std::cout << "Jacobi numerisch \n" << ls.jacobi(lambda_n, coeff, he, 0.) - ls.finite(lambda_n, coeff, he, .0000001);

        std::cout << "\n";
      });
  return 0;
}
