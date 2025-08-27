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
                               const param_float_t = 0.)
  {
    return p[0] * p[0] + p[1];
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& p,
                                       const param_float_t = 0.)
  {
    return p[0] * p[0] + p[1];
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
  static constexpr param_float_t kappa=-1.;
  static param_float_t tau_f(param_float_t arg)
  {
    return 1.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }

};  

int main() {
  typedef LocalSolver::Chkp<2, 1, 3, ChkpParameters> lst;
  SmallVec<2, unsigned int> top_con(3U);
  HDGHyperGraph<lst::n_glob_dofs_per_node(),
                Topology::Cubic<2, 2>,
                Geometry::UnitCube<2, 2>,
                NodeDescriptor::Cubic<2, 2>,
                lst::data_type>
    hg(top_con);
  lst ls(*(hg.begin()));
  std::array< std::array< double, 6 >, 4> lambda_n;
  std::vector<double> xv;
  for(unsigned int i = 0; i < hg.n_global_dofs(); i++)
    xv.push_back( (double) i);
  SmallVec<28> coeff(1.);
  SmallVec<4, unsigned int> hyEdge_hyNodes;
  std::for_each(hg.begin(), hg.end(), [&](auto he)
      {
        hyEdge_hyNodes = he.topology.get_hyNode_indices();
        for (unsigned int n = 0; n < 4; ++n)
          hg.hyNode_factory().get_dof_values(hyEdge_hyNodes[n], xv, lambda_n[n]);
        ls.make_initial_skeleton(lambda_n, he);
        std::cout << "local lambda\n";
        for (unsigned int n = 0; n < 4; ++n)
        {
          std::for_each(lambda_n[n].begin(), lambda_n[n].end(), [](auto i){std::cout << i <<"\t";});
          std::cout << "\n";
        }
        std::cout << "u_old:\n";
        std::cout << he.data.u_old;
        //SmallVec<28> res = ls.get_residual(lambda_n, coeff, he, 0.);
        //int i = 0;
        //std::for_each(res.begin(), res.end(), [&i](double e) {std::cout << i++ << "\t" << e << "\n";});
        //std::cout << "Jacobi analytisch \n" << ls.jacobi(lambda_n, coeff, he, 0.);
        //std::cout << "Jacobi numerisch \n" << ls.jacobi(lambda_n, coeff, he, 0.) - ls.finite(lambda_n, coeff, he, .01);
        //std::cout << ls.newton(lambda_n, coeff, he, 0.);
      });
  return 0;
}
