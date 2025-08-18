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

int main() {
  typedef LocalSolver::Chkp<2, 1, 3> lst;
  SmallVec<2, unsigned int> top_con(1U);
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
        ls.make_initial(he);
        hyEdge_hyNodes = he.topology.get_hyNode_indices();
        for (unsigned int n = 0; n < 4; ++n)
        {
          hg.hyNode_factory().get_dof_values(hyEdge_hyNodes[n], xv, lambda_n[n]);
          std::for_each(lambda_n[n].begin(), lambda_n[n].end(), [](auto i){std::cout << i <<"\t";});
          std::cout << "\n";
        }
        SmallVec<28> res = ls.get_residual(lambda_n, coeff, he, 0.);
        int i = 0;
        std::for_each(res.begin(), res.end(), [&i](double e) {std::cout << i++ << "\t" << e << "\n";});
        //std::cout << "Jacobi analytisch \n" << ls.jacobi(lambda_n, coeff, he, 0.);
        //std::cout << "Jacobi numerisch \n" << ls.jacobi(lambda_n, coeff, he, 0.) - ls.finite(lambda_n, coeff, he, .01);
        std::cout << ls.newton(lambda_n, coeff, he, 0.);
      });
  return 0;
}
