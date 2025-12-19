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

#include <HyperHDG/dense_la.hxx>
#include "../reproducibles_python/parameters/chkp.hxx"

#include <chrono>


int main() 
{
  const unsigned int poly_deg = 1;
  const unsigned int n_sf = (poly_deg + 1) * (poly_deg + 1);
  typedef LocalSolver::Chkp<2, poly_deg, 3 * poly_deg, ChkpParameters> lst;
  HDGHyperGraph<lst::n_glob_dofs_per_node(),
                Topology::File<2, 2>,
                Geometry::File<2, 2>,
                NodeDescriptor::File<2, 2>,
                lst::data_type>
    hg("domains/square.geo");
  hg.set_refinement(4);
  lst ls ;
  std::array< std::array< double, 3 * (poly_deg + 1)>, 4> lambda_n, res_flux, dir, out;
  std::vector<double> xv;
  for(unsigned int i = 0; i < hg.n_global_dofs(); i++)
    xv.push_back( (double) i);
  SmallVec<7 * n_sf> coeff(1.);
  SmallVec<4, unsigned int> hyEdge_hyNodes;
  auto it = hg.begin();
  for(unsigned int i = 0; i < 5; ++i)
    it++;
  auto he = *it;
      {
        hyEdge_hyNodes = he.topology.get_hyNode_indices();
        for (unsigned int n = 0; n < 4; ++n)
          hg.hyNode_factory().get_dof_values(hyEdge_hyNodes[n], xv, lambda_n[n]);
        ls.make_initial(lambda_n, he);
        const auto start_lr = std::chrono::high_resolution_clock::now();
        SmallVec<7 * n_sf> res(0.);
        for (unsigned int i = 0; i < 100; ++i)
          res = ls.get_residual(lambda_n, coeff, res, he, 1.);
        const auto end_lr = std::chrono::high_resolution_clock::now();
        std::cout << "Local residual took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_lr - start_lr).count() / 100. << "ms\n";
        
        const auto start_j = std::chrono::high_resolution_clock::now();
        for (unsigned int i = 0; i < 100; ++i)
          auto jac = ls.jacobi(lambda_n, coeff, he, 1.);
        const auto end_j = std::chrono::high_resolution_clock::now();
        std::cout << "Jacobian took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_j - start_j).count() / 100. << "ms\n";

        const auto start_n = std::chrono::high_resolution_clock::now();
        unsigned int newton_iter = 0;
        for (unsigned int i = 0; i < 100; ++i)
        {
          newton_iter = ls.newton(lambda_n, coeff, he, 1.);
          coeff *= 0;
        }
        const auto end_n = std::chrono::high_resolution_clock::now();
        std::cout << "Newton took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_n - start_n).count() / 100. << "ms for " << newton_iter << " iterations\n";
          
        const auto start_c = std::chrono::high_resolution_clock::now();
        for (unsigned int i = 0; i < 100; ++i)
          ls.coupling_function(lambda_n, coeff, out, he, 1.);
        const auto end_c = std::chrono::high_resolution_clock::now();
        std::cout << "Coupling took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_c - start_c).count() / 100. << "ms\n";

        const auto start_rf = std::chrono::high_resolution_clock::now();
        for (unsigned int i = 0; i < 100; ++i)
          ls.residual_flux(lambda_n, res_flux, he, 1.);
        const auto end_rf = std::chrono::high_resolution_clock::now();
        std::cout << "Residual flux took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_rf - start_rf).count() / 100. << "ms\n";

        const auto start_t = std::chrono::high_resolution_clock::now();
        for (unsigned int i = 0; i < 100; ++i)
          ls.trace_to_flux(lambda_n, out, res_flux, he, 1.);
        const auto end_t = std::chrono::high_resolution_clock::now();
        std::cout << "Trace to flux took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_t - start_t).count() / 100. << "ms\n";


        std::cout << "\n";
      }
  return 0;
}
