#include <vector>
#include <HyperHDG/global_loop/ch-kp.hxx>
#include <HyperHDG/local_solver/ch-kp.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/sparse_la.hxx>
#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include "../reproducibles_python/parameters/chkp.hxx"

#include <chrono>


int main()
{
  SmallVec<2, unsigned int> top_con(32);
  std::vector<double> def_con={1., 3., 4., 1., 1., 1., 1., 1.,-3., -2.} ;
  GlobalLoop::Nonlinear<Topology::Cubic<2, 2>,
                       Geometry::UnitCube<2, 2, double>,
                       NodeDescriptor::Cubic<2, 2>,
                       LocalSolver::Chkp<2, 1, 3, ChkpParametersAcc1> >
    problem(top_con, def_con);
  std::vector<double> v = problem.make_initial(problem.zero_vector(), 0.);
  const auto start_rf = std::chrono::high_resolution_clock::now();
  std::vector<double> w = problem.residual_flux(v, 0.);
  const auto end_rf = std::chrono::high_resolution_clock::now();  
  std::cout << "Residual flux took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_rf - start_rf).count() << "ms \n";
  /*
  const auto start_ttf = std::chrono::high_resolution_clock::now();
  sparse_mat<std::vector<double>> m = problem.trace_to_flux_mat(v, 0.);
  const auto end_ttf = std::chrono::high_resolution_clock::now();
  std::cout << "Matrix took " << std::chrono::duration_cast<std::chrono::milliseconds>(end_ttf - start_ttf).count() << "ms \n";
  */
  //std::for_each(w.begin(), w.end(), [](double u) {std::cout << u << "\n";});
  //std::cout << "\n";                    
  //std::for_each(v.begin(), v.end(), [](double u) {std::cout << u << "\n";});
  return 0;
}
