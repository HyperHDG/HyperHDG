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
#include "../reproducibles_python/parameters/chkp.hxx"

#include <array>
#include <HyperHDG/dense_la.hxx>


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
  if (bdr == UP)
    std::cout << "v^:\t" << lambda[4] << "\t" << lambda[5] << "\n";
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
  typedef LocalSolver::Chkp<2, 1, 3, ChkpParameters> lst;
  HDGHyperGraph<lst::n_glob_dofs_per_node(),
                Topology::File<2, 2>,
                Geometry::File<2, 2>,
                NodeDescriptor::File<2, 2>,
                lst::data_type>
    hg("domains/square.geo");
  hg.set_refinement(3);
  lst ls ;
  std::array< std::array< double, 6 >, 4> lambda_n, res_flux, dir, out;
  std::vector<double> xv;
  for(unsigned int i = 0; i < hg.n_global_dofs(); i++)
    xv.push_back( (double) i);
  SmallVec<28> coeff(0.), res(0.);
  SmallVec<4, unsigned int> hyEdge_hyNodes;
  std::for_each(hg.begin(), hg.end(), [&](auto he)
      {
        hyEdge_hyNodes = he.topology.get_hyNode_indices();
        for (unsigned int n = 0; n < 4; ++n)
          hg.hyNode_factory().get_dof_values(hyEdge_hyNodes[n], xv, lambda_n[n]);
        ls.make_initial(lambda_n, he);
        /*
        for (unsigned int n = 0; n < 4; ++n)
        {
          lambda_n[n][0] = 1.;
        }
        */
        //ls.make_skeleton(lambda_n, he, 1.);
        std::cout << "local lambda\n";
        for (unsigned int n = 0; n < 4; ++n)
        {
          std::for_each(lambda_n[n].begin(), lambda_n[n].end(), [](auto i){std::cout << i <<"\t";});
          std::cout << "\n";
          std::cout << he.data.uh_old[n];
        }
        std::cout << "u_old:\t";
        std::cout << he.data.u_old;
        coeff[0] = 1.;
        coeff[1] = 2.;
        coeff[2] = 3.;
        coeff[3] = 4.;
        std::cout << ls.get_residual(lambda_n, coeff, res, he, 1.);
        
        std::cout << ls.newton(lambda_n, coeff, he, 1.) << std::endl;
          
        SmallVec<28> res = ls.get_residual(lambda_n, coeff, res, he, 1.);
        std::cout << "Residuen:\n" << res;
        print_coeff(coeff);
        for (unsigned int n = 0; n < 4; ++n)
        {
          std::cout << he.node_descriptor[n] << "\n";
          print_bdr_values(lambda_n[n], coeff, get_bdr(he, n));
        }
        ls.residual_flux(lambda_n, res_flux, he, 1.);
        std::cout << "Kopplungsbeitrag:\n";
        for (unsigned int n = 0; n < 4; ++n)
        {
          std::for_each(res_flux[n].begin(), res_flux[n].end(), [](auto i){std::cout << i <<"\t";});
          std::cout << "\n";
        }
        for (unsigned int n = 0; n < 4; ++n)
        {
          for (unsigned int i = 0; i < 6; ++i)
          {
            dir[n][i] = 1.;
            res_flux[n][i] = lambda_n[n][i] + 0.001;
          }
        }
        
        SmallVec<28> f0, f1;
        res = ls.residual_lambda_directional_derivative(lambda_n, coeff, dir, he, 1.);
        int i = 0;
        //std::for_each(res.begin(), res.end(), [&i](double e) {std::cout << i++ << "\t" << e << "\n";});
        f0 = ls.get_residual(lambda_n, coeff, f0, he, 1.);
        f1 = ls.get_residual(res_flux, coeff, f1, he, 1.);
        f0 = f1 - f0;
        f0 = 1000. * f0;
        i = 0;
        //std::for_each(f0.begin(), f0.end(), [&i](double e) {std::cout << i++ << "\t" << e << "\n";});

        res = res - f0;
        i = 0;
        std::for_each(res.begin(), res.end(), [&i](double e) {std::cout << i++ << "\t" << e << "\n";});
       
        
        
        for (unsigned int n = 0; n < 4; ++n)
        {
          for (unsigned int i = 0; i < 6; ++i)
          {
            dir[n][i] = 0.;
            out[n][i] = 0.;
            res_flux[n][i] = lambda_n[n][i] + 0.001;
          }
        }

        out = ls.residual_flux(lambda_n, out, he, 1.);
        dir = ls.residual_flux(res_flux, dir, he, 1.);
        for(int n = 0; n < 4; ++n)
        {
          for(int i = 0; i < 6; ++i)
          {
            out[n][i] = 1000. * (dir[n][i] - out[n][i]);
            dir[n][i] = 1.;
            res_flux[n][i] = 0.;
          }
        }
        res_flux = ls.trace_to_flux(lambda_n, dir, res_flux, he, 1.);
        for(int n = 0; n < 4; ++n)
          for(int i = 0; i < 6; ++i)
            std::cout << n << " " << i << " " << res_flux[n][i] - out[n][i] << "\n";
        

        
        for(int n = 0; n < 4; ++n)
        {
          for(int i = 0; i < 6; ++i)
          {
            out[n][i] = 0.;
            dir[n][i] = 0.;
            res_flux[n][i] = 0.;
          }
        }
        out = ls.coupling_function(lambda_n, coeff, out, he, 1.);
        dir = ls.coupling_function(lambda_n, coeff + 0.0001, dir, he, 1.);
        for(int n = 0; n < 4; ++n)
        {
          for(int i = 0; i < 6; ++i)
          {
            out[n][i] = 10000. * (dir[n][i] - out[n][i]);
            res_flux[n][i] = 0.;
          }
        }
        SmallVec<28> coeff_dir(1.);
        res_flux = ls.coupling_coeff_directional_derivative(lambda_n, coeff, coeff_dir, res_flux, he, 1.);
        for(int n = 0; n < 4; ++n)
          for(int i = 0; i < 6; ++i)
            std::cout << n << " " << i << " " << res_flux[n][i] - out[n][i] << "\n";
        
          
        std::cout << "\n";
      });
  return 0;
}
