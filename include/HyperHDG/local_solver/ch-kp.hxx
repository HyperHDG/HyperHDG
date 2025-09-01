#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

#include <algorithm>
#include <vector>
#include <tuple>
#include <iostream>
#include <cmath>

namespace LocalSolver
{
/*!*************************************************************************************************
 * \brief   Default parameters for the CH-KP equation.
 *
 * \authors   Ruben Gutendorf, Saarland University, 2025.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct ChkpParametersDefault
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> dirichlet_nodes{};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Inverse diffusion coefficient in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t initial(const Point<space_dimT, param_float_t>&,
                               const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
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
    return 4.;
  }
  static param_float_t tau_df(param_float_t arg)
  {
    return 0.;
  }

};  

/*!*************************************************************************************************
 * \brief   Local solver for parabolic diffusion equation on hypergraph.
 *
 * \note    Theta must be one, since only implicit Euler has been implemented at the moment.
 *
 * This class contains the local solver for an isotropic diffusion equation, i.e.,
 * \f[
 *  \partial_t u - \nabla \cdot ( d \nabla u ) = f \quad \text{ in } \Omega, \qquad
 *  u = u_\text D \quad \text{ on } \partial \Omega_\text D, \qquad
 *  - d \nabla u \cdot \nu = g_\text N \quad \text{ on } \partial \Omega_\text N
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dimT = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dimT = 2) or volumes (\c hyEdge_dimT = 3) or hypervolumes (in
 * case of \c hyEdge_dimT > 3). \f$f\f$ and \f$d\f$ are scalars defined in the whole domain, the
 * Dirichlet and Neumann boundary data needs to be defined on their respective hypernodes.
 *
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 * \tparam  parametersT   Struct depending on templates \c space_dimTP and \c lSol_float_TP that
 *                        contains static parameter functions.
 *                        Defaults to above functions included in \c DiffusionParametersDefault.
 * \tparam  lSol_float_t  The floating point type calculations are executed in. Defaults to double.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT = ChkpParametersDefault,
          typename lSol_float_t = double>
class Chkp
{
 public:
  // -----------------------------------------------------------------------------------------------
  // Public, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Dimension of hyper edge type that this object solves on.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Evaluate amount of global degrees of freedom per hypernode.
   *
   * This number must be equal to HyperNodeFactory::n_dofs_per_node() of the HyperNodeFactory
   * cooperating with this object.
   *
   * \retval  n_dofs        Number of global degrees of freedom per hypernode.
   ************************************************************************************************/
  static constexpr unsigned int n_glob_dofs_per_node()
  {
    return 3 * Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
 private:
  // -----------------------------------------------------------------------------------------------
  // Private, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Number of local shape functions (with respect to all spatial dimensions).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1) / 3;
  /*!***********************************************************************************************
   * \brief   Number of local  shape functions (with respect to a face / hypernode).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node() / 3;
  /*!***********************************************************************************************
   * \brief   Number of (local) degrees of freedom per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_dofs_ = 7 * n_shape_fct_;
  /*!***********************************************************************************************
   * \brief   Find out whether a node is of Dirichlet type.
   ************************************************************************************************/
  template <typename parameters>
  static constexpr bool is_dirichlet(const unsigned int node_type)
  {
    return std::find(parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                     node_type) != parameters::dirichlet_nodes.end();
  }

  // -----------------------------------------------------------------------------------------------
  // Private, const members: Parameters and auxiliaries that help assembling matrices, etc.
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Time step size.
   ************************************************************************************************/
  const lSol_float_t delta_t_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_ppu_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_mpu_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_mpv_;

  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_pzu_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_mzu_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_mzv_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_pvu_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_uqq_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_yvu_;
  /*!***********************************************************************************************
   * \brief   Parameter theta that defines the one-step theta scheme.
   ************************************************************************************************/
  const lSol_float_t theta_ = 1.;
  /*!***********************************************************************************************
   * \brief   An integrator helps to easily evaluate integrals (e.g. via quadrature).
   ************************************************************************************************/
  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT> >,
    lSol_float_t>
    integrator;

   
  public:
  // -----------------------------------------------------------------------------------------------
  // Public functions (and one typedef) to be utilized by external functions.
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  struct data_type
  {
    SmallVec<n_shape_fct_, lSol_float_t> u_old = SmallVec<n_shape_fct_, lSol_float_t>(0.);
    std::array<SmallVec<n_shape_bdr_, lSol_float_t>, 2 * hyEdge_dim()> uh_old;
  };
  /*!***********************************************************************************************
   *  \brief  Define type of node elements, especially with respect to nodal shape functions.
   ************************************************************************************************/
  struct node_element
  {
    typedef std::tuple<TPP::ShapeFunction<
      TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT - 1> > >
      functions;
  };
  /*!***********************************************************************************************
   *  \brief  Define how errors are evaluated.
   ************************************************************************************************/
  struct error_def
  {
    /*!*********************************************************************************************
     *  \brief  Define the typename returned by function errors.
     **********************************************************************************************/
    typedef std::array<lSol_float_t, 1U> error_t;
    /*!*********************************************************************************************
     *  \brief  Define how initial error is generated.
     **********************************************************************************************/
    static error_t initial_error()
    {
      std::array<lSol_float_t, 1U> summed_error;
      summed_error.fill(0.);
      return summed_error;
    }
    /*!*********************************************************************************************
     *  \brief  Define how local errors should be accumulated.
     **********************************************************************************************/
    static error_t sum_error(error_t& summed_error, const error_t& new_error)
    {
      for (unsigned int k = 0; k < summed_error.size(); ++k)
        summed_error[k] += new_error[k];
      return summed_error;
    }
    /*!*********************************************************************************************
     *  \brief  Define how global errors should be postprocessed.
     **********************************************************************************************/
    static error_t postprocess_error(error_t& summed_error)
    {
      for (unsigned int k = 0; k < summed_error.size(); ++k)
        summed_error[k] = std::sqrt(summed_error[k]);
      return summed_error;
    }
  };
   /*!***********************************************************************************************
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef std::vector<lSol_float_t> constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   constru       Constructor object.
   ************************************************************************************************/
  Chkp(const constructor_value_type& constru = std::vector<double>({1., 3., 4., 1., 1., 1., 1., 1., -3., -2.}))
  : delta_t_(constru[0]), tau_ppu_(constru[1]), tau_mpu_(constru[2]), tau_mpv_(constru[3]), 
  tau_pzu_(constru[4]), tau_mzu_(constru[5]), tau_mzv_(constru[6]), tau_pvu_(constru[7]), 
  tau_uqq_(constru[8]), tau_yvu_(constru[9])
  {
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallSquareMat<n_loc_dofs_, lSol_float_t> jacobi(const SmallMatT& lambda_values,
                                                          const SmallVec<n_loc_dofs_, lSol_float_t> ca,
                                                           hyEdgeT& hyper_edge,
                                                          const lSol_float_t time) const
  {
  SmallSquareMat<n_shape_fct_, lSol_float_t> mass, mdx, mdy;
  std::array<SmallMat<n_shape_fct_, n_shape_bdr_, lSol_float_t>, 2 * hyEdge_dim()> bdr_sh_sk;
  std::array<SmallSquareMat<n_shape_fct_, lSol_float_t>, 2 * hyEdge_dim()> bdr_sh_sh;
  std::array<SmallVec<hyEdge_dim(), lSol_float_t>, 2 * hyEdge_dim()> loc_normal;
  std::array<SmallSquareMat<n_shape_fct_, lSol_float_t>, n_shape_fct_> tv_shx_sh_sh;
  std::array<std::array<SmallSquareMat<n_shape_fct_, lSol_float_t>, 2 * hyEdge_dim()>, n_shape_fct_> tb_sh_sh_sh;
  std::array<std::array<SmallMat<n_shape_fct_, n_shape_bdr_, lSol_float_t>, 2 * hyEdge_dim()>, n_shape_fct_> tb_sh_sh_sk;
  std::array<std::array<SmallSquareMat<n_shape_bdr_, lSol_float_t>, 2 * hyEdge_dim()>, n_shape_fct_> tb_sh_sk_sk;
 {
   for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j) 
      {
        mass(i, j) = integrator::template integrate_vol_phiphi<decltype(hyEdgeT::geometry)>(
            i, j, hyper_edge.geometry);
        mdx(i, j) = integrator::template integrate_vol_phiDphi<decltype(hyEdgeT::geometry)>(
            j, i, 0, hyper_edge.geometry);
        mdy(i, j) = integrator::template integrate_vol_phiDphi<decltype(hyEdgeT::geometry)>(
            j, i, 1, hyper_edge.geometry);
        for (unsigned int k = 0; k < n_shape_fct_; ++k)
        {
          tv_shx_sh_sh[i].operator()(j, k) = integrator::template integrate_vol_phiphiDphi<decltype(hyEdgeT::geometry)>(
            j, k, i, 0, hyper_edge.geometry);
        }
      }
    }
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr) 
    {
      loc_normal[bdr] = hyper_edge.geometry.local_normal(bdr);
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        for (unsigned int j = 0; j < n_shape_bdr_; j++) 
        {
          bdr_sh_sk[bdr].operator()(i, j) = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
            i, j, bdr, hyper_edge.geometry);
          for (unsigned int k = 0; k < n_shape_fct_; ++k)
          {
            tb_sh_sh_sk[i][bdr].operator()(k, j) = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
              i, k, j, bdr, hyper_edge.geometry);
          }
          for (unsigned int k = 0; k < n_shape_bdr_; ++k)
          {
            tb_sh_sk_sk[i][bdr].operator()(k, j) = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
              i, k, j, bdr, hyper_edge.geometry);
          }
        }
        for (unsigned int j = 0; j < n_shape_fct_; j++) 
        {
          bdr_sh_sh[bdr].operator()(i, j) = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
            i, j, bdr, hyper_edge.geometry);
          for (unsigned int k = 0; k < n_shape_fct_; ++k)
          {
            tb_sh_sh_sh[i][bdr].operator()(k, j) = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, k, j, bdr, hyper_edge.geometry);
          }
        }
      }
    }

  }
    //rearrange coefficients
    std::array<SmallVec<n_shape_fct_, lSol_float_t>, 7> coeff;
    for (unsigned int i = 0; i < 7; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
        coeff[i](j, 0) = ca(i * n_shape_fct_ + j, 0);
    }
    SmallSquareMat<n_loc_dofs_, lSol_float_t> ret(0.);
    using parameters = parametersT<hyEdge_dim(), lSol_float_t>;
    //rearrange lambda values
    std::array<SmallVec<n_shape_bdr_, lSol_float_t>, 2 * hyEdge_dim()> u_hat, q_hat, v_hat;
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr) 
    {
      for (unsigned int i = 0; i < n_shape_bdr_; i++)
      {
        u_hat[bdr][i] = lambda_values[bdr][i];
        q_hat[bdr][i] = lambda_values[bdr][n_shape_bdr_ + i];
        v_hat[bdr][i] = lambda_values[bdr][2 * n_shape_bdr_ + i];
      }
    }

    const lSol_float_t eps = ldexp(1., -40);

    SmallVec<n_shape_fct_, lSol_float_t> helper;
    //first equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        ret(i, n_shape_fct_ + j) = mass(i, j);
        ret(i, j) = mdx(i, j);
      }
    }

    //second equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        ret(n_shape_fct_ + i, 2 * n_shape_fct_ + j) = mass(i, j);
        ret(n_shape_fct_ + i, j) = (tv_shx_sh_sh[i] * coeff[1])(j, 0);
        ret(n_shape_fct_ + i, n_shape_fct_ + j) = (tv_shx_sh_sh[i] * coeff[0])(j, 0);
        for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
        {
          ret(n_shape_fct_ + i, j) -= (0.5 * tb_sh_sh_sh[i][bdr] * coeff[1] 
              + 0.5 * tb_sh_sh_sk[i][bdr] * q_hat[bdr])(j, 0) * loc_normal[bdr][0];
          ret(n_shape_fct_ + i, n_shape_fct_ + j) -= (0.5 * tb_sh_sh_sh[i][bdr] * coeff[0])(j, 0)
           * loc_normal[bdr][0]; 
          ret(n_shape_fct_ + i, n_shape_fct_ + j) += tau_uqq_ * bdr_sh_sh[bdr](i, j)
            * loc_normal[bdr][0] * loc_normal[bdr][0];
        }
      }
    }

    //third equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        ret(2 * n_shape_fct_ + i, 3 * n_shape_fct_ + j) = mass(i, j);
        ret(2 * n_shape_fct_ + i, j) = mdy(i, j);
      }
    }

    //fourth equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        ret(3 * n_shape_fct_ + i, 3 * n_shape_fct_ + j) = mass(i, j);
        ret(3 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) = mdx(i, j);
        for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
        {
          if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] < 0)
          {
            ret(3 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) -= bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0];
            ret(3 * n_shape_fct_ + i, j) = tau_pvu_ * bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0] * loc_normal[bdr][0];
          }
        }
      }
    }

    //fifth equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        ret(4 * n_shape_fct_ + i, 5 * n_shape_fct_ + j) = mass(i, j);
        ret(4 * n_shape_fct_ + i, 6 * n_shape_fct_ + j) = mdx(i, j);
      }
    }

    //sixth equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        ret(5 * n_shape_fct_ + i, 6 * n_shape_fct_ + j) = mass(i, j);
        ret(5 * n_shape_fct_ + i, 5 * n_shape_fct_ + j) = mdx(i, j);
        //third product
        //f-part
        ret(5 * n_shape_fct_ + i, j) = 2 * parameters::kappa * mdx(i, j);
        ret(5 * n_shape_fct_ + i, j) += 3 * (tv_shx_sh_sh[i] * coeff[0])(j, 0);
        //p-part
        ret(5 * n_shape_fct_ + i, 2 * n_shape_fct_ + j) = -mdx(i, j);
        //q-part
        ret(5 * n_shape_fct_ +i, n_shape_fct_ + j) = (tv_shx_sh_sh[i] * coeff[1])(j, 0);
       
        //fourth product
        for(unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
        {
          //z^ and p^
          if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] < 0) //V+
          {
            ret(5 * n_shape_fct_ + i, 5 * n_shape_fct_ + j) -= bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, 2 * n_shape_fct_ + j) -= -bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, j) -= -tau_pzu_ * bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0] * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, j) -= tau_ppu_ * bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0] * loc_normal[bdr][0];
          }
          else if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] > 0) //V+
          {
            ret(5 * n_shape_fct_ + i, 5 * n_shape_fct_ + j) -= bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, 2 * n_shape_fct_ + j) -= -bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, j) -= -tau_mzu_ * bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0] * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, j) -= tau_mpu_ * bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0] * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) -= -tau_mzv_ * bdr_sh_sh[bdr](i, j) 
              * loc_normal[bdr][0] * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) -= tau_mpv_ * bdr_sh_sh[bdr](i, j) 
              * loc_normal[bdr][0] * loc_normal[bdr][0];
          }
          //f^
          if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps) //on V
          {
            //f
            ret(5 * n_shape_fct_ + i, j) -= 2 * parameters::kappa * bdr_sh_sh[bdr](i, j) * loc_normal[bdr][0];
            ret(5 * n_shape_fct_ + i, j) -= 3 * (tb_sh_sh_sh[i][bdr] * coeff[0])(j, 0) * loc_normal[bdr][0];
            //tau_f      
            std::array<lSol_float_t, n_shape_fct_> u;
            for(unsigned int j = 0; j < n_shape_fct_; ++j)
              u[j] = coeff[0][j];
            std::array<lSol_float_t, n_shape_bdr_> uh;
            for(unsigned int j = 0; j < n_shape_bdr_; ++j)
              uh[j] = u_hat[bdr][j];
            ret(5 * n_shape_fct_ + i, j) -= integrate_bdr_phiphicompfun<decltype(hyEdgeT::geometry), parameters::tau_df>
              (i, j, uh, u, bdr, hyper_edge.geometry) * loc_normal[bdr][0] * loc_normal[bdr][0];
          }
        }

        //(v, phi_y)
        ret(5 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) = mdy(i, j);

        //<v^ n_y, phi>
        for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
        {
          ret(5 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) -= bdr_sh_sh[bdr](i, j) * loc_normal[bdr][1];
          ret(5 * n_shape_fct_ + i, j) -= -tau_yvu_ * bdr_sh_sh[bdr](i, j) * loc_normal[bdr][1] * loc_normal[bdr][1];
        }
      }
    }

    //seventh equation
    for(unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        ret(6 * n_shape_fct_ + i, j) = (1. / delta_t_) * mass(i, j);
        ret(6 * n_shape_fct_ + i, 6 * n_shape_fct_ + j) =mass(i, j);
      }
    }
   
    return ret;
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallSquareMat<n_loc_dofs_, lSol_float_t> finite(const SmallMatT& lambda_values,
                                                          const SmallVec<n_loc_dofs_, lSol_float_t> ca,
                                                           hyEdgeT& hyper_edge,
                                                          const lSol_float_t h) const
  {
    SmallSquareMat<n_loc_dofs_, lSol_float_t> res;
    SmallVec<n_loc_dofs_, lSol_float_t> cph, helper;
    for (int j = 0; j < n_loc_dofs_; j++) {
      cph = ca;
      cph(j, 0) += h;
      helper = (1. / h) * (get_residual(lambda_values, cph, hyper_edge, 0.) 
          - get_residual(lambda_values, ca, hyper_edge, 0.));
      for (int i = 0; i <n_loc_dofs_; i++)
        res(i, j) = helper(i, 0);
    }
    return res;
  }

  /*!***********************************************************************************************
   * \brief   Solve local problem (with right-hand side from skeletal).
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatT     The data type of the \c lambda_values.
   * \param   lambda_values Encodes uh, qh, vh.
   * \param   coeff         Coefficients of u, q, p, s, v, z, r to be tested.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Point of time the problem is solved.
   * \retval  residual      Residual (should be zero).
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> get_residual(const SmallMatT& lambda_values,
                                                            const SmallVec<n_loc_dofs_, lSol_float_t> ca,
                                                            hyEdgeT& hyper_edge,
                                                            const lSol_float_t time) const
  {
    static_assert(std::is_same<typename SmallMatT::value_type::value_type, lSol_float_t>::value,
        "Lambda values ...");
    hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
        "The size ...");
    for (unsigned int i = 0; i < hyEdge_dimT; ++i)
      hy_assert(lambda_values[i].size() == 3 * n_shape_bdr_,
          "The size of ...");

    std::array<lSol_float_t, n_loc_dofs_> residual;
    residual.fill(0.);
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    const lSol_float_t eps = ldexp(1., -40);

    //calculate frequently used coefficients
    std::array<std::array<lSol_float_t, n_shape_fct_>, 2 * hyEdge_dim()> flux_ux, flux_uy, flux_q, flux_v, trace_v, trace_r;
    std::array<SmallVec<hyEdge_dim(), lSol_float_t>, 2 * hyEdge_dim()> loc_normal;
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      loc_normal[bdr] = hyper_edge.geometry.local_normal(bdr);
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        lSol_float_t h = integrate_bdr_phifunb(i, bdr, lambda_values[bdr].begin(), hyper_edge.geometry)
            - integrate_bdr_phifunv(i, bdr, ca.begin(), hyper_edge.geometry);
        flux_ux[bdr][i] = h * loc_normal[bdr][0];
        flux_uy[bdr][i] = h * loc_normal[bdr][1];
        if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] > 0)
        {
          flux_q[bdr][i] = integrate_bdr_phifunb(i, bdr, lambda_values[bdr].begin() + n_shape_bdr_, hyper_edge.geometry)
              - integrate_bdr_phifunv(i, bdr, ca.begin() + n_shape_fct_, hyper_edge.geometry);
          flux_v[bdr][i] = integrate_bdr_phifunb(i, bdr, lambda_values[bdr].begin() + 2 * n_shape_bdr_, hyper_edge.geometry)
              - integrate_bdr_phifunv(i, bdr, ca.begin() + 4 * n_shape_fct_, hyper_edge.geometry);
          trace_v[bdr][i] = 0;
          trace_r[bdr][i] = -(1. / delta_t_) * (integrate_bdr_phifunb(i, bdr, lambda_values[bdr].begin(), hyper_edge.geometry)
              - integrate_bdr_phifunb(i, bdr, hyper_edge.data.uh_old[bdr].begin(), hyper_edge.geometry));
        } else if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] < 0)
        {
          flux_q[bdr][i] = integrate_bdr_phifunb(i, bdr, lambda_values[bdr].begin() + n_shape_bdr_, hyper_edge.geometry)
              - integrate_bdr_phifunv(i, bdr, ca.begin() + n_shape_fct_, hyper_edge.geometry);
          flux_v[bdr][i] = 0;
          trace_v[bdr][i] = integrate_bdr_phifunv(i, bdr, ca.begin() + 4 * n_shape_fct_, hyper_edge.geometry)
              + tau_pvu_ * flux_ux[bdr][i];
          trace_r[bdr][i] = -(1. / delta_t_) * (integrate_bdr_phifunb(i, bdr, lambda_values[bdr].begin(), hyper_edge.geometry)
              - integrate_bdr_phifunb(i, bdr, hyper_edge.data.uh_old[bdr].begin(), hyper_edge.geometry));
        } else if (loc_normal[bdr][0] * loc_normal[bdr][0] < eps) {
          flux_q[bdr][i] = 0;
          flux_v[bdr][i] = 0;
          trace_v[bdr][i] = integrate_bdr_phifunv(i, bdr, ca.begin() + 4 * n_shape_fct_, hyper_edge.geometry)
              + tau_yvu_ * flux_uy[bdr][i];
          trace_r[bdr][i] = 0;
        } else {
          flux_q[bdr][i] = 0;
          flux_v[bdr][i] = 0;
          trace_v[bdr][i] = 0;
          trace_r[bdr][i] = 0;
        }
      }
    }

    //first eq.
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      residual[i] = integrate_vol_phifun(i, ca.begin() + n_shape_fct_, hyper_edge.geometry);
      residual[i] += integrate_vol_dphifun(i, 0, ca.begin(), hyper_edge.geometry);
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
        residual[i] -= integrate_bdr_phifunb(i, bdr, lambda_values[bdr].begin(), hyper_edge.geometry) * loc_normal[bdr][0];
    }
    
    //second eq
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      residual[n_shape_fct_ + i] = integrate_vol_phifun(i, ca.begin() + 2 * n_shape_fct_, hyper_edge.geometry);
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        for (unsigned int k = 0; k < n_shape_fct_; ++k)
        { 
          residual[n_shape_fct_ + i] += integrator::template integrate_vol_phiphiDphi<decltype(hyEdgeT::geometry)>(
              j, k, i, 0, hyper_edge.geometry) * ca[j] * ca[n_shape_fct_ + k];
        }
      }
    }
    
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        { 
          for (unsigned int k = 0; k < n_shape_fct_; ++k)
          {
            residual[n_shape_fct_ + i] -= 0.5 * integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, j, k, bdr, hyper_edge.geometry) * ca[j] * ca[n_shape_fct_ + k] * loc_normal[bdr][0];
          }
          for (unsigned int k = 0; k < n_shape_bdr_; ++k)
          {
            residual[n_shape_fct_ + i] -= 0.5 * integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, j, k, bdr, hyper_edge.geometry) * ca[j] * lambda_values[bdr][n_shape_bdr_ + k] * loc_normal[bdr][0];
          }
        }
        residual[n_shape_fct_ + i] -= tau_uqq_ * flux_q[bdr][i];
      }
    }

    //third eq
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      residual[2 * n_shape_fct_ + i] = integrate_vol_phifun(i, ca.begin() + 3 * n_shape_fct_, hyper_edge.geometry);
      residual[2 * n_shape_fct_ + i] += integrate_vol_dphifun(i, 1, ca.begin(), hyper_edge.geometry);
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        residual[2 * n_shape_fct_ + i] -= integrate_bdr_phifunb(i, bdr, lambda_values[bdr].begin(), hyper_edge.geometry)
            * loc_normal[bdr][1];
      }
    }
    
    //fourth equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      residual[3 * n_shape_fct_ + i] = integrate_vol_phifun(i, ca.begin() + 3 * n_shape_fct_, hyper_edge.geometry);
      residual[3 * n_shape_fct_ + i] += integrate_vol_dphifun(i, 0, ca.begin() + 4 * n_shape_fct_, hyper_edge.geometry);
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        residual[3 * n_shape_fct_ + i] -= trace_v[bdr][i] * loc_normal[bdr][0];
      }
    }
    
    //fifth equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      residual[4 * n_shape_fct_ + i] = integrate_vol_phifun(i, ca.begin() + 5 * n_shape_fct_, hyper_edge.geometry);
      residual[4 * n_shape_fct_ + i] += integrate_vol_dphifun(i, 0, ca.begin() + 6 * n_shape_fct_, hyper_edge.geometry);
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        residual[4 * n_shape_fct_ + i] -= trace_r[bdr][i] * loc_normal[bdr][0];
      }
    }

    //sixth equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      residual[5 * n_shape_fct_ + i] = integrate_vol_phifun(i, ca.begin() + 6 * n_shape_fct_, hyper_edge.geometry);
      residual[5 * n_shape_fct_ + i] += integrate_vol_dphifun(i, 0, ca.begin() + 5 * n_shape_fct_, hyper_edge.geometry);
      lSol_float_t f_int, q2_int;
      f_int = 2 * parameters::kappa * integrate_vol_dphifun(i, 0, ca.begin(), hyper_edge.geometry);
      q2_int = 0;
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        for (unsigned int k = 0; k < n_shape_fct_; ++k)
        {
          lSol_float_t h = integrator::template integrate_vol_phiphiDphi<decltype(hyEdgeT::geometry)>(
              j, k, i, 0, hyper_edge.geometry);
          f_int += 1.5 * h * ca[j] * ca[k];
          q2_int += h * ca[n_shape_fct_ + j] * ca[n_shape_fct_ + k];
        }
      }
      residual[5 * n_shape_fct_ + i] += f_int - integrate_vol_dphifun(i, 0, ca.begin() + 2 * n_shape_fct_, hyper_edge.geometry) 
          + 0.5 * q2_int;
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        //calculate f^ and (q^)^2
        lSol_float_t f_intb, q2_intb;
        f_intb = 2 * parameters::kappa * integrate_bdr_phifunv(i, bdr, ca.begin(), hyper_edge.geometry);
        q2_intb = 0;
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          for (unsigned int k = 0; k < n_shape_fct_; ++k)
          {
            f_intb += 1.5 * integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
                j, k, i, bdr, hyper_edge.geometry) * ca[j] * ca[k];
          }
        }
        for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        {
          for (unsigned int k = 0; k < n_shape_bdr_; ++k)
          {
            q2_intb += integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
                i, j, k, bdr, hyper_edge.geometry) * lambda_values[bdr][n_shape_bdr_ + j] 
                * lambda_values[bdr][n_shape_bdr_ + k];
          }
        }
        std::array<lSol_float_t, n_shape_fct_> u;
        for(unsigned int j = 0; j < n_shape_fct_; ++j)
          u[j] = ca[j];
        std::array<lSol_float_t, n_shape_bdr_> uh;
        for(unsigned int j = 0; j < n_shape_bdr_; ++j)
          uh[j] = lambda_values[bdr][j];
        f_intb -= integrate_bdr_phicompfun<decltype(hyEdgeT::geometry), parameters::tau_f>(
            i, uh, u, bdr, hyper_edge.geometry) * loc_normal[bdr][0];
        
        residual[5 * n_shape_fct_ + i] -= (f_intb + 0.5 * q2_intb) * loc_normal[bdr][0];
        //p^ and z^
        if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] < 0)  //left
        {
          residual[5 * n_shape_fct_ + i] -= (integrate_bdr_phifunv(i, bdr, ca.begin() + 5 * n_shape_fct_, hyper_edge.geometry)
              + tau_pzu_ * flux_ux[bdr][i]) * loc_normal[bdr][0];
          residual[5 * n_shape_fct_ + i] -= -(integrate_bdr_phifunv(i, bdr, ca.begin() + 2 * n_shape_fct_, hyper_edge.geometry)
              + tau_ppu_ * flux_ux[bdr][i]) * loc_normal[bdr][0];
        } else if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] > 0) //right
        {
          residual[5 * n_shape_fct_ + i] -= (integrate_bdr_phifunv(i, bdr, ca.begin() + 5 * n_shape_fct_, hyper_edge.geometry)
              + tau_mzu_ * flux_ux[bdr][i] + tau_mzv_ * flux_v[bdr][i]) * loc_normal[bdr][0];
          residual[5 * n_shape_fct_ + i] -= -(integrate_bdr_phifunv(i, bdr, ca.begin() + 2 * n_shape_fct_, hyper_edge.geometry)
              + tau_mpu_ * flux_ux[bdr][i] + tau_mpv_ * flux_v[bdr][i]) * loc_normal[bdr][0];
        }
      }
      residual[5 * n_shape_fct_ + i] += integrate_vol_dphifun(i, 1, ca.begin() + 4 * n_shape_fct_, hyper_edge.geometry);
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        residual[5 * n_shape_fct_ + i] -= trace_v[bdr][i] * loc_normal[bdr][1];
      }
    }
      
    //seventh equation
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      residual[6 * n_shape_fct_ + i] = (1. / delta_t_) * (integrate_vol_phifun(i, ca.begin(), hyper_edge.geometry)
          - integrate_vol_phifun(i, hyper_edge.data.u_old.begin(), hyper_edge.geometry));
      residual[6 * n_shape_fct_ + i] += integrate_vol_phifun(i, ca.begin() + 6 * n_shape_fct_, hyper_edge.geometry);
    }

    return residual;
  }

  /*!***********************************************************************************************
   * \brief   Solve local problem (with right-hand side from skeletal).
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatT     The data type of the \c lambda_values.
   * \param   lambda_values Encodes uh, qh, vh.
   * \param   coeff         Coefficients of u, q, p, s, v, z, r to be tested.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Point of time the problem is solved.
   * \retval  residual      Residual (should be zero).
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& residual_flux(const SmallMatInT& lambda_values_in_uc,
                                                           SmallMatOutT& lambda_values_out,
                                                           hyEdgeT& hyper_edge,
                                                           const lSol_float_t time) const
  {
    //ensure dirichlet conditions are met
    SmallMatInT lambda_values_in = lambda_values_in_uc;
    make_skeleton(lambda_values_in, hyper_edge, time);
    //calculate integral coefficients
    //set_coeff(hyper_edge);

    std::array<SmallVec<hyEdge_dim(), lSol_float_t>, 2 * hyEdge_dim()> loc_normal;
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr) 
      loc_normal[bdr] = hyper_edge.geometry.local_normal(bdr);
    SmallVec<n_loc_dofs_, lSol_float_t> coeff(0.);
    newton(lambda_values_in, coeff, hyper_edge, time);
    const lSol_float_t eps = ldexp(1., -30);
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      using parameters = parametersT<hyEdge_dim(), lSol_float_t>;
      if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
	      {
      	  lambda_values_out[bdr][i] = 0;
  	      lambda_values_out[bdr][n_shape_bdr_ + i] = 0;
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] = 0;
  	      lSol_float_t uh_int = 0, u_int = 0, qh_int = 0, q_int = 0, vh_int = 0, v_int = 0;
    	    for (unsigned int j = 0; j < n_shape_bdr_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
	          uh_int += c * lambda_values_in[bdr][j];
      	    qh_int += c * lambda_values_in[bdr][n_shape_bdr_ + j];
      	    vh_int += c * lambda_values_in[bdr][2 * n_shape_bdr_ + j];
      	  }
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
      	    u_int += c * coeff[j];
	          q_int += c * coeff[n_shape_fct_ + j];
      	    v_int += c * coeff[4 * n_shape_fct_ + j];
      	  }
      	  //z_hat - p_hat
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j) 
	        {
            lSol_float_t c = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
      	    lambda_values_out[bdr][i] += c * (coeff[5 * n_shape_fct_ + j] - coeff[2 * n_shape_fct_ + j]);
      	  }
      	  lambda_values_out[bdr][i] += tau_mzu_ * (uh_int - u_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] -= tau_mpu_ * (uh_int - u_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] += tau_mzv_ * (vh_int - v_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] -= tau_mpv_ * (vh_int - v_int) * loc_normal[bdr][0];
      	  //f_hat
      	  lambda_values_out[bdr][i] += 2 * parameters::kappa * u_int;
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
      	  {
            for (unsigned int k = 0; k < n_shape_fct_; ++k) 
	          {
              lSol_float_t c = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
      	      lambda_values_out[bdr][i] += 1.5 * c * coeff[j] * coeff[k];
	          }
	        }
      	  std::array<lSol_float_t, n_shape_bdr_> uh_arr;
      	  std::array<lSol_float_t, n_shape_fct_> u_arr;
      	  for (unsigned int j = 0; j < n_shape_bdr_; ++j)
            uh_arr[j] = lambda_values_in[bdr][j];
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
            u_arr[j] = lambda_values_in[bdr][j];
          lambda_values_out[bdr][i] -= integrate_bdr_psicompfun<decltype(hyEdgeT::geometry), parameters::tau_f>(
		    	  i, uh_arr, u_arr, bdr, hyper_edge.geometry) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] *= loc_normal[bdr][0];
    	    //uq^
	        //0.5 u q
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
	        {
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
            {
              lSol_float_t c = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
	            lambda_values_out[bdr][n_shape_bdr_ + i] += 0.5 * c * coeff[j] * coeff[n_shape_fct_ + k];
	          }
      	  }
      	  //0.5 u q^
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
      	  {
            for (unsigned int k = 0; k < n_shape_bdr_; ++k)
            {
              lSol_float_t c = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
	            lambda_values_out[bdr][n_shape_bdr_ + i] += 0.5 * c * coeff[j] * lambda_values_in[bdr][n_shape_bdr_ + k];
      	    }
      	  }
      	  lambda_values_out[bdr][n_shape_bdr_ + i] += tau_uqq_ * (qh_int - q_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][n_shape_bdr_ + i] *= loc_normal[bdr][0];
      	  //v^
      	  //TODO allow for v_R
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] += v_int;
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] += tau_pvu_ * (uh_int - u_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] *= loc_normal[bdr][0];
      	}
      }
      if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] > 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
	      {
      	  lambda_values_out[bdr][i] = 0;
  	      lambda_values_out[bdr][n_shape_bdr_ + i] = 0;
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] = 0;
  	      lSol_float_t uh_int = 0, u_int = 0, qh_int = 0, q_int = 0, vh_int = 0, v_int = 0;
    	    for (unsigned int j = 0; j < n_shape_bdr_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
	          uh_int += c * lambda_values_in[bdr][j];
      	    qh_int += c * lambda_values_in[bdr][n_shape_bdr_ + j];
      	    vh_int += c * lambda_values_in[bdr][2 * n_shape_bdr_ + j];
      	  }
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
      	    u_int += c * coeff[j];
	          q_int += c * coeff[n_shape_fct_ + j];
      	    v_int += c * coeff[4 * n_shape_fct_ + j];
      	  }
      	  //z_hat - p_hat
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j) 
	        {
            lSol_float_t c = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
      	    lambda_values_out[bdr][i] += c * (coeff[5 * n_shape_fct_ + j] - coeff[2 * n_shape_fct_ + j]);
      	  }
      	  lambda_values_out[bdr][i] += tau_mzu_ * (uh_int - u_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] -= tau_mpu_ * (uh_int - u_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] += tau_mzv_ * (vh_int - v_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] -= tau_mpv_ * (vh_int - v_int) * loc_normal[bdr][0];
      	  //f_hat
      	  lambda_values_out[bdr][i] += 2 * parameters::kappa * u_int;
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
      	  {
            for (unsigned int k = 0; k < n_shape_fct_; ++k) 
	          {
              lSol_float_t c = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
      	      lambda_values_out[bdr][i] += 1.5 * c * coeff[j] * coeff[k];
	          }
	        }
      	  std::array<lSol_float_t, n_shape_bdr_> uh_arr;
      	  std::array<lSol_float_t, n_shape_fct_> u_arr;
      	  for (unsigned int j = 0; j < n_shape_bdr_; ++j)
            uh_arr[j] = lambda_values_in[bdr][j];
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
            u_arr[j] = lambda_values_in[bdr][j];
          lambda_values_out[bdr][i] -= integrate_bdr_psicompfun<decltype(hyEdgeT::geometry), parameters::tau_f>(
		    	  i, uh_arr, u_arr, bdr, hyper_edge.geometry) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] *= loc_normal[bdr][0];
    	    //uq^
	        //0.5 u q
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
	        {
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
            {
              lSol_float_t c = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
	            lambda_values_out[bdr][n_shape_bdr_ + i] += 0.5 * c * coeff[j] * coeff[n_shape_fct_ + k];
	          }
      	  }
      	  //0.5 u q^
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
      	  {
            for (unsigned int k = 0; k < n_shape_bdr_; ++k)
            {
              lSol_float_t c = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
	            lambda_values_out[bdr][n_shape_bdr_ + i] += 0.5 * c * coeff[j] * lambda_values_in[bdr][n_shape_bdr_ + k];
      	    }
      	  }
      	  lambda_values_out[bdr][n_shape_bdr_ + i] += tau_uqq_ * (qh_int - q_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][n_shape_bdr_ + i] *= loc_normal[bdr][0];
      	  //v^
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] = lambda_values_in[bdr][2 * n_shape_bdr_ + i];
      	}
      }
      if (loc_normal[bdr][0] * loc_normal[bdr][0] < eps && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
      	{
  	      lambda_values_out[bdr][i] = 0;
      	  lSol_float_t uh_int = 0, u_int = 0, v_int = 0;
      	  for (unsigned int j = 0; j < n_shape_bdr_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
	          uh_int += c * lambda_values_in[bdr][j];
      	  }
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
      	    u_int += c * coeff[j];
	          v_int += c * coeff[4 * n_shape_fct_ + j];
      	  }
      	  lambda_values_out[bdr][i] += v_int;
      	  lambda_values_out[bdr][i] += tau_yvu_ * (uh_int - u_int) * loc_normal[bdr][1];
      	  lambda_values_out[bdr][i] *= loc_normal[bdr][1];
      	  lambda_values_out[bdr][n_shape_bdr_ + i] = lambda_values_in[bdr][n_shape_bdr_ + i];
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] = lambda_values_in[bdr][2 * n_shape_bdr_ + i];
      	}
      }
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
  	      lambda_values_out[bdr][i] = 0;
  	      lambda_values_out[bdr][n_shape_bdr_ + i] = 0;
  	      lambda_values_out[bdr][2 * n_shape_bdr_ + i] = 0;
        }
      }
    }
    return lambda_values_out;
  }

  template <typename hyEdgeT, typename SmallMatT>
  lSol_float_t newton(const SmallMatT& lambda_values, SmallVec<n_loc_dofs_, lSol_float_t>& coeff,
                      hyEdgeT& hyper_edge, const lSol_float_t time) const
  {
    const lSol_float_t eps = ldexp(1., -40);
    SmallVec<n_loc_dofs_, lSol_float_t> res = get_residual(lambda_values, coeff, hyper_edge, time);
    lSol_float_t ra = norm_2(res);
    for (unsigned int i = 0; ra > eps && i < 100; ++i)
    {
      lSol_float_t rn;
      lSol_float_t stepsize = 1.;
      SmallVec<n_loc_dofs_, lSol_float_t> cn = coeff;
      SmallVec<n_loc_dofs_, lSol_float_t> step = res / jacobi(lambda_values, coeff, hyper_edge, time);
      ra = norm_2(res);
      do
      {
        //std::cout << "aktuelles residuum: " << ra << "\n";
        cn = coeff - stepsize * (step);
        res = get_residual(lambda_values, cn, hyper_edge, time);
        rn = norm_2(res);
        //std::cout << "neues residuum: " << rn << "\n";
        stepsize *= .5;
      } while (ra < rn);
      coeff = cn;
      ra = rn;
    }
    return ra;
  }

  /*********************************************************************************
   * Fills projection of intitial u into data.u_old
   * Computes local contribution to data.uh_old from initial and dirichlet
  *********************************************************************************/
  template <typename hyEdgeT, typename SmallMatT>
  SmallMatT& make_initial(SmallMatT& lambda_values,
      hyEdgeT& hyper_edge, const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<hyEdge_dim(), lSol_float_t>;
    //project initial to u
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      hyper_edge.data.u_old[i] = integrator::template integrate_volUni_phifunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::initial, Point<hyEdge_dimT, lSol_float_t> > (i, hyper_edge.geometry, 0.);
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      for (unsigned int i = 0; i < n_shape_bdr_; ++i)
      {
        if (is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
        {
          lambda_values[bdr][i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::dirichlet_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);

        } else 
        {
          lambda_values[bdr][i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::initial, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
        }
        hyper_edge.data.uh_old[bdr][i] = lambda_values[bdr][i];
        lambda_values[bdr][n_shape_bdr_ + i] = 0.;
        lambda_values[bdr][2 * n_shape_bdr_ + i] = 0.;
      }
    }
    set_skeleton_data(lambda_values, hyper_edge);
    return lambda_values;
  }

  template <typename hyEdgeT, typename SmallMatT>
  static inline void set_skeleton_data(const SmallMatT& lambda_values,
      hyEdgeT& hyper_edge)
  {
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        hyper_edge.data.uh_old[bdr][i] = lambda_values[bdr][i];
    }
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline void set_bulk_data(const SmallMatT& lambda_values, 
      hyEdgeT& hyper_edge, const lSol_float_t time) const
  {
    SmallVec<n_loc_dofs_, lSol_float_t> coeff;
    newton(lambda_values, coeff, hyper_edge, time);
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      hyper_edge.data.u_old[i] = coeff[i];
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline void set_data(const SmallMatT& lambda_values, 
      hyEdgeT& hyper_edge, const lSol_float_t time) const
  {
    set_bulk_data(lambda_values, hyper_edge, time);
    set_skeleton_data(lambda_values, hyper_edge);
  }

  /*********************************************************************************
   * Updates uh at the boundary to new dirichlet. Leaves other skeleton values as is.
  *********************************************************************************/
  template <typename hyEdgeT, typename SmallMatOutT>
  static inline void make_skeleton(SmallMatOutT& lambda_values_out,
      hyEdgeT& hyper_edge, const lSol_float_t time = 0.)
  {
    using parameters = parametersT<hyEdge_dim(), lSol_float_t>;
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
          lambda_values_out[bdr][i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::dirichlet_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
      }
    }
  }
  /*!***********************************************************************************************
   * \brief   Evaluate squared local L2 error.
   *
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hy_edge           The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which error is evaluated.
   * \retval  err               Local squared L2 error.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<lSol_float_t, 1U> errors(const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>,
                                                       2 * hyEdge_dimT>& lambda_values,
                                      hyEdgeT& hy_edge,
                                      const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    return std::array<lSol_float_t, 1U>({integrator::template integrate_vol_diffsquare_discana<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result, Point<hyEdge_dimT, lSol_float_t> >(hy_edge.data.u_old.data(),
                                                                      hy_edge.geometry, time)});
  }

private:
  template <typename geom_t, lSol_float_t fun(const lSol_float_t)>
  static lSol_float_t integrate_bdr_phicompfun(const unsigned int i,
                                        const std::array<lSol_float_t, n_shape_bdr_>& coeff1,
                                        const std::array<lSol_float_t, n_shape_fct_>& coeff2,
                                        const unsigned int bdr,
                                        geom_t& geom)
  {
    //recover private members
    lSol_float_t result = 0;
    //const unsigned int dim = 2;
    const auto &qp = integrator::quad_points;
    const unsigned int n_points = qp.size();
    const std::array<lSol_float_t, n_points> &qw = integrator::quad_weights;
    const auto &shape_fcts = integrator::shape_fcts_at_quad;
    const unsigned int n_fct = shape_fcts.size();
    const auto &shape_bdr = integrator::shape_fcts_at_bdr;

    //determine boundary
    const unsigned int bdr_d = bdr/2, bdr_i = bdr % 2;

    std::array<std::array<unsigned int, hyEdge_dim()>, n_shape_fct_> phi_ind;
    for(unsigned int j = 0; j < n_shape_fct_; ++j)
      phi_ind[j] = TPP::Hypercube<hyEdge_dim()>::index_decompose(j, n_fct);
    //for every point
    //specialized for 1d boundary of 2d square
    for(unsigned int p = 0; p < n_points; ++p)
    {
      //evaluate inner difference
     lSol_float_t ws = 0; 
     for(unsigned int j = 0; j < n_shape_fct_; ++j)
     {
       std::array<unsigned int, 2> ind_j = phi_ind[j];
       lSol_float_t wsj = 1;
       for(unsigned int k = 0; k < 2; k++)
       {
         if (k == bdr_d)
           wsj *= shape_bdr[ind_j[k]][bdr_i];
         else
           wsj *= shape_fcts[ind_j[k]][p];
       }
       ws -= coeff2[j] * wsj;
     }
     for (unsigned int j = 0; j < n_shape_bdr_; ++j)
     {
       ws += coeff1[j] * shape_fcts[j][p];
     }
     //evaluate phi_i
     lSol_float_t wpi = 1.;
     for (unsigned int k = 0; k < 2; ++k)
     {
       if (k == bdr_d)
         wpi *= shape_bdr[phi_ind[i][k]][bdr_i];
       else
         wpi *= shape_fcts[phi_ind[i][k]][p];
     }

     result += qw[p] * fun(ws) * wpi;
    }
    return result * geom.face_area(bdr);
  }


  template <typename geom_t, lSol_float_t fun(const lSol_float_t)>
  static lSol_float_t integrate_bdr_psicompfun(const unsigned int i,
                                        const std::array<lSol_float_t, n_shape_bdr_>& coeff1,
                                        const std::array<lSol_float_t, n_shape_fct_>& coeff2,
                                        const unsigned int bdr,
                                        geom_t& geom)
  {
    //recover private members
    lSol_float_t result = 0;
    //const unsigned int dim = 2;
    const auto &qp = integrator::quad_points;
    const unsigned int n_points = qp.size();
    const std::array<lSol_float_t, n_points> &qw = integrator::quad_weights;
    const auto &shape_fcts = integrator::shape_fcts_at_quad;
    const unsigned int n_fct = shape_fcts.size();
    const auto &shape_bdr = integrator::shape_fcts_at_bdr;

    //determine boundary
    const unsigned int bdr_d = bdr/2, bdr_i = bdr % 2;

    std::array<std::array<unsigned int, hyEdge_dim()>, n_shape_fct_> phi_ind;
    for(unsigned int j = 0; j < n_shape_fct_; ++j)
      phi_ind[j] = TPP::Hypercube<hyEdge_dim()>::index_decompose(j, n_fct);
    //works for hyEdge_dim >= 2
    std::array<unsigned int, hyEdge_dim()-1> psi_ind = TPP::Hypercube<hyEdge_dim() - 1>::index_decompose(i, n_fct);
    //for every point
    //specialized for 1d boundary of 2d square
    for(unsigned int p = 0; p < n_points; ++p)
    {
      //evaluate inner difference
     lSol_float_t ws = 0; 
     for(unsigned int j = 0; j < n_shape_fct_; ++j)
     {
       std::array<unsigned int, 2> ind_j = phi_ind[j];
       lSol_float_t wsj = 1;
       for(unsigned int k = 0; k < 2; k++)
       {
         if (k == bdr_d)
           wsj *= shape_bdr[ind_j[k]][bdr_i];
         else
           wsj *= shape_fcts[ind_j[k]][p];
       }
       ws -= coeff2[j] * wsj;
     }
     for (unsigned int j = 0; j < n_shape_bdr_; ++j)
     {
       ws += coeff1[j] * shape_fcts[j][p];
     }
     //evaluate psi_i
     lSol_float_t wpi = 1.;
     for (unsigned int k = 0; k < hyEdge_dim() - 1; ++k)
       wpi *= shape_fcts[psi_ind[k]][p];
     result += qw[p] * fun(ws) * wpi;
    }
    return result * geom.face_area(bdr);
  }


  template <typename geom_t, lSol_float_t fun(const lSol_float_t)>
  static lSol_float_t integrate_bdr_phiphicompfun(const unsigned int i,
                                        const unsigned int j,
                                        const std::array<lSol_float_t, n_shape_bdr_>& coeff1,
                                        const std::array<lSol_float_t, n_shape_fct_>& coeff2,
                                        const unsigned int bdr,
                                        geom_t& geom)
  {
    //recover private members
    lSol_float_t result = 0;
    //const unsigned int dim = 2;
    const auto &qp = integrator::quad_points;
    const unsigned int n_points = qp.size();
    const std::array<lSol_float_t, n_points> &qw = integrator::quad_weights;
    const auto &shape_fcts = integrator::shape_fcts_at_quad;
    const unsigned int n_fct = shape_fcts.size();
    const auto &shape_bdr = integrator::shape_fcts_at_bdr;

    //determine boundary
    const unsigned int bdr_d = bdr/2, bdr_i = bdr % 2;

    std::array<std::array<unsigned int, hyEdge_dim()>, n_shape_fct_> phi_ind;
    for(unsigned int j = 0; j < n_shape_fct_; ++j)
      phi_ind[j] = TPP::Hypercube<hyEdge_dim()>::index_decompose(j, n_fct);
    //for every point
    //specialized for 1d boundary of 2d square
    for(unsigned int p = 0; p < n_points; ++p)
    {
      //evaluate inner difference
     lSol_float_t ws = 0; 
     for(unsigned int j = 0; j < n_shape_fct_; ++j)
     {
       std::array<unsigned int, 2> ind_j = phi_ind[j];
       lSol_float_t wsj = 1;
       for(unsigned int k = 0; k < 2; k++)
       {
         if (k == bdr_d)
           wsj *= shape_bdr[ind_j[k]][bdr_i];
         else
           wsj *= shape_fcts[ind_j[k]][p];
       }
       ws -= coeff2[j] * wsj;
     }
     for (unsigned int j = 0; j < n_shape_bdr_; ++j)
     {
       ws += coeff1[j] * shape_fcts[j][p];
     }
     //evaluate phi_i and phi_j
     lSol_float_t wpi = 1., wpj = 1.;
     for (unsigned int k = 0; k < 2; ++k)
     {
       if (k == bdr_d)
       {
         wpi *= shape_bdr[phi_ind[i][k]][bdr_i];
         wpj *= shape_bdr[phi_ind[j][k]][bdr_i];
       }
       else
       {
         wpi *= shape_fcts[phi_ind[i][k]][p];
         wpj *= shape_fcts[phi_ind[j][k]][p];
       }
     }

     result += qw[p] * fun(ws) * wpi * wpj;
    }
    return result * geom.face_area(bdr);
  }

  template <typename geom_t, typename SmallVecItT>
  static lSol_float_t integrate_vol_phifun(const unsigned int i,
                                            SmallVecItT funcv_it,
                                            geom_t& geom)
  {
    lSol_float_t r = 0.;
    for (unsigned int j = 0; j < n_shape_fct_; ++j, ++funcv_it)
      r += *funcv_it * integrator::template integrate_vol_phiphi<geom_t>(
          i, j, geom);
    return r;
  }

  template <typename geom_t, typename SmallVecItT>
  static lSol_float_t integrate_vol_dphifun(const unsigned int i,
                                            const unsigned int dim_der,
                                            SmallVecItT funcv_it,
                                            geom_t& geom)
  {
    lSol_float_t r = 0.;
    for (unsigned int j = 0; j < n_shape_fct_; ++j, ++funcv_it)
      r += *funcv_it * integrator::template integrate_vol_phiDphi<geom_t>(
          j, i, dim_der, geom);
    return r;
  }
  
  template <typename geom_t, typename SmallVecItT>
  static lSol_float_t integrate_bdr_phifunb(const unsigned int i,
                                            const unsigned int bdr,
                                            SmallVecItT funcb_it,
                                            geom_t& geom)
  {
    lSol_float_t r = 0.;
    for (unsigned int j = 0; j < n_shape_bdr_; ++j, ++funcb_it)
      r += *funcb_it * integrator::template integrate_bdr_phipsi<geom_t>(
            i, j, bdr, geom);
    return r;
  }
  
  template <typename geom_t, typename SmallVecItT>
  static lSol_float_t integrate_bdr_phifunv(const unsigned int i,
                                            const unsigned int bdr,
                                            SmallVecItT funcv_it,
                                            geom_t& geom)
  {
    lSol_float_t r = 0.;
    for (unsigned int j = 0; j < n_shape_fct_; ++j, ++funcv_it)
      r += *funcv_it * integrator::template integrate_bdr_phiphi<geom_t>(
            i, j, bdr, geom);
    return r;
  }

};

}
