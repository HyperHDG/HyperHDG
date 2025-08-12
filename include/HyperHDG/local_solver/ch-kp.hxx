#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

#include <algorithm>
#include <tuple>
#include <iostream>

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
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1.;
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
    return 2 * Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
 private:
  // -----------------------------------------------------------------------------------------------
  // Private, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Number of local shape functions (with respect to all spatial dimensions).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_fct_ = n_glob_dofs_per_node() * (poly_deg + 1) / 2;
  /*!***********************************************************************************************
   * \brief   Number of local  shape functions (with respect to a face / hypernode).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node() / 2;
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
  const lSol_float_t tau_f_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_uqq_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_yvu_;
//TODO: Add support for non-constant tau_f. This requires new integration routines.
  /*!***********************************************************************************************
   * \brief   Parameter theta that defines the one-step theta scheme.
   ************************************************************************************************/
  const lSol_float_t theta_ = 1.;
  /*!***********************************************************************************************
   * \brief   Time step size.
   ************************************************************************************************/
  const lSol_float_t delta_t_;
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
    SmallVec<n_shape_fct_, lSol_float_t> u_old;
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
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef std::vector<lSol_float_t> constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   constru       Constructor object.
   ************************************************************************************************/
  Chkp(const constructor_value_type& constru = std::vector(11, 1.))
  : delta_t_(constru[0]), tau_ppu_(constru[1]), tau_mpu_(constru[2]), tau_mpv_(constru[3]), 
  tau_pzu_(constru[4]), tau_mzu_(constru[5]), tau_mzv_(constru[6]), tau_pvu_(constru[7]), 
  tau_f_(constru[8]), tau_uqq_(constru[9]), tau_yvu_(constru[10])
  {
  }
  /*!***********************************************************************************************
   * \brief   Solve local problem (with right-hand side from skeletal).
   *
   * \tparam  hyEdgeT       The geometry type / typename of the considered hyEdge's geometry.
   * \tparam  SmallMatT     The data type of the \c lambda_values.
   * \param   lambda_values Encodes uh, qh.
   * \param   coeff         Coefficients of u, q, p, s, v, z, r to be tested.
   * \param   hyper_edge    The geometry of the considered hyperedge (of typename GeomT).
   * \param   time          Point of time the problem is solved.
   * \retval  residual      Residual (should be zero).
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatT>
  inline std::array<lSol_float_t, n_loc_dofs_> get_residual(const SmallMatT& lambda_values,
                                                            const std::array<SmallVec<n_shape_fct_, lSol_float_t>, 7> coeff,
                                                            hyEdgeT& hyper_edge,
                                                            const lSol_float_t time) const
  {
    static_assert(std::is_same<typename SmallMatT::value_type::value_type, lSol_float_t>::value,
        "Lambda values ...");
    hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
        "The size ...");
    for (unsigned int i = 0; i < hyEdge_dimT; ++i)
      hy_assert(lambda_values[i].size() == 2 * n_shape_bdr_,
          "The size of ...");

    std::array<lSol_float_t, n_loc_dofs_> residual;
    residual.fill(0.);
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    //rearrange lambda values
    std::array<SmallVec<n_shape_bdr_, lSol_float_t>, 2 * hyEdge_dim()> u_hat, q_hat;
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr) 
    {
      for (unsigned int i = 0; i < n_shape_bdr_; i++)
      {
        u_hat[bdr][i] = lambda_values[bdr][i];
        q_hat[bdr][i] = lambda_values[bdr][n_shape_bdr_ + i];
      }
    }

    const lSol_float_t eps = 1./ 1073741824.;

    //calculate frequently used coefficients
    SmallSquareMat<n_shape_fct_, lSol_float_t> mass, mdx, mdy;
    std::array<SmallMat<n_shape_fct_, n_shape_bdr_, lSol_float_t>, 2 * hyEdge_dim()> bdr_sh_sk;
    std::array<SmallSquareMat<n_shape_fct_, lSol_float_t>, 2 * hyEdge_dim()> bdr_sh_sh;
    std::array<SmallVec<hyEdge_dim(), lSol_float_t>, 2 * hyEdge_dim()> loc_normal;
    std::array<SmallSquareMat<n_shape_fct_, lSol_float_t>, n_shape_fct_> tv_shx_sh_sh;
    std::array<std::array<SmallSquareMat<n_shape_fct_, lSol_float_t>, 2 * hyEdge_dim()>, n_shape_fct_> tb_sh_sh_sh;
    std::array<std::array<SmallMat<n_shape_fct_, n_shape_bdr_, lSol_float_t>, 2 * hyEdge_dim()>, n_shape_fct_> tb_sh_sh_sk;
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
    std::cout << bdr_sh_sk[0] << std::endl;
    std::cout << u_hat[0] << std::endl;
    SmallVec<n_shape_fct_, lSol_float_t> helper;

    //first eq.
    helper = mass * coeff[1] + mdx * coeff[0];
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr) 
      helper -= (bdr_sh_sk[bdr] * u_hat[bdr]) * loc_normal[bdr][0];
    std::cout << helper << std::endl;
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      residual[i] += helper[i];
    }

    //second eq
    helper = mass * coeff[2];
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr) 
    {
      helper -= tau_uqq_ * (bdr_sh_sk[bdr] * q_hat[bdr] - bdr_sh_sh[bdr] * coeff[1]) * loc_normal[bdr][0] * loc_normal[bdr][0];
    }
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
      residual[n_shape_fct_ + i] = helper[i];
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      helper = tv_shx_sh_sh[i] * coeff[1];
      residual[n_shape_fct_ + i] += scalar_product(coeff[0], helper);
      helper = SmallVec<n_shape_fct_, lSol_float_t>(0.);
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr) 
      {
        helper += 0.5 * (tb_sh_sh_sk[i][bdr] * q_hat[bdr] + tb_sh_sh_sh[i][bdr] * coeff[1]) * loc_normal[bdr][0];
      }
      residual[n_shape_fct_ + i] -= scalar_product(coeff[0], helper);
    }

    //third eq
    helper = mass * coeff[3] + mdy * coeff[0];
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      helper -= bdr_sh_sk[bdr] * u_hat[bdr] * loc_normal[bdr][1];
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
      residual[2 * n_shape_fct_ + i] = helper[i];

    //fourth equation
    helper = mass * coeff[3] + mdx * coeff[4];
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      helper -= bdr_sh_sh[bdr] * coeff[4] * loc_normal[bdr][0];
      if(loc_normal[bdr][1] * loc_normal[bdr][1] > eps && loc_normal[bdr][0] * loc_normal[bdr][0] < eps) //on H
        helper -= tau_yvu_ * (bdr_sh_sk[bdr] * u_hat[bdr] - bdr_sh_sh[bdr] * coeff[0]) 
          * loc_normal[bdr][1] * loc_normal[bdr][0];
      else if(loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] < 0) //on V_left
        helper -= tau_pvu_ * (bdr_sh_sk[bdr] * u_hat[bdr] - bdr_sh_sh[bdr] * coeff[0]) 
          * loc_normal[bdr][0] * loc_normal[bdr][0];
    }
    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
      residual[3 * n_shape_fct_ + i] = helper[i];

    return residual;
  }

private:
  template <typename geom_t, lSol_float_t fun(const lSol_float_t)>
  static lSol_float_t integrate_bdr_phicompfun(const unsigned int i,
                                        const std::array<lSol_float_t, n_shape_fct_>& coeff1,
                                        const std::array<lSol_float_t, n_shape_bdr_>& coeff2,
                                        const unsigned int bdr,
                                        geom_t& geom)
  {
    //recover private members
    lSol_float_t result = 0;
    const unsigned int dim = 2;
    const auto &qp = integrator::quad_points;
    const unsigned int n_points = qp.size();
    const std::array<lSol_float_t, n_points> &qw = integrator::quad_weights;
    const auto &shape_fcts = integrator::shape_fcts_at_quad;
    const unsigned int n_fct = shape_fcts.size();
    const auto &shape_bdr = integrator::shape_fcts_at_bdr;

    //determine boundary
    const unsigned int bdr_d = bdr/2, bdr_i = bdr % 2;

    std::array<std::array<unsigned int, hyEdge_dim()>, n_shape_fct_> phi_ind;
    for(int j = 0; j < n_shape_fct_; ++j)
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
       ws += coeff1[j] * wsj;
     }
     for (unsigned int j = 0; j < n_shape_bdr_; ++j)
     {
       ws -= coeff2[j] * shape_fcts[j][p];
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
};
}
