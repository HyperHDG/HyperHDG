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
  static constexpr param_float_t tau_fr = 1.;

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
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return 1; }
  
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
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   *
   * This allows to the use of this quantity as template parameter in member functions.
   ************************************************************************************************/
  static constexpr unsigned int system_dim = system_dimension();
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   *
   * This allows to the use of this quantity as template parameter in member functions.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dim = node_system_dimension();
  /*!***********************************************************************************************
   * \brief   Find out whether a node is of Dirichlet type.
   ************************************************************************************************/
  template <typename parameters>
  static constexpr bool is_dirichlet(const unsigned int node_type)
  {
    return std::find(parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                     node_type) != parameters::dirichlet_nodes.end();
  }

  template <typename parameters>
  static constexpr bool is_neumann(const unsigned int node_type)
  {
    return std::find(parameters::neumann_nodes.begin(), parameters::neumann_nodes.end(),
                     node_type) != parameters::neumann_nodes.end();
  }

  template <typename parameters>
  static constexpr bool is_right(const unsigned int node_type)
  {
    return std::find(parameters::right_nodes.begin(), parameters::right_nodes.end(),
                     node_type) != parameters::right_nodes.end();
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
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_f_;
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

  const lSol_float_t onepointfive = 1.5; 
   
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
    SmallVec<n_shape_fct_, lSol_float_t> q_old = SmallVec<n_shape_fct_, lSol_float_t>(0.);
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
    typedef std::array<lSol_float_t, 2U> error_t;
    /*!*********************************************************************************************
     *  \brief  Define how initial error is generated.
     **********************************************************************************************/
    static error_t initial_error()
    {
      std::array<lSol_float_t, 2U> summed_error;
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
  Chkp(const constructor_value_type& constru = std::vector<double>({1., 3., 4., 1., 1., 1., 1., 1., -3., -2., 4.}))
  : delta_t_(constru[0]), tau_ppu_(constru[1]), tau_mpu_(constru[2]), tau_mpv_(constru[3]), 
  tau_pzu_(constru[4]), tau_mzu_(constru[5]), tau_mzv_(constru[6]), tau_pvu_(constru[7]), 
  tau_uqq_(constru[8]), tau_yvu_(constru[9]), tau_f_(constru[10])
  {
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallSquareMat<n_loc_dofs_, lSol_float_t> jacobi(const SmallMatT& lambda_values,
                                                          const SmallVec<n_loc_dofs_, lSol_float_t>& ca,
                                                          hyEdgeT& hyper_edge,
                                                          const lSol_float_t time) const
  {
    SmallSquareMat<n_loc_dofs_, lSol_float_t> ret(0.);
    SmallVec<hyEdge_dim(), lSol_float_t> grad;
    const lSol_float_t eps = ldexp(1., -40);
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
      //mass components
        lSol_float_t mass_ij = integrator::template integrate_vol_phiphi<decltype(hyEdgeT::geometry)>(
            i, j, hyper_edge.geometry);
        ret(i, n_shape_fct_ + j) = mass_ij;
        ret(n_shape_fct_ + i, 2 * n_shape_fct_ + j) = mass_ij;
        ret(2 * n_shape_fct_ + i, 3 * n_shape_fct_ + j) = mass_ij;
        ret(3 * n_shape_fct_ + i, 3 * n_shape_fct_ + j) = mass_ij;
        ret(4 * n_shape_fct_ + i, 5 * n_shape_fct_ + j) = mass_ij;
        ret(5 * n_shape_fct_ + i, 6 * n_shape_fct_ + j) = mass_ij;
        ret(6 * n_shape_fct_ + i, j) = (1. / delta_t_) * mass_ij;
        ret(6 * n_shape_fct_ + i, 6 * n_shape_fct_ + j) = mass_ij;

        //bulk derivatives
        grad = integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
            i, j, hyper_edge.geometry);
        lSol_float_t mdx_ij = grad[0];
        lSol_float_t mdy_ij = grad[1];
        ret(i, j) = mdx_ij;
        ret(3 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) = mdx_ij;
        ret(4 * n_shape_fct_ + i, 6 * n_shape_fct_ + j) = mdx_ij;
        ret(5 * n_shape_fct_ + i, 5 * n_shape_fct_ + j) = mdx_ij;
        ret(5 * n_shape_fct_ + i, j) = 2 * parameters::kappa * mdx_ij;
        ret(5 * n_shape_fct_ + i, 2 * n_shape_fct_ + j) = -mdx_ij;
        ret(2 * n_shape_fct_ + i, j) = mdy_ij;
        ret(5 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) = mdy_ij;

        //bulk triple products
        lSol_float_t hu = 0, hq = 0;
        for (unsigned int k = 0; k < j; ++k)
        {
          grad = integrator::template integrate_vol_nablaphiphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
              i, j, k, hyper_edge.geometry);
          lSol_float_t tv_shx_sh_sh_ijk = grad[0];
          hq += tv_shx_sh_sh_ijk * ca[n_shape_fct_ + k];
          hu += tv_shx_sh_sh_ijk * ca[k];
          ret(n_shape_fct_ + i, k) += tv_shx_sh_sh_ijk * ca[n_shape_fct_ + j];
          ret(n_shape_fct_ + i, n_shape_fct_ + k) += tv_shx_sh_sh_ijk * ca[j];
        }
        {
          grad = integrator::template integrate_vol_nablaphiphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
              i, j, j, hyper_edge.geometry);
          lSol_float_t tv_shx_sh_sh_ijj = grad[0];
          hq += tv_shx_sh_sh_ijj * ca[n_shape_fct_ + j];
          hu += tv_shx_sh_sh_ijj * ca[j];
        }
        ret(n_shape_fct_ + i, j) += hq;
        ret(n_shape_fct_ + i, n_shape_fct_ + j) += hu;
      }
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        ret(5 * n_shape_fct_ + i, j) += 2 * onepointfive * ret(n_shape_fct_ + i, n_shape_fct_ + j);
        ret(5 * n_shape_fct_ +i, n_shape_fct_ + j) = ret(n_shape_fct_ + i, j);
      }
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        SmallVec<hyEdge_dim(), lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          lSol_float_t b_sh_sh_ij = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
            i, j, bdr, hyper_edge.geometry);
          ret(n_shape_fct_ + i, n_shape_fct_ + j) += tau_uqq_ * b_sh_sh_ij
            * normal[0] * normal[0];
          ret(5 * n_shape_fct_ + i, 5 * n_shape_fct_ + j) -= b_sh_sh_ij * normal[0];
          ret(5 * n_shape_fct_ + i, 2 * n_shape_fct_ + j) -= -b_sh_sh_ij * normal[0];
          ret(5 * n_shape_fct_ + i, j) -= 2 * parameters::kappa * b_sh_sh_ij * normal[0];
          ret(5 * n_shape_fct_ + i, j) -= tau_f_ * b_sh_sh_ij * normal[0] * normal[0];

          if (normal[1] * normal[1] < eps && normal[0] < 0) //V+
          {
            ret(3 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) -= b_sh_sh_ij * normal[0];
            ret(5 * n_shape_fct_ + i, j) -= -(tau_pzu_ - tau_ppu_) * b_sh_sh_ij * normal[0] * normal[0];
          }
          else if (normal[1] * normal[1] < eps && normal[0] > 0) //V+
          {
            ret(5 * n_shape_fct_ + i, j) -= -(tau_mzu_ - tau_mpu_) * b_sh_sh_ij * normal[0] * normal[0];
            ret(5 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) -= -(tau_mzv_ - tau_mpv_) * b_sh_sh_ij 
              * normal[0] * normal[0];
          }
          else if (normal[0] * normal[0] < eps && normal[1] < 0) //H+
          {
            ret(5 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) -= b_sh_sh_ij * normal[1];
          }
          else if (normal[0] * normal[0] < eps && normal[1] > 0) //H-
          {
            ret(2 * n_shape_fct_ + i, j) -= b_sh_sh_ij * normal[1];
          }
        }
      }
    }

    for (unsigned int bdr = 0; bdr < 2; ++bdr)
    {
      SmallVec<hyEdge_dim(), lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
      std::array<lSol_float_t, n_shape_fct_ * n_shape_fct_> u_int, q_int;
      u_int.fill(0.);
      q_int.fill(0.);
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        //i > j > k
        for (unsigned int j = 0; j < i; ++j)
        {
          for (unsigned int k = 0; k < j; ++k)
          {
            lSol_float_t tb_sh_sh_sh_ijk = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, k, j, bdr, hyper_edge.geometry);
            q_int[n_shape_fct_ * i + j] += tb_sh_sh_sh_ijk * ca[n_shape_fct_ + k] * normal[0];
            q_int[n_shape_fct_ * i + k] += tb_sh_sh_sh_ijk * ca[n_shape_fct_ + j] * normal[0];
            q_int[n_shape_fct_ * j + k] += tb_sh_sh_sh_ijk * ca[n_shape_fct_ + i] * normal[0];
            q_int[n_shape_fct_ * j + i] += tb_sh_sh_sh_ijk * ca[n_shape_fct_ + k] * normal[0];
            q_int[n_shape_fct_ * k + i] += tb_sh_sh_sh_ijk * ca[n_shape_fct_ + j] * normal[0];
            q_int[n_shape_fct_ * k + j] += tb_sh_sh_sh_ijk * ca[n_shape_fct_ + i] * normal[0];

            u_int[n_shape_fct_ * i + j] += tb_sh_sh_sh_ijk * ca[k] * normal[0];
            u_int[n_shape_fct_ * i + k] += tb_sh_sh_sh_ijk * ca[j] * normal[0];
            u_int[n_shape_fct_ * j + k] += tb_sh_sh_sh_ijk * ca[i] * normal[0];
            u_int[n_shape_fct_ * j + i] += tb_sh_sh_sh_ijk * ca[k] * normal[0];
            u_int[n_shape_fct_ * k + i] += tb_sh_sh_sh_ijk * ca[j] * normal[0];
            u_int[n_shape_fct_ * k + j] += tb_sh_sh_sh_ijk * ca[i] * normal[0];
          }

        //i > j = k
          {
            lSol_float_t tb_sh_sh_sh_ijj = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, j, j, bdr, hyper_edge.geometry);
            q_int[n_shape_fct_ * i + j] += tb_sh_sh_sh_ijj * ca[n_shape_fct_ + j] * normal[0];
            q_int[n_shape_fct_ * j + i] += tb_sh_sh_sh_ijj * ca[n_shape_fct_ + j] * normal[0];
            q_int[n_shape_fct_ * j + j] += tb_sh_sh_sh_ijj * ca[n_shape_fct_ + i] * normal[0];

            u_int[n_shape_fct_ * i + j] += tb_sh_sh_sh_ijj * ca[j] * normal[0];
            u_int[n_shape_fct_ * j + i] += tb_sh_sh_sh_ijj * ca[j] * normal[0];
            u_int[n_shape_fct_ * j + j] += tb_sh_sh_sh_ijj * ca[i] * normal[0];
          }
        }
        //i = j > k
        for (unsigned int k = 0; k < i; ++k)
        {
          {
            lSol_float_t tb_sh_sh_sh_iik = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, i, k, bdr, hyper_edge.geometry);
            q_int[n_shape_fct_ * i + i] += tb_sh_sh_sh_iik * ca[n_shape_fct_ + k] * normal[0];
            q_int[n_shape_fct_ * i + k] += tb_sh_sh_sh_iik * ca[n_shape_fct_ + i] * normal[0];
            q_int[n_shape_fct_ * k + i] += tb_sh_sh_sh_iik * ca[n_shape_fct_ + i] * normal[0];

            u_int[n_shape_fct_ * i + i] += tb_sh_sh_sh_iik * ca[k] * normal[0];
            u_int[n_shape_fct_ * i + k] += tb_sh_sh_sh_iik * ca[i] * normal[0];
            u_int[n_shape_fct_ * k + i] += tb_sh_sh_sh_iik * ca[i] * normal[0];
          }
        }
        //i = j = k
        {
          {
            lSol_float_t tb_sh_sh_sh_iii = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, i, i, bdr, hyper_edge.geometry);
            q_int[n_shape_fct_ * i + i] += tb_sh_sh_sh_iii * ca[n_shape_fct_ + i] * normal[0];
            u_int[n_shape_fct_ * i + i] += tb_sh_sh_sh_iii * ca[i] * normal[0];
          }
        }
      }
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          ret(n_shape_fct_ + i, j) -= 0.5 * q_int[n_shape_fct_ * i + j];
          ret(n_shape_fct_ + i, n_shape_fct_ + j) -= 0.5 * u_int[n_shape_fct_ * i + j];
          ret(5 * n_shape_fct_ + i, j) -= 2 * onepointfive * u_int[n_shape_fct_ * i + j];
        }

        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          lSol_float_t hq = 0;
          for (unsigned int k = 0; k < n_shape_bdr_; ++k)
          {
            lSol_float_t tb_sh_sh_sk_ijk = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
              i, j, k, bdr, hyper_edge.geometry);
            hq += tb_sh_sh_sk_ijk * lambda_values[bdr][n_shape_bdr_ + k] * normal[0];
          }
          ret(n_shape_fct_ + i, j) -= 0.5 * hq;
        }
      }
    }
    return ret;
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> residual_lambda_directional_derivative(const SmallMatT& lambda_values,
                                                                                    const SmallVec<n_loc_dofs_, lSol_float_t>& coeff,
                                                                                    const SmallMatT& lambda_dir,
                                                                                    hyEdgeT& hyper_edge,
                                                                                    const lSol_float_t time) const
  {
    SmallVec<n_loc_dofs_, lSol_float_t> grad(0.);
    const lSol_float_t eps = ldexp(1., -40);
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    for (unsigned int bdr = 0; bdr < 4; ++bdr)
    {
      SmallVec<hyEdge_dim(), lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
      if (normal[1] * normal[1] < eps && normal[0] > 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //V^-
      {
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            const lSol_float_t b_shsk_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            //first eq
            grad[i] -= normal[0] * b_shsk_ij * lambda_dir[bdr][j];
            //second eq cf below
            //third eq does not contribute here
            //fourth eq
            grad[3 * n_shape_fct_ + i] -= normal[0] * b_shsk_ij * lambda_dir[bdr][2 * n_shape_bdr_ + j];
            //fifth eq
            grad[4 * n_shape_fct_ + i] += normal[0] * (1. / delta_t_) * b_shsk_ij * lambda_dir[bdr][j];
            //z^, f^ and p^ of sixth eq
            grad[5 * n_shape_fct_ + i] -= normal[0] * normal[0] * (tau_mzu_ - tau_mpu_ - tau_f_) * b_shsk_ij * lambda_dir[bdr][j];
            grad[5 * n_shape_fct_ + i] -= normal[0] * normal[0] * (tau_mzv_ - tau_mpv_) * b_shsk_ij * lambda_dir[bdr][2 * n_shape_bdr_ + j];
            //seventh eq is independent of lambda
            //contribution of uq^
            lSol_float_t uqc = normal[0] * normal[0] * tau_uqq_ * b_shsk_ij * lambda_dir[bdr][n_shape_bdr_ + j];
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
            {
              const lSol_float_t b_shshsk_ikj = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
                  i, k, j, bdr, hyper_edge.geometry);
              uqc += normal[0] * 0.5 * b_shshsk_ikj * coeff[k] * lambda_dir[bdr][n_shape_bdr_ + j];
            }
            grad[n_shape_fct_ + i] -= uqc;
            //contribution of q^2^
            lSol_float_t q2c = 0;
            for (unsigned int k = 0; k < n_shape_bdr_; ++k)
            {
              const lSol_float_t b_shsksk_ikj = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
                  i, k, j, bdr, hyper_edge.geometry);
              q2c += normal[0] * b_shsksk_ikj * lambda_values[bdr][n_shape_bdr_ + k] * lambda_dir[bdr][n_shape_bdr_ + j];
            }
            grad[5 * n_shape_fct_ + i] -= q2c;
          }
        }
      }
      else if (normal[1] * normal[1] < eps && normal[0] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //V^+
      {
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            const lSol_float_t b_shsk_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            //first eq
            grad[i] -= normal[0] * b_shsk_ij * lambda_dir[bdr][j];
            //second eq cf below
            //third eq does not contribute here
            //fourth eq does not contribute here (v^ = v)
            //fifth eq
            grad[4 * n_shape_fct_ + i] += normal[0] * (1. / delta_t_) * b_shsk_ij * lambda_dir[bdr][j];
            //z^, f^ and p^ of sixth eq
            grad[5 * n_shape_fct_ + i] -= normal[0] * normal[0] * ((tau_pzu_ - tau_ppu_) - tau_f_) * b_shsk_ij * lambda_dir[bdr][j];
            //seventh eq is independent of lambda
            //contribution of uq^
            lSol_float_t uqc = normal[0] * normal[0] * tau_uqq_ * b_shsk_ij * lambda_dir[bdr][n_shape_bdr_ + j];
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
            {
              const lSol_float_t b_shshsk_ikj = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
                  i, k, j, bdr, hyper_edge.geometry);
              uqc += normal[0] * 0.5 * b_shshsk_ikj * coeff[k] * lambda_dir[bdr][n_shape_bdr_ + j];
            }
            grad[n_shape_fct_ + i] -= uqc;
            //contribution of q^2^
            lSol_float_t q2c = 0;
            for (unsigned int k = 0; k < n_shape_bdr_; ++k)
            {
              const lSol_float_t b_shsksk_ikj = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
                  i, k, j, bdr, hyper_edge.geometry);
              q2c += normal[0] * b_shsksk_ikj * lambda_values[bdr][n_shape_bdr_ + k] * lambda_dir[bdr][n_shape_bdr_ + j];
            }
            grad[5 * n_shape_fct_ + i] -= q2c;
          }
        }
      }
      else if (normal[0] * normal[0] < eps && normal[1] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //H+
      {
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            const lSol_float_t b_shsk_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            //first two equations don't contribute here
            //third equation
            grad[2 * n_shape_fct_ + i] -= normal[1] * b_shsk_ij * lambda_dir[bdr][j];
            //fourth and fifth equations don't contribute here
            //sixth equation doesn't contribute here (v^ = v)
            //seventh equation is independent of lambda
          }
        }
      }
      else if (normal[0] * normal[0] < eps && normal[1] > 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //H-
      {
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            const lSol_float_t b_shsk_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            //first two equations don't contribute here
            //third equation does not contribute here (u^ = u)
            //fourth and fifth equations don't contribute here
            //sixth equation
            grad[5 * n_shape_fct_ + i] -= normal[1] * b_shsk_ij * lambda_dir[bdr][2 * n_shape_bdr_ + j];
            //seventh equation is independent of lambda
          }
        }
      }
    }
    return grad;
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
  inline SmallVec<n_loc_dofs_, lSol_float_t>& get_residual(const SmallMatT& lambda_values,
                                                            const SmallVec<n_loc_dofs_, lSol_float_t>& ca,
                                                            SmallVec<n_loc_dofs_, lSol_float_t>& residual,
                                                            hyEdgeT& hyper_edge,
                                                            const lSol_float_t time,
                                                            const bool print=false) const
  {
    static_assert(std::is_same<typename SmallMatT::value_type::value_type, lSol_float_t>::value,
        "Lambda values ...");
    hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
        "The size ...");
    for (unsigned int i = 0; i < hyEdge_dimT; ++i)
      hy_assert(lambda_values[i].size() == 3 * n_shape_bdr_,
          "The size of ...");

    residual *= 0.;
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    SmallVec<hyEdge_dim(), lSol_float_t> grad;

    const lSol_float_t eps = ldexp(1., -40);

    for (unsigned int i = 0; i < n_shape_fct_; ++i) 
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        //mass integrals
        const lSol_float_t mass_ij = integrator::template integrate_vol_phiphi<decltype(hyEdgeT::geometry)>
          (i, j, hyper_edge.geometry);
        residual[i] += mass_ij * ca[n_shape_fct_ + j];
        residual[n_shape_fct_ + i] += mass_ij * ca[2 * n_shape_fct_ + j];
        residual[2 * n_shape_fct_ + i] += mass_ij * ca[3 * n_shape_fct_ + j];
        residual[3 * n_shape_fct_ + i] += mass_ij * ca[3 * n_shape_fct_ + j];
        residual[4 * n_shape_fct_ + i] += mass_ij * ca[5 * n_shape_fct_ + j];
        residual[5 * n_shape_fct_ + i] += mass_ij * ca[6 * n_shape_fct_ + j];
        residual[6 * n_shape_fct_ + i] += (1. / delta_t_) * mass_ij * (ca[j] - hyper_edge.data.u_old[j]);
        residual[6 * n_shape_fct_ + i] += mass_ij * ca[6 * n_shape_fct_ + j];
        
        //integrals involving derivatives
        grad =  integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, 
             decltype(hyEdgeT::geometry)>(i, j, hyper_edge.geometry);
        const lSol_float_t mdx_ij = grad[0], mdy_ij = grad[1];
        residual[i] += mdx_ij * ca[j];
        residual[2 * n_shape_fct_ + i] += mdy_ij * ca[j];
        residual[3 * n_shape_fct_ + i] += mdx_ij * ca[4 * n_shape_fct_ + j];
        residual[4 * n_shape_fct_ + i] += mdx_ij * ca[6 * n_shape_fct_ + j];
        residual[5 * n_shape_fct_ + i] += mdx_ij * ca[5 * n_shape_fct_ + j];
        residual[5 * n_shape_fct_ + i] += 2 * parameters::kappa * mdx_ij * ca[j];
        residual[5 * n_shape_fct_ + i] += -mdx_ij * ca[2 * n_shape_fct_ + j];
        residual[5 * n_shape_fct_ + i] += mdy_ij * ca[4 * n_shape_fct_ + j];

        //bulk triple products
        lSol_float_t q2_int = 0, u2_int = 0, uq_int = 0;
        for (unsigned int k = 0; k < j; ++k)
        {
          const lSol_float_t tv_ix_jk = integrator::template integrate_vol_nablaphiphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
              i, j, k, hyper_edge.geometry)[0];
          q2_int += 2 * tv_ix_jk * ca[n_shape_fct_ + j] * ca[n_shape_fct_ + k];
          uq_int += tv_ix_jk * (ca[n_shape_fct_ + j] * ca[k] + ca[n_shape_fct_ + k] * ca[j]);
          u2_int += 2 * tv_ix_jk * ca[j] * ca[k];
        }
        const lSol_float_t tv_ix_jj = integrator::template integrate_vol_nablaphiphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
            i, j, j, hyper_edge.geometry)[0];
        q2_int += tv_ix_jj * ca[n_shape_fct_ + j] * ca[n_shape_fct_ + j];
        uq_int += tv_ix_jj * ca[n_shape_fct_ + j] * ca[j];
        u2_int += tv_ix_jj * ca[j] * ca[j];
        residual[n_shape_fct_ + i] += uq_int;
        residual[5 * n_shape_fct_ + i] += onepointfive * u2_int;
        residual[5 * n_shape_fct_ + i] += 0.5 * q2_int;
      }

      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        //fixed integrals for notational convenience
        SmallVec<hyEdge_dim(), lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
        lSol_float_t u_int = 0, q_int = 0, uo_int = 0, p_int = 0, v_int = 0, z_int = 0;
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          lSol_float_t h = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
              i, j, bdr, hyper_edge.geometry);
          u_int += h * ca[j];
          uo_int += h * hyper_edge.data.u_old[j];
          q_int += h * ca[1 * n_shape_fct_ + j];
          p_int += h * ca[2 * n_shape_fct_ + j];
          v_int += h * ca[4 * n_shape_fct_ + j];
          z_int += h * ca[5 * n_shape_fct_ + j];
        }
        lSol_float_t uh_int = 0, uho_int = 0, qh_int = 0, vh_int = 0;
        for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        {
          lSol_float_t h = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
              i, j, bdr, hyper_edge.geometry);
          uh_int += h * lambda_values[bdr][j];
          uho_int += h * hyper_edge.data.uh_old[bdr][j];
          qh_int += h * lambda_values[bdr][n_shape_bdr_ + j];
          vh_int += h * lambda_values[bdr][2 * n_shape_bdr_ + j];
        }
        lSol_float_t flux_ux = (uh_int - u_int) * normal[0];
        lSol_float_t flux_uy = (uh_int - u_int) * normal[1];
        lSol_float_t flux_q = (qh_int - q_int) * normal[0];
        lSol_float_t flux_v = (vh_int - v_int) * normal[0];
        lSol_float_t trace_f = 0, trace_q2 = 0, trace_uq = 0;
        
        trace_f += 2 * parameters::kappa * u_int - tau_f_ * flux_ux;
        trace_uq += tau_uqq_ * flux_q;
        lSol_float_t trace_p = 0, trace_vv = 0, trace_vh = 0, trace_z = 0, trace_r = 0;
        lSol_float_t trace_zp = 0., trace_uh = 0.;
        if (normal[1] * normal[1] < eps && normal[0] < 0) //V+
        {
          trace_zp = z_int - p_int + (tau_pzu_ - tau_ppu_) * flux_ux;
          trace_vv = v_int;
          trace_r = -1. / delta_t_ * (uh_int - uho_int);
        }
        else if (normal[1] * normal[1] < eps && normal[0] > 0) //V-
        {
          trace_zp = z_int - p_int + (tau_mzu_ - tau_mpu_) * flux_ux + (tau_mzv_ - tau_mpv_) * flux_v;
          trace_vv = vh_int;
          trace_r = -1. / delta_t_ * (uh_int - uho_int);
        }
        else if (normal[0] * normal[0] < eps && normal[1] < 0) //H+
        {
          trace_vh = v_int;
          trace_uh = uh_int;
        }
        else if (normal[0] * normal[0] < eps && normal[1] > 0) //H-
        {
          trace_vh = vh_int;
          trace_uh = u_int;
        }

        residual[i] -= uh_int * normal[0];
        residual[n_shape_fct_ + i] -= trace_uq * normal[0];
        residual[2 * n_shape_fct_ + i] -= trace_uh * normal[1];
        residual[3 * n_shape_fct_ + i] -= trace_vv * normal[0];
        residual[4 * n_shape_fct_ + i] -= trace_r * normal[0];
        residual[5 * n_shape_fct_ + i] -= (trace_zp + trace_f + 0.5 * trace_q2) * normal[0];
        residual[5 * n_shape_fct_ + i] -= trace_vh * normal[1];
      }

      lSol_float_t rhs = integrator::template integrate_vol_phifunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side, Point<hyEdge_dimT, lSol_float_t> > (i, hyper_edge.geometry, time);
      residual[5 * n_shape_fct_ + i] += rhs;
    }
    for (unsigned int bdr = 0; bdr < 2; ++bdr)
    {
      const int nx = -1 + 2 * bdr;
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        
        for (unsigned int j = 0; j < i; ++j)
        {
          for (unsigned int k = 0; k < j; ++k)
          {
            const lSol_float_t sss_ijk = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
                i, j, k, bdr, hyper_edge.geometry);
            residual[5 * n_shape_fct_ + i] -= onepointfive * 2 * sss_ijk * ca[j] * ca[k] * nx;
            residual[5 * n_shape_fct_ + j] -= onepointfive * 2 * sss_ijk * ca[i] * ca[k] * nx;
            residual[5 * n_shape_fct_ + k] -= onepointfive * 2 * sss_ijk * ca[j] * ca[i] * nx;
            residual[1 * n_shape_fct_ + i] -= .5 * sss_ijk * (ca[j] * ca[n_shape_fct_ + k] + ca[k] * ca[n_shape_fct_ + j]) * nx;
            residual[1 * n_shape_fct_ + j] -= .5 * sss_ijk * (ca[i] * ca[n_shape_fct_ + k] + ca[k] * ca[n_shape_fct_ + i]) * nx;
            residual[1 * n_shape_fct_ + k] -= .5 * sss_ijk * (ca[j] * ca[n_shape_fct_ + i] + ca[i] * ca[n_shape_fct_ + j]) * nx;
          }
          const lSol_float_t sss_ijj = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
                i, j, j, bdr, hyper_edge.geometry);
          residual[5 * n_shape_fct_ + i] -= onepointfive * sss_ijj * ca[j] * ca[j] * nx;
          residual[5 * n_shape_fct_ + j] -= onepointfive * 2 * sss_ijj * ca[i] * ca[j] * nx;
          residual[1 * n_shape_fct_ + i] -= .5 * sss_ijj * ca[j] * ca[n_shape_fct_ + j] * nx;
          residual[1 * n_shape_fct_ + j] -= .5 * sss_ijj * (ca[i] * ca[n_shape_fct_ + j] + ca[j] * ca[n_shape_fct_ + i]) * nx;
        }
        for (unsigned int k = 0; k < i; ++k)
        {
          const lSol_float_t sss_iik = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
                i, i, k, bdr, hyper_edge.geometry);
          residual[5 * n_shape_fct_ + k] -= onepointfive * sss_iik * ca[i] * ca[i] * nx;
          residual[5 * n_shape_fct_ + i] -= onepointfive * 2 * sss_iik * ca[i] * ca[k] * nx;
          residual[1 * n_shape_fct_ + k] -= .5 * sss_iik * ca[i] * ca[n_shape_fct_ + i] * nx;
          residual[1 * n_shape_fct_ + i] -= .5 * sss_iik * (ca[i] * ca[n_shape_fct_ + k] + ca[k] * ca[n_shape_fct_ + i]) * nx;
        }
        const lSol_float_t sss_iii = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, i, i, bdr, hyper_edge.geometry);
        residual[5 * n_shape_fct_ + i] -= onepointfive * sss_iii * ca[i] * ca[i] * nx;
        residual[1 * n_shape_fct_ + i] -= .5 * sss_iii * ca[i] * ca[n_shape_fct_ + i] * nx;
        
        for (unsigned int j = 0; j < i; ++j)
        {
          for (unsigned int k = 0; k < n_shape_bdr_; ++k)
          {
            const lSol_float_t ssb_ijk = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
                i, j, k, bdr, hyper_edge.geometry);
            residual[1 * n_shape_fct_ + i] -= .5 * ssb_ijk * ca[j] * lambda_values[bdr][n_shape_bdr_ + k] * nx;
            residual[1 * n_shape_fct_ + j] -= .5 * ssb_ijk * ca[i] * lambda_values[bdr][n_shape_bdr_ + k] * nx;
          }
        }
        for (unsigned int k = 0; k < n_shape_bdr_; ++k)
        {
          const lSol_float_t ssb_iik = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
              i, i, k, bdr, hyper_edge.geometry);
          residual[1 * n_shape_fct_ + i] -= .5 * ssb_iik * ca[i] * lambda_values[bdr][n_shape_bdr_ + k] * nx;
        }
        for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        {
          for (unsigned int k = 0; k < j; ++k)
          {
            const lSol_float_t sbb_ijk = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
                i, j, k, bdr, hyper_edge.geometry);
            residual[5 * n_shape_fct_ + i] -= 0.5 * 2 * sbb_ijk * lambda_values[bdr][n_shape_bdr_ + j] * lambda_values[bdr][n_shape_bdr_ + k] * nx;
          }
          const lSol_float_t sbb_ijj = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
              i, j, j, bdr, hyper_edge.geometry);
          residual[5 * n_shape_fct_ + i] -= 0.5 *  sbb_ijj * lambda_values[bdr][n_shape_bdr_ + j] * lambda_values[bdr][n_shape_bdr_ + j] * nx;
        }
      }
    }
 
    return residual;
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallMatT coupling_lambda_directional_derivative(const SmallMatT& lambda_values,
                                                          const SmallVec<n_loc_dofs_, lSol_float_t> coeff,
                                                          const SmallMatT& lambda_dir,
                                                          SmallMatT& out,
                                                          hyEdgeT& hyper_edge,
                                                          const lSol_float_t time) const
  {
    lSol_float_t eps = ldexp(1., -40);
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      SmallVec<2, lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
      using parameters = parametersT<hyEdge_dim(), lSol_float_t>;
      //delete
      //out[bdr].fill(0.);

      if (normal[1] * normal[1] < eps && normal[0] > 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            lSol_float_t b_int = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                j, i, bdr, hyper_edge.geometry);
            out[bdr][i] += normal[0] * normal[0] * (tau_mzu_ - tau_mpu_ - tau_f_) * b_int * lambda_dir[bdr][j];
            out[bdr][i] += normal[0] * normal[0] * (tau_mzv_ - tau_mpv_) * b_int * lambda_dir[bdr][2 * n_shape_bdr_ + j];
            out[bdr][n_shape_bdr_ + i] += normal[0] * normal[0] * tau_uqq_ * b_int * lambda_dir[bdr][n_shape_bdr_ + j];
            out[bdr][2 * n_shape_bdr_ + i] += normal[0] * b_int * lambda_dir[bdr][2 * n_shape_bdr_ + j];
            lSol_float_t h = 0;
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
            {
              lSol_float_t t_int = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
                  k, j, i, bdr, hyper_edge.geometry);
              h += t_int * coeff[k];
            }
            out[bdr][n_shape_bdr_ + i] += .5 * h * lambda_dir[bdr][n_shape_bdr_ + j] * normal[0];
          }
        }
      }
      else if (normal[1] * normal[1] < eps && normal[0] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            lSol_float_t b_int = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                j, i, bdr, hyper_edge.geometry);
            out[bdr][i] += normal[0] * normal[0] * (tau_pzu_ - tau_ppu_ - tau_f_ ) * b_int * lambda_dir[bdr][j];
            out[bdr][n_shape_bdr_ + i] += normal[0] * normal[0] * tau_uqq_ * b_int * lambda_dir[bdr][n_shape_bdr_ + j];
            lSol_float_t h = 0;
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
            {
              lSol_float_t t_int = integrator::template integrate_bdr_phipsipsi<decltype(hyEdgeT::geometry)>(
                  k, j, i, bdr, hyper_edge.geometry);
              h += t_int * coeff[k];
            }
            out[bdr][n_shape_bdr_ + i] += .5 * h * lambda_dir[bdr][n_shape_bdr_ + j] * normal[0];
          }
        }
      }
      else if (normal[0] * normal[0] < eps && normal[1] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //H+
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            lSol_float_t b_int = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                j, i, bdr, hyper_edge.geometry);
            out[bdr][n_shape_bdr_ + i] += normal[1] * b_int * lambda_dir[bdr][j];
          }
        }
      }
      else if (normal[0] * normal[0] < eps && normal[1] > 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //H-
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            lSol_float_t b_int = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                j, i, bdr, hyper_edge.geometry);
            out[bdr][i] += normal[1] * b_int * lambda_dir[bdr][2 * n_shape_bdr_ + j];
          }
        }
      }

    }
    return out;
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallMatT coupling_coeff_directional_derivative(const SmallMatT& lambda_values,
                                                         const SmallVec<n_loc_dofs_, lSol_float_t>& coeff,
                                                         const SmallVec<n_loc_dofs_, lSol_float_t>& coeff_dir,
                                                         SmallMatT& out,
                                                         hyEdgeT& hyper_edge,
                                                         const lSol_float_t time) const
  {
    lSol_float_t eps = ldexp(1., -40);
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      SmallVec<2, lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
      using parameters = parametersT<hyEdge_dim(), lSol_float_t>;
      //delete
      //out[bdr].fill(0.);
      if (normal[1] * normal[1] < eps && normal[0] > 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //V-
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_fct_; ++j)
          {
            lSol_float_t b_sk_sh_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
              j, i, bdr, hyper_edge.geometry);
            //contributions of p and z
            out[bdr][i] += normal[0] * b_sk_sh_ij * (coeff_dir[5 * n_shape_fct_ + j] - coeff_dir[2 * n_shape_fct_ + j]);
            //contributions of the fluxes
            out[bdr][i] -= normal[0] * normal[0] * b_sk_sh_ij * (tau_mzu_ - tau_mpu_ - tau_f_) * coeff_dir[j];
            out[bdr][i] -= normal[0] * normal[0] * b_sk_sh_ij * (tau_mzv_ - tau_mpv_) * coeff_dir[4 * n_shape_fct_ + j];
            out[bdr][n_shape_bdr_ + i] -= normal[0] * normal[0] * b_sk_sh_ij * tau_uqq_ * coeff_dir[n_shape_fct_ + j];
            //contribution of f and uq
            lSol_float_t fc = 2 * parameters::kappa * b_sk_sh_ij * coeff_dir[j], uqc = 0;
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
            {
              lSol_float_t b_sk_sh_sh_ikj = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
                j, k, i, bdr, hyper_edge.geometry);
              fc += 2 * onepointfive * b_sk_sh_sh_ikj * coeff[j] * coeff_dir[k];
              uqc += 0.5 * b_sk_sh_sh_ikj * (coeff[n_shape_fct_ + k] * coeff_dir[j] + coeff[k] * coeff_dir[n_shape_fct_ + j]);
            }
            for (unsigned int k = 0; k < n_shape_bdr_; ++k)
            {
              lSol_float_t b_sk_sk_sh_ikj = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
                j, k, i, bdr, hyper_edge.geometry);
              uqc += 0.5 * b_sk_sk_sh_ikj * lambda_values[bdr][n_shape_bdr_ + k] * coeff_dir[j];
            }
            out[bdr][i] += normal[0] * fc;
            out[bdr][n_shape_bdr_ + i] += normal[0] * uqc;
          }
        }
      }
      else if (normal[1] * normal[1] < eps && normal[0] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //V+
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_fct_; ++j)
          {
            lSol_float_t b_sk_sh_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
              j, i, bdr, hyper_edge.geometry);
            //contributions of p and z
            out[bdr][i] += normal[0] * b_sk_sh_ij * (coeff_dir[5 * n_shape_fct_ + j] - coeff_dir[2 * n_shape_fct_ + j]);
            //contribution of v
            out[bdr][2 * n_shape_bdr_ + i] += normal[0] * b_sk_sh_ij * coeff_dir[4 * n_shape_fct_ + j];
            //contributions of the fluxes
            out[bdr][i] -= normal[0] * normal[0] * b_sk_sh_ij * (tau_pzu_ - tau_ppu_ - tau_f_) * coeff_dir[j];
            out[bdr][n_shape_bdr_ + i] -= normal[0] * normal[0] * b_sk_sh_ij * tau_uqq_ * coeff_dir[n_shape_fct_ + j];
            //contribution of f and uq
            lSol_float_t fc = 2 * parameters::kappa * b_sk_sh_ij * coeff_dir[j], uqc = 0;
            for (unsigned int k = 0; k < n_shape_fct_; ++k)
            {
              lSol_float_t b_sk_sh_sh_ikj = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
                j, k, i, bdr, hyper_edge.geometry);
              fc += 2 * onepointfive * b_sk_sh_sh_ikj * coeff[j] * coeff_dir[k];
              uqc += 0.5 * b_sk_sh_sh_ikj * (coeff[n_shape_fct_ + k] * coeff_dir[j] + coeff[k] * coeff_dir[n_shape_fct_ + j]);
            }
            for (unsigned int k = 0; k < n_shape_bdr_; ++k)
            {
              lSol_float_t b_sk_sk_sh_ikj = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
                j, k, i, bdr, hyper_edge.geometry);
              uqc += 0.5 * b_sk_sk_sh_ikj * lambda_values[bdr][n_shape_bdr_ + k] * coeff_dir[j];
            }
            out[bdr][i] += normal[0] * fc;
            out[bdr][n_shape_bdr_ + i] += normal[0] * uqc;
          }
        }
      }
      else if (normal[0] * normal[0] < eps && normal[1] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //H+
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_fct_; ++j)
          {
            lSol_float_t b_sk_sh_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
              j, i, bdr, hyper_edge.geometry);
            //contribution of v
            out[bdr][i] += normal[1] * b_sk_sh_ij * coeff_dir[4 * n_shape_fct_ + j];
          }
        }
      }
      else if (normal[0] * normal[0] < eps && normal[1] > 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr])) //H-
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_fct_; ++j)
          {
            lSol_float_t b_sk_sh_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
              j, i, bdr, hyper_edge.geometry);
            //contribution of u
            out[bdr][n_shape_bdr_ + i] += normal[1] * b_sk_sh_ij * coeff_dir[j];
          }
        }
      }
    }
    return out;
  }

  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& trace_to_flux(const SmallMatInT& lambda_values_in_uc,
                              const SmallMatInT& lambda_values_dir,
                              SmallMatOutT& lambda_values_out,
                              hyEdgeT& hyper_edge,
                              const lSol_float_t time) const
  {
    //ensure dirichlet conditions are met
    SmallMatInT lambda_values_in = lambda_values_in_uc;
    make_skeleton(lambda_values_in, hyper_edge, time);

    //compute local coefficients
    SmallVec<n_loc_dofs_, lSol_float_t> coeff(0.);
    newton(lambda_values_in, coeff, hyper_edge, time);

    //compute derivative of newton
    SmallVec<n_loc_dofs_, lSol_float_t> implicit_derivative = -1. * residual_lambda_directional_derivative(lambda_values_in,
        coeff, lambda_values_dir, hyper_edge, time) / jacobi(lambda_values_in, coeff, hyper_edge, time);
    
    //compute derivative
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      lambda_values_out[bdr].fill(0.);
    }
    coupling_coeff_directional_derivative(lambda_values_in, coeff, implicit_derivative, lambda_values_out, hyper_edge, time);
    coupling_lambda_directional_derivative(lambda_values_in, coeff, lambda_values_dir, lambda_values_out, hyper_edge, time);
    
    return lambda_values_out;
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
    //calculate coefficients
    SmallVec<n_loc_dofs_, lSol_float_t> coeff(0.);
    newton(lambda_values_in, coeff, hyper_edge, time);
    //call residual function
    return coupling_function(lambda_values_in, coeff, lambda_values_out, hyper_edge, time);
  }

  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& coupling_function(const SmallMatInT& lambda_values_in,
                                  const SmallVec<n_loc_dofs_, lSol_float_t>& coeff,
                                  SmallMatOutT& lambda_values_out,
                                  hyEdgeT& hyper_edge,
                                  const lSol_float_t time) const
  {
    std::array<SmallVec<hyEdge_dim(), lSol_float_t>, 2 * hyEdge_dim()> loc_normal;
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr) 
      loc_normal[bdr] = hyper_edge.geometry.local_normal(bdr);
    const lSol_float_t eps = ldexp(1., -30);
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      using parameters = parametersT<hyEdge_dim(), lSol_float_t>;
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
      	  lambda_values_out[bdr][i] += (tau_mzu_ - tau_mpu_) * (uh_int - u_int) * loc_normal[bdr][0];
      	  lambda_values_out[bdr][i] += (tau_mzv_ - tau_mpv_) * (vh_int - v_int) * loc_normal[bdr][0];
      	  //f_hat
      	  lambda_values_out[bdr][i] += 2 * parameters::kappa * u_int;
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
      	  {
            for (unsigned int k = 0; k < n_shape_fct_; ++k) 
	          {
              lSol_float_t c = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
      	      lambda_values_out[bdr][i] += onepointfive * c * coeff[j] * coeff[k];
	          }
	        }
          /* 
      	  std::array<lSol_float_t, n_shape_bdr_> uh_arr;
      	  std::array<lSol_float_t, n_shape_fct_> u_arr;
      	  for (unsigned int j = 0; j < n_shape_bdr_; ++j)
            uh_arr[j] = lambda_values_in[bdr][j];
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
            u_arr[j] = lambda_values_in[bdr][j];
          lambda_values_out[bdr][i] -= integrate_bdr_psicompfun<decltype(hyEdgeT::geometry), parameters::tau_f>(
		    	  i, uh_arr, u_arr, bdr, hyper_edge.geometry) * loc_normal[bdr][0];
            */
          lambda_values_out[bdr][i] -= tau_f_ * (uh_int - u_int) * loc_normal[bdr][0];
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
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] += vh_int * loc_normal[bdr][0];
      	}
      }
      if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
        //V^+ on the left
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
	      {
      	  lambda_values_out[bdr][i] = 0;
  	      lambda_values_out[bdr][n_shape_bdr_ + i] = 0;
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] = 0;
  	      lSol_float_t uh_int = 0, u_int = 0, qh_int = 0, q_int = 0, v_int = 0;
    	    for (unsigned int j = 0; j < n_shape_bdr_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
	          uh_int += c * lambda_values_in[bdr][j];
      	    qh_int += c * lambda_values_in[bdr][n_shape_bdr_ + j];
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
      	  lambda_values_out[bdr][i] += (tau_pzu_ - tau_ppu_) * (uh_int - u_int) * loc_normal[bdr][0];
      	  //f_hat
      	  lambda_values_out[bdr][i] += 2 * parameters::kappa * u_int;
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
      	  {
            for (unsigned int k = 0; k < n_shape_fct_; ++k) 
	          {
              lSol_float_t c = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
      	      lambda_values_out[bdr][i] += onepointfive * c * coeff[j] * coeff[k];
	          }
	        }
          /* 
      	  std::array<lSol_float_t, n_shape_bdr_> uh_arr;
      	  std::array<lSol_float_t, n_shape_fct_> u_arr;
      	  for (unsigned int j = 0; j < n_shape_bdr_; ++j)
            uh_arr[j] = lambda_values_in[bdr][j];
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
            u_arr[j] = lambda_values_in[bdr][j];
          lambda_values_out[bdr][i] -= integrate_bdr_psicompfun<decltype(hyEdgeT::geometry), parameters::tau_f>(
		    	  i, uh_arr, u_arr, bdr, hyper_edge.geometry) * loc_normal[bdr][0];
            */
          lambda_values_out[bdr][i] -= tau_f_ * (uh_int - u_int) * loc_normal[bdr][0];
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
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] *= loc_normal[bdr][0];
      	}
      }
      if (loc_normal[bdr][0] * loc_normal[bdr][0] < eps && loc_normal[bdr][1] < 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
      	{
  	      lambda_values_out[bdr][i] = 0;
      	  lambda_values_out[bdr][n_shape_bdr_ + i] = 0.;
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] = 0.;
      	  lSol_float_t uh_int = 0, v_int = 0;
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
	          v_int += c * coeff[4 * n_shape_fct_ + j];
      	  }
      	  lambda_values_out[bdr][i] += v_int;
      	  lambda_values_out[bdr][i] *= loc_normal[bdr][1];
      	  lambda_values_out[bdr][n_shape_bdr_ + i] += uh_int;
      	  lambda_values_out[bdr][n_shape_bdr_ + i] *= loc_normal[bdr][1];
      	}
      }
      if (loc_normal[bdr][0] * loc_normal[bdr][0] < eps && loc_normal[bdr][1] > 0 && !is_dirichlet<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
      	{
  	      lambda_values_out[bdr][i] = 0;
      	  lambda_values_out[bdr][n_shape_bdr_ + i] = 0.;
      	  lambda_values_out[bdr][2 * n_shape_bdr_ + i] = 0.;
      	  lSol_float_t u_int = 0, vh_int = 0;
      	  for (unsigned int j = 0; j < n_shape_bdr_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
	          vh_int += c * lambda_values_in[bdr][2 * n_shape_bdr_ + j];
      	  }
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
	          u_int += c * coeff[j];
      	  }
      	  lambda_values_out[bdr][i] += vh_int;
      	  lambda_values_out[bdr][i] *= loc_normal[bdr][1];
      	  lambda_values_out[bdr][n_shape_bdr_ + i] += u_int;
      	  lambda_values_out[bdr][n_shape_bdr_ + i] *= loc_normal[bdr][1];
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
  unsigned int newton(const SmallMatT& lambda_values, SmallVec<n_loc_dofs_, lSol_float_t>& coeff,
                      hyEdgeT& hyper_edge, const lSol_float_t time) const
  {
    const lSol_float_t eps = ldexp(1., -40);
    SmallVec<n_loc_dofs_, lSol_float_t> res(0.);
    res = get_residual(lambda_values, coeff, res, hyper_edge, time);
    lSol_float_t ra = norm_2(res);
    SmallSquareMat<n_loc_dofs_, lSol_float_t> jac = jacobi(lambda_values, coeff, hyper_edge, time);
    unsigned int i = 0;
    for (i = 0; ra > eps && i < 100; ++i)
    {
      lSol_float_t rn;
      lSol_float_t stepsize = 1.;
      SmallVec<n_loc_dofs_, lSol_float_t> cn = coeff;
      SmallVec<n_loc_dofs_, lSol_float_t> step = res / jac;
      do
      {
        cn = coeff - stepsize * step;
        res = get_residual(lambda_values, cn, res, hyper_edge, time);
        rn = norm_2(res);
        stepsize *= .5;
        ++i;
      } while (ra < rn && i < 100);
      coeff = cn;
      ra = rn;
      jac = jacobi(lambda_values, coeff, hyper_edge, time);
    }
    //std::cout << i << "\n";
    return i;
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
    //hy_assert(false, "Test");
    //project initial to u
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      hyper_edge.data.u_old[i] = integrator::template integrate_volUni_phifunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::initial, Point<hyEdge_dimT, lSol_float_t> > (i, hyper_edge.geometry, time);
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
        //hyper_edge.data.uh_old[bdr][i] = lambda_values[bdr][i];	-> unecessary bc. of set_skeleton_data
        if (is_neumann<parameters>(hyper_edge.node_descriptor[bdr]))
        {
          lambda_values[bdr][n_shape_bdr_ + i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::neumann_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);

        } else {
          lambda_values[bdr][n_shape_bdr_ + i] = 0.;
        }
        if (is_right<parameters>(hyper_edge.node_descriptor[bdr]))
        {
          lambda_values[bdr][2 * n_shape_bdr_ + i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::reference_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);

        } else {
          lambda_values[bdr][2 * n_shape_bdr_ + i] = 0.;
        }
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
    {
      hyper_edge.data.u_old[i] = coeff[i];
      hyper_edge.data.q_old[i] = coeff[n_shape_fct_ + i];
    }
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline void set_data(const SmallMatT& lambda_values_uc, 
      hyEdgeT& hyper_edge, const lSol_float_t time) const
  {
    //ensure dirichlet conditions
    SmallMatT lambda_values = lambda_values_uc;
    make_skeleton(lambda_values, hyper_edge, time);
    
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
        {
          lambda_values_out[bdr][i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::dirichlet_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
        }
      }
      if (is_neumann<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
          lambda_values_out[bdr][n_shape_bdr_ + i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::neumann_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
      }
      if (is_right<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
          lambda_values_out[bdr][2 * n_shape_bdr_ + i]  = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::reference_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
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
  std::array<lSol_float_t, 2U> errors(const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>,
                                                       2 * hyEdge_dimT>& lambda_values,
                                      hyEdgeT& hy_edge,
                                      const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    return std::array<lSol_float_t, 2U>({integrator::template integrate_vol_diffsquare_discana<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result, Point<hyEdge_dimT, lSol_float_t> >(hy_edge.data.u_old.data(),
                                                                      hy_edge.geometry, time),
      integrator::template integrate_vol_diffsquare_discana<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::neumann_value, Point<hyEdge_dimT, lSol_float_t> >(hy_edge.data.q_old.data(),
                                                                      hy_edge.geometry, time)});
  }
  
  /*!***********************************************************************************************
   * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
   *
   * \tparam  abscissa_float_t  Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Input array type.
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which function is plotted.
   * \retval  func_vals         Array of function values.
   ************************************************************************************************/
  template <typename abscissa_float_t,
            std::size_t abscissas_sizeT,
            class input_array_t,
            class hyEdgeT>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>, system_dim>
  bulk_values(
    const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
    const input_array_t& lambda_values_uc,
    hyEdgeT& hyper_edge,
    const lSol_float_t time) const
  {
    //ensure dirichlet conditions
    input_array_t lambda_values = lambda_values_uc;
    make_skeleton(lambda_values, hyper_edge, time);
    
    SmallVec<n_loc_dofs_, lSol_float_t> coefficients(0.);
    newton(lambda_values, coefficients, hyper_edge, time);
    SmallVec<n_shape_fct_, lSol_float_t> coeffs;
    SmallVec<static_cast<unsigned int>(abscissas_sizeT), abscissa_float_t> helper(abscissas);

    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>,
             Chkp<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::system_dim>
      point_vals;

    for (unsigned int d = 0; d < system_dim; ++d)
    {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[i + d * n_shape_fct_];
      for (unsigned int pt = 0; pt < Hypercube<hyEdge_dimT>::pow(abscissas_sizeT); ++pt)
        point_vals[d][pt] = integrator::shape_fun_t::template lin_comb_fct_val<float>(
          coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT> >(pt, helper));
    }

    return point_vals;
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
};

}
