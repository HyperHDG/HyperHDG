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
struct ZKParametersDefault
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> dirichlet_nodes{};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> right_nodes{};
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
  static param_float_t uh_value(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   q values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t qh_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   s values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t sh_value(const Point<space_dimT, param_float_t>&,
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
          template <unsigned int, typename> typename parametersT = ZKParametersDefault,
          typename lSol_float_t = double>
class ZK
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
  static constexpr unsigned int n_loc_dofs_ = 5 * n_shape_fct_;
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
  const lSol_float_t tau_mpq_;
  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_pqu_;

  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_psu_;

  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_ru_;
    /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_f_;
  /*!***********************************************************************************************
   * \brief   An integrator helps to easily evaluate integrals (e.g. via quadrature).
   ************************************************************************************************/
  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT> >,
    lSol_float_t>
    integrator;

  const lSol_float_t onehalf = 0.5; 
   
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
  ZK(const constructor_value_type& constru = std::vector<double>({1., -1., -1., 1., 1., 1.,  -3., 4.}))
  : delta_t_(constru[0]), tau_ppu_(constru[1]), tau_mpu_(constru[2]), tau_mpq_(constru[3]),
    tau_pqu_(constru[4]), tau_psu_(constru[5]), tau_ru_(constru[6]), tau_f_(constru[7])
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
        ret(i, j) = (1. / delta_t_) * mass_ij;
        ret(1 * n_shape_fct_ + i, 1 * n_shape_fct_ + j) = mass_ij;
        ret(2 * n_shape_fct_ + i, 2 * n_shape_fct_ + j) = mass_ij;
        ret(3 * n_shape_fct_ + i, 3 * n_shape_fct_ + j) = mass_ij;
        ret(4 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) = mass_ij;

        //bulk derivatives
        grad = integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
            i, j, hyper_edge.geometry);
        lSol_float_t mdx_ij = grad[0];
        lSol_float_t mdy_ij = grad[1];
        ret(i, 1 * n_shape_fct_ + j) = -mdx_ij;
        ret(i, 3 * n_shape_fct_ + j) = -mdy_ij;
        ret(1 * n_shape_fct_ + i, 2 * n_shape_fct_ + j) = mdx_ij;
        ret(2 * n_shape_fct_ + i, j) = mdx_ij;
        ret(3 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) = mdx_ij;
        ret(4 * n_shape_fct_ + i, j) = mdy_ij;

        //bulk triple products
        lSol_float_t hu = 0;
        for (unsigned int k = 0; k < j; ++k)
        {
          grad = integrator::template integrate_vol_nablaphiphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
              i, j, k, hyper_edge.geometry);
          lSol_float_t tv_shx_sh_sh_ijk = grad[0];
          hu += tv_shx_sh_sh_ijk * ca[k];
          ret(i, k) -= 2. * onehalf * tv_shx_sh_sh_ijk * ca[j];        }
        {
          grad = integrator::template integrate_vol_nablaphiphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
              i, j, j, hyper_edge.geometry);
          lSol_float_t tv_shx_sh_sh_ijj = grad[0];
          hu += tv_shx_sh_sh_ijj * ca[j];
        }
        ret(i, j) -= 2. * onehalf * hu;
      }
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        SmallVec<hyEdge_dim(), lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          lSol_float_t b_sh_sh_ij = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
            i, j, bdr, hyper_edge.geometry);
          lSol_float_t dj_flux_ux = -b_sh_sh_ij * normal[0] * normal[0];
          //flux part of f^
          ret(i, j) += -tau_f_ * dj_flux_ux;
          //main part of p^
          ret(i, n_shape_fct_ + j) += b_sh_sh_ij * normal[0];
          //r^
          ret(i, 3 * n_shape_fct_ + j) += b_sh_sh_ij * normal[1];
          ret(i, j) += tau_ru_ * -b_sh_sh_ij * normal[1] * normal[1];

          if (normal[1] * normal[1] < eps && normal[0] < 0) //left
          {
            //p^
            ret(i, j) += tau_ppu_ * dj_flux_ux;
            //q^
            ret(n_shape_fct_ + i, 2 * n_shape_fct_ + j) -= b_sh_sh_ij * normal[0];
            ret(n_shape_fct_ + i, j) -= tau_pqu_ * dj_flux_ux;
            //s^
            ret(3 * n_shape_fct_ + i, 4 * n_shape_fct_ + j) -= b_sh_sh_ij * normal[0];
            ret(3 * n_shape_fct_ + i, j) -= tau_psu_ * dj_flux_ux;
          }
          else if (normal[1] * normal[1] < eps && normal[0] > 0) //right
          {
            //p^
            ret(i, j) += tau_mpu_ * dj_flux_ux;
            ret(i, 2 * n_shape_fct_ + j) += tau_mpq_ * -b_sh_sh_ij * normal[0] * normal[0];
          }
        }
      }
    }

    //main part of f^
    for (unsigned int bdr = 0; bdr < 2; ++bdr)
    {
      SmallVec<hyEdge_dim(), lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
      std::array<lSol_float_t, n_shape_fct_ * n_shape_fct_> u_int, q_int;
      u_int.fill(0.);
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        //i > j > k
        for (unsigned int j = 0; j < i; ++j)
        {
          for (unsigned int k = 0; k < j; ++k)
          {
            lSol_float_t tb_sh_sh_sh_ijk = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, k, j, bdr, hyper_edge.geometry);
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
            u_int[n_shape_fct_ * i + i] += tb_sh_sh_sh_iii * ca[i] * normal[0];
          }
        }
      }
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          ret(i, j) += 2 * onehalf * u_int[n_shape_fct_ * i + j];
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
    for (unsigned int bdr = 0; bdr < 4; ++bdr)
    {
      SmallVec<hyEdge_dim(), lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
      if (normal[1] * normal[1] < eps && normal[0] > 0 && hyper_edge.node_descriptor[bdr] == 0) //rechts
      {
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            const lSol_float_t b_shsk_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            //first eq
            //f^
            grad[i] += normal[0] * normal[0] * (-tau_f_) * b_shsk_ij * lambda_dir[bdr][j];
            //p^
            grad[i] += normal[0] * normal[0] * tau_mpu_ * b_shsk_ij * lambda_dir[bdr][j];
            grad[i] += normal[0] * normal[0] * tau_mpq_ * b_shsk_ij * lambda_dir[bdr][n_shape_bdr_ + j];
            //r^ doesn't contribute here

            //second eq
            grad[n_shape_fct_ + i] -= normal[0] * b_shsk_ij * lambda_dir[bdr][n_shape_bdr_ + j];
            //third eq
            grad[2 * n_shape_fct_ + i] -= normal[0] * b_shsk_ij * lambda_dir[bdr][j];
            //fourth eq
            grad[3 * n_shape_fct_ + i] -= normal[0] * b_shsk_ij * lambda_dir[bdr][2 * n_shape_bdr_ + j];
            //fifth eq doesn't contribute here
          }
        }
      }
      else if (normal[1] * normal[1] < eps && normal[0] < 0 && hyper_edge.node_descriptor[bdr] == 0) //left
      {
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            const lSol_float_t b_shsk_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            //first eq
            //f^
            grad[i] += normal[0] * normal[0] * (-tau_f_) * b_shsk_ij * lambda_dir[bdr][j];
            //p^
            grad[i] += normal[0] * normal[0] * tau_ppu_ * b_shsk_ij * lambda_dir[bdr][j];
            //r^ doesn't contribute here

            //second eq
            grad[n_shape_fct_ + i] -= normal[0] * normal[0] * tau_pqu_ * b_shsk_ij * lambda_dir[bdr][j];
            //third eq
            grad[2 * n_shape_fct_ + i] -= normal[0] * b_shsk_ij * lambda_dir[bdr][j];
            //fourth eq
            grad[3 * n_shape_fct_ + i] -= normal[0] * normal[0] * tau_psu_ * b_shsk_ij * lambda_dir[bdr][j];
            //fifth eq doesn't contribute here
          }
        }
      }
      else if (normal[0] * normal[0] < eps && hyper_edge.node_descriptor[bdr] == 0) //horizontal
      {
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            const lSol_float_t b_shsk_ij = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            //first equation
            grad[i] += normal[1] * normal[1] * tau_ru_ * b_shsk_ij * lambda_dir[bdr][j];
            //second, third and fourth equation don't contribute here
            //fifth equation
            grad[4 * n_shape_fct_ + i] -= normal[1] * b_shsk_ij * lambda_dir[bdr][j];
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
                                                            const lSol_float_t time) const
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
        residual[i] += (1. / delta_t_) * mass_ij * (ca[j] - hyper_edge.data.u_old[j]);
        residual[n_shape_fct_ + i] += mass_ij * ca[n_shape_fct_ + j];
        residual[2 * n_shape_fct_ + i] += mass_ij * ca[2 * n_shape_fct_ + j];
        residual[3 * n_shape_fct_ + i] += mass_ij * ca[3 * n_shape_fct_ + j];
        residual[4 * n_shape_fct_ + i] += mass_ij * ca[4 * n_shape_fct_ + j];
        
        //integrals involving derivatives
        grad =  integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, 
             decltype(hyEdgeT::geometry)>(i, j, hyper_edge.geometry);
        const lSol_float_t mdx_ij = grad[0], mdy_ij = grad[1];
        residual[i] -= mdx_ij * ca[n_shape_fct_ + j];
        residual[i] -= mdy_ij * ca[3 * n_shape_fct_ + j];
        residual[n_shape_fct_ + i] += mdx_ij * ca[2 * n_shape_fct_ + j];
        residual[2 * n_shape_fct_ + i] += mdx_ij * ca[j];
        residual[3 * n_shape_fct_ + i] += mdx_ij * ca[4 * n_shape_fct_ + j];
        residual[4 * n_shape_fct_ + i] += mdy_ij * ca[j];

        //bulk triple products
        lSol_float_t u2_int = 0;
        for (unsigned int k = 0; k < j; ++k)
        {
          const lSol_float_t tv_ix_jk = integrator::template integrate_vol_nablaphiphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
              i, j, k, hyper_edge.geometry)[0];
          u2_int += 2 * tv_ix_jk * ca[j] * ca[k];
        }
        const lSol_float_t tv_ix_jj = integrator::template integrate_vol_nablaphiphiphi<SmallVec<hyEdge_dim(), lSol_float_t>, decltype(hyEdgeT::geometry)>(
            i, j, j, hyper_edge.geometry)[0];
        u2_int += tv_ix_jj * ca[j] * ca[j];
        residual[i] -= onehalf * u2_int;
      }

      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
      {
        //fixed integrals for notational convenience
        SmallVec<hyEdge_dim(), lSol_float_t> normal = hyper_edge.geometry.local_normal(bdr);
        lSol_float_t u_int = 0, p_int = 0, q_int = 0, r_int = 0, s_int = 0;
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
        {
          lSol_float_t h = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
              i, j, bdr, hyper_edge.geometry);
          u_int += h * ca[j];
          p_int += h * ca[1 * n_shape_fct_ + j];
          q_int += h * ca[2 * n_shape_fct_ + j];
          r_int += h * ca[3 * n_shape_fct_ + j];
          s_int += h * ca[4 * n_shape_fct_ + j];
        }
        lSol_float_t uh_int = 0, qh_int = 0, sh_int = 0;
        for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        {
          lSol_float_t h = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
              i, j, bdr, hyper_edge.geometry);
          uh_int += h * lambda_values[bdr][j];
          qh_int += h * lambda_values[bdr][n_shape_bdr_ + j];
          sh_int += h * lambda_values[bdr][2 * n_shape_bdr_ + j];
        }
        lSol_float_t flux_ux = (uh_int - u_int) * normal[0];
        lSol_float_t flux_uy = (uh_int - u_int) * normal[1];

        lSol_float_t trace_f = - tau_f_ * flux_ux;
        lSol_float_t trace_p = 0., trace_q = 0., trace_r = 0., trace_s = 0.;
        if (normal[1] * normal[1] < eps && normal[0] < 0) //left
        {
          trace_p = p_int + tau_ppu_ * flux_ux;
          trace_q = q_int + tau_pqu_ * flux_ux;
          trace_s = s_int + tau_psu_ * flux_ux;
        }
        else if (normal[1] * normal[1] < eps && normal[0] > 0) //right
        {
          trace_p = p_int + tau_mpu_ * flux_ux + tau_mpq_ * (qh_int - q_int) * normal[0];
          trace_q = qh_int;
          trace_s = sh_int;
        }
        else if (normal[0] * normal[0] < eps) //horizontal
        {
          trace_r = r_int + tau_ru_ * flux_uy;
        }

        residual[i] += (trace_f + trace_p) * normal[0];
        residual[i] += trace_r * normal[1];
        residual[n_shape_fct_ + i] -= trace_q * normal[0];
        residual[2 * n_shape_fct_ + i] -= uh_int * normal[0];
        residual[3 * n_shape_fct_ + i] -= trace_s * normal[0];
        residual[4 * n_shape_fct_ + i] -= uh_int * normal[1];
      }

      /*
      lSol_float_t rhs = integrator::template integrate_vol_phifunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side, Point<hyEdge_dimT, lSol_float_t> > (i, hyper_edge.geometry, time);
      residual[0 * n_shape_fct_ + i] += rhs;
      */
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
            residual[i] += onehalf * 2 * sss_ijk * ca[j] * ca[k] * nx;
            residual[j] += onehalf * 2 * sss_ijk * ca[i] * ca[k] * nx;
            residual[k] += onehalf * 2 * sss_ijk * ca[j] * ca[i] * nx;
          }
          const lSol_float_t sss_ijj = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
                i, j, j, bdr, hyper_edge.geometry);
          residual[i] += onehalf * sss_ijj * ca[j] * ca[j] * nx;
          residual[j] += onehalf * 2 * sss_ijj * ca[i] * ca[j] * nx;
        }
        for (unsigned int k = 0; k < i; ++k)
        {
          const lSol_float_t sss_iik = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
                i, i, k, bdr, hyper_edge.geometry);
          residual[k] += onehalf * sss_iik * ca[i] * ca[i] * nx;
          residual[i] += onehalf * 2 * sss_iik * ca[i] * ca[k] * nx;
        }
        const lSol_float_t sss_iii = integrator::template integrate_bdr_phiphiphi<decltype(hyEdgeT::geometry)>(
              i, i, i, bdr, hyper_edge.geometry);
        residual[i] += onehalf * sss_iii * ca[i] * ca[i] * nx;
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
      //delete
      //out[bdr].fill(0.);

      if (normal[1] * normal[1] < eps && normal[0] > 0 && hyper_edge.node_descriptor[bdr] == 0) //right
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            lSol_float_t b_int = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                j, i, bdr, hyper_edge.geometry);
            out[bdr][i] += normal[0] * normal[0] * (tau_mpu_ - tau_f_) * b_int * lambda_dir[bdr][j];
            out[bdr][i] += normal[0] * normal[0] * tau_mpq_ * b_int * lambda_dir[bdr][n_shape_bdr_ + j];
            out[bdr][n_shape_bdr_ + i] += normal[0] * b_int * lambda_dir[bdr][n_shape_bdr_ + j];
            out[bdr][2 * n_shape_bdr_ + i] += normal[0] * b_int * lambda_dir[bdr][2 * n_shape_bdr_ + j];
          }
        }
      }
      else if (normal[1] * normal[1] < eps && normal[0] < 0 && hyper_edge.node_descriptor[bdr] == 0) //left
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            lSol_float_t b_int = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                j, i, bdr, hyper_edge.geometry);
            out[bdr][i] += normal[0] * normal[0] * (tau_ppu_ - tau_f_ ) * b_int * lambda_dir[bdr][j];
            out[bdr][n_shape_bdr_ + i] += normal[0] * normal[0] * tau_pqu_ * b_int * lambda_dir[bdr][j];
            out[bdr][n_shape_bdr_ + i] += normal[0] * normal[0] * tau_psu_ * b_int * lambda_dir[bdr][j];
          }
        }
      }
      else if (normal[0] * normal[0] < eps && hyper_edge.node_descriptor[bdr] == 0) //horizontal
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          for (unsigned int j = 0; j < n_shape_bdr_; ++j)
          {
            lSol_float_t b_int = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                j, i, bdr, hyper_edge.geometry);
            out[bdr][i] += normal[1] * normal[1] * tau_ru_ * b_int * lambda_dir[bdr][j];
          }
        }
      }
    }
    return out;
  }

  /*
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
      if (normal[1] * normal[1] < eps && normal[0] > 0 && hyper_edge.node_descriptor[bdr] == 0) //V-
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
            out[bdr][i] -= normal[0] * normal[0] * b_sk_sh_ij * (tau_mzpu_ - tau_f_) * coeff_dir[j];
            out[bdr][i] -= normal[0] * normal[0] * b_sk_sh_ij * (tau_mzpv_) * coeff_dir[4 * n_shape_fct_ + j];
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
      else if (normal[1] * normal[1] < eps && normal[0] < 0 && hyper_edge.node_descriptor[bdr] == 0) //V+
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
            out[bdr][i] -= normal[0] * normal[0] * b_sk_sh_ij * (tau_pzpu_ - tau_f_) * coeff_dir[j];
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
      else if (normal[0] * normal[0] < eps && normal[1] < 0 && hyper_edge.node_descriptor[bdr] == 0) //H+
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
      else if (normal[0] * normal[0] < eps && normal[1] > 0 && hyper_edge.node_descriptor[bdr] == 0) //H-
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
*/
  /*
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
*/
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
 /*
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
*/

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
      if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] > 0 && hyper_edge.node_descriptor[bdr] == 0) //right
      {
        const lSol_float_t normal_comp = loc_normal[bdr][0];
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
	      {
          //required integrals
  	      lSol_float_t uh_int = 0, qh_int = 0, sh_int = 0;
  	      lSol_float_t u_int = 0, p_int = 0, q_int = 0;
    	    for (unsigned int j = 0; j < n_shape_bdr_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_psipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
	          uh_int += c * lambda_values_in[bdr][j];
      	    qh_int += c * lambda_values_in[bdr][n_shape_bdr_ + j];
      	    sh_int += c * lambda_values_in[bdr][2 * n_shape_bdr_ + j];
      	  }
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j) 
          {
            lSol_float_t c = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                           j, i, bdr, hyper_edge.geometry);
      	    u_int += c * coeff[j];
	          p_int += c * coeff[n_shape_fct_ + j];
      	    q_int += c * coeff[2 * n_shape_fct_ + j];
      	  }
          lSol_float_t trace_f = 0, trace_p = 0, trace_q = 0, trace_s = 0;
          //f^
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
      	  {
            for (unsigned int k = 0; k < n_shape_fct_; ++k) 
	          {
              lSol_float_t c = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
      	      trace_f += onehalf * c * coeff[j] * coeff[k];
	          }
	        }
          trace_f -= tau_f_ * (uh_int - u_int) * normal_comp;
          //p^
          trace_p = p_int + normal_comp * (tau_mpu_ * (uh_int - u_int) + tau_mpq_ * (qh_int - q_int));
          //q^
          trace_q = qh_int;
          //s^
          trace_s = sh_int;
          //output
          lambda_values_out[bdr][i] = (trace_f + trace_p) * normal_comp;
          lambda_values_out[bdr][n_shape_bdr_ + i] = trace_q * normal_comp;
          lambda_values_out[bdr][2 * n_shape_bdr_ + i] = trace_s * normal_comp;
      	}
      }
      if (loc_normal[bdr][1] * loc_normal[bdr][1] < eps && loc_normal[bdr][0] < 0 && hyper_edge.node_descriptor[bdr] == 0) //left
      {
        const lSol_float_t normal_comp = loc_normal[bdr][0];
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
	      {
          //required integrals
  	      lSol_float_t uh_int = 0;
  	      lSol_float_t u_int = 0, p_int = 0, q_int = 0, s_int = 0;
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
	          p_int += c * coeff[n_shape_fct_ + j];
      	    q_int += c * coeff[2 * n_shape_fct_ + j];
      	    s_int += c * coeff[4 * n_shape_fct_ + j];
      	  }
          lSol_float_t trace_f = 0, trace_p = 0, trace_q = 0, trace_s = 0;
          //f^
      	  for (unsigned int j = 0; j < n_shape_fct_; ++j)
      	  {
            for (unsigned int k = 0; k < n_shape_fct_; ++k) 
	          {
              lSol_float_t c = integrator::template integrate_bdr_phiphipsi<decltype(hyEdgeT::geometry)>(
			      j, k, i, bdr, hyper_edge.geometry);
      	      trace_f += onehalf * c * coeff[j] * coeff[k];
	          }
	        }
          trace_f -= tau_f_ * (uh_int - u_int) * normal_comp;
          //p^
          trace_p = p_int + normal_comp * tau_ppu_ * (uh_int - u_int);
          //q^
          trace_q = q_int + tau_pqu_ * (uh_int - u_int) * normal_comp;
          //s^
          trace_s = s_int + tau_psu_ * (uh_int - u_int) * normal_comp;
          //output
          lambda_values_out[bdr][i] = (trace_f + trace_p) * normal_comp;
          lambda_values_out[bdr][n_shape_bdr_ + i] = trace_q * normal_comp;
          lambda_values_out[bdr][2 * n_shape_bdr_ + i] = trace_s * normal_comp;
      	}
      }
      if (loc_normal[bdr][0] * loc_normal[bdr][0] < eps && hyper_edge.node_descriptor[bdr] == 0)
      {
        const lSol_float_t normal_comp = loc_normal[bdr][1];
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
	      {
          //required integrals
  	      lSol_float_t uh_int = 0;
  	      lSol_float_t u_int = 0, r_int = 0;
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
      	    r_int += c * coeff[3 * n_shape_fct_ + j];
          }
          //r^
          lSol_float_t trace_r = r_int + tau_ru_ * (uh_int - u_int) * normal_comp;
          //output
          lambda_values_out[bdr][i] = trace_r * normal_comp;
          lambda_values_out[bdr][n_shape_bdr_ + i] = 0.;
          lambda_values_out[bdr][2 * n_shape_bdr_ + i] = 0.;
        }
      }
      if (hyper_edge.node_descriptor[bdr] != 0)
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
            parameters::uh_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);

        } else 
        {
          lambda_values[bdr][i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::initial, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
        }
        //hyper_edge.data.uh_old[bdr][i] = lambda_values[bdr][i];	-> unecessary bc. of set_skeleton_data
        if (is_right<parameters>(hyper_edge.node_descriptor[bdr]))
        {
          lambda_values[bdr][n_shape_bdr_ + i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::qh_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
          lambda_values[bdr][2 * n_shape_bdr_ + i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::sh_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);

        } else {
          lambda_values[bdr][n_shape_bdr_ + i] = 0.;
          lambda_values[bdr][2 * n_shape_bdr_ + i] = 0.;
        }
      }
    }
    return lambda_values;
  }

  //TODO: delete asap
  template <typename hyEdgeT, typename SmallMatT>
  static inline void set_skeleton_data(const SmallMatT& lambda_values,
      hyEdgeT& hyper_edge)
  {
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dim(); ++bdr)
    {
      for (unsigned int i = 0; i < n_shape_bdr_; ++i)
      {
      }
    }
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline void set_data(const SmallMatT& lambda_values, 
      hyEdgeT& hyper_edge, const lSol_float_t time) const
  {
    SmallVec<n_loc_dofs_, lSol_float_t> coeff;
    newton(lambda_values, coeff, hyper_edge, time);
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      hyper_edge.data.u_old[i] = coeff[i];
    }
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
            parameters::uh_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
        }
      }
      if (is_right<parameters>(hyper_edge.node_descriptor[bdr]))
      {
        for (unsigned int i = 0; i < n_shape_bdr_; ++i)
        {
          lambda_values_out[bdr][n_shape_bdr_ + i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::qh_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
          lambda_values_out[bdr][2 * n_shape_bdr_ + i] = integrator::template integrate_bdrUni_psifunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::sh_value, Point<hyEdge_dimT, lSol_float_t> > (i, bdr, hyper_edge.geometry, time);
        }
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
  /*
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
             ZK<hyEdge_dimT, poly_deg, quad_deg, parametersT, lSol_float_t>::system_dim>
      point_vals;

    for (unsigned int d = 0; d < system_dim-1; ++d)
    {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[i + d * n_shape_fct_];
      for (unsigned int pt = 0; pt < Hypercube<hyEdge_dimT>::pow(abscissas_sizeT); ++pt)
        point_vals[d][pt] = integrator::shape_fun_t::template lin_comb_fct_val<float>(
          coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT> >(pt, helper));
    }
    for (unsigned int d = system_dim-1; d < system_dim; ++d)
    {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[i + (d + 2) * n_shape_fct_];
      for (unsigned int pt = 0; pt < Hypercube<hyEdge_dimT>::pow(abscissas_sizeT); ++pt)
        point_vals[d][pt] = integrator::shape_fun_t::template lin_comb_fct_val<float>(
          coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT> >(pt, helper));
    }

    return point_vals;
  }
  */


private:

};

}
