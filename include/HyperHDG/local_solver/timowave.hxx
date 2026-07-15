#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <iostream>
#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

#include <tuple>

template <unsigned int space_dimT, typename param_float_t = double>
struct TestTimoWave0
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 10U> dirichlet_nodes{1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Neumann boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return 0;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    return analytic_result_u(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    return analytic_result_phi(point, normal, time);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t time = 0.)
  {
    auto res = initial_u(point, time);
    return scalar_product(res, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t time = 0.)
  {
    auto res = initial_r(point, time);
    return scalar_product(res, normal);
  }

  static SmallVec<space_dimT, param_float_t> initial_u(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(1.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_v(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_s(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_r(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_n(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }

  static SmallVec<space_dimT, param_float_t> initial_m(const Point<space_dimT, param_float_t>& point, const param_float_t time = 0.) {
    SmallVec<space_dimT, param_float_t> res(0.);
    return res;
  }
};  // end of struct TestTimoWave0

namespace LocalSolver
{

/*!*************************************************************************************************
 * \brief   Local solver for the Timoshenko beam-network wave equation (dual-mixed HDG).
 *
 * Second-order-in-time elastic beam system, reduced to first order via the velocity/momentum
 * field and advanced by fully implicit Gauss collocation (hdg_gauss.pdf): one mass-shifted local
 * solve per stage. n_stages = 1 is the implicit midpoint rule, which reproduces the classic
 * Crank-Nicolson trajectory (stage values = CN midpoint averages; see HORK_GAUSS_PLAN.md).
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT = TestTimoWave0,
          typename lSol_float_t = double,
          unsigned int n_stages = 1>
class TimoshenkoWave
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of node elements, especially with respect to nodal shape functions.
   ************************************************************************************************/
  struct node_element
  {
    typedef std::tuple<TPP::ShapeFunction<
      TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT - 1>>>
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
   * \brief   Return template parameter \c hyEdge_dimT.
   *
   * \retval  hyEdge_dimT    Dimension of hypergraph's hyperedges.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Evaluate amount of global degrees of freedom per hypernode.
   *
   * This number must be equal to HyperNodeFactory::n_glob_dofs_per_node()() of the HyperNodeFactory
   * cooperating with this object.
   *
   * \retval  n_dofs        Number of global degrees of freedom per hypernode.
   ************************************************************************************************/
  static constexpr unsigned int n_glob_dofs_per_node()
  {
    return 2 * space_dim * Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return 6*space_dim; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return 6*space_dim; }

 private:
  // -----------------------------------------------------------------------------------------------
  // Private, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Number of local shape functions (with respect to all spatial dimensions).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_fct_ = Hypercube<hyEdge_dimT>::pow(poly_deg + 1);
  /*!***********************************************************************************************
   * \brief   Number oflocal  shape functions (with respect to a face / hypernode).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_bdr_ = Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  /*!***********************************************************************************************
   * \brief   Number of (local) degrees of freedom per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_dofs_ = 6 * space_dim * n_shape_fct_;
  /*!***********************************************************************************************
   * \brief   Block-diagonal decomposition of the local matrix into independent components.
   *
   * The local matrix decouples into 4 connected components (see assemble_loc_matrix). Each
   * variable group dim in 0..2*space_dim-1 forms a self-coupled triplet of variable-blocks
   * {sigma, w, wdot} = {dim, 2*space_dim+dim, 4*space_dim+dim}; the only links between triplets
   * are the t x n cross-product terms, which (for space_dim==3) couple the transverse shear
   * force with the perpendicular bending rotation. The components are, in units of variable-
   * blocks of width n_shape_fct_:
   *   A = {n_x, u_x, v_x}            (axial)        size 3
   *   B = {m_x, r_x, s_x}            (torsion)      size 3
   *   C = {n_z,u_z,v_z, m_y,r_y,s_y} (shear+bend)   size 6
   *   D = {n_y,u_y,v_y, m_z,r_z,s_z} (shear+bend)   size 6
   * Level 2 (assemble_schur / solve_local_problem) further eliminates each triplet's diagonal
   * mass blocks (sigma, wdot) explicitly, leaving only the displacement w: A,B reduce to an
   * nw_ x nw_ system, C,D to an n2w_ x n2w_ system.
   ************************************************************************************************/
  static_assert(space_dim == 3, "Block decomposition of the local matrix assumes space_dim==3 "
                                "(the t x n cross-product coupling is 3D-specific).");
  static constexpr unsigned int nw_ = n_shape_fct_;       // single displacement block
  static constexpr unsigned int n2w_ = 2 * n_shape_fct_;  // cross-coupled displacement pair
  // Variable-block index (in units of n_shape_fct_) of the (sigma, w, wdot) of triplet `dim`:
  static constexpr unsigned int sig_blk(unsigned int dim) { return dim; }
  static constexpr unsigned int w_blk(unsigned int dim) { return 2 * space_dim + dim; }
  static constexpr unsigned int wdot_blk(unsigned int dim) { return 4 * space_dim + dim; }
  // Triplet dims forming the four independent components: A,B single; C,D couple a force group
  // dim_n with a rotation group dim_r via the cross product (sign +1 for C, -1 for D).
  static constexpr unsigned int compA_dim = 0;
  static constexpr unsigned int compB_dim = 3;
  static constexpr unsigned int compC_dim_n = 2, compC_dim_r = 4;
  static constexpr unsigned int compD_dim_n = 1, compD_dim_r = 5;
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   *
   * This allows to the use of this quantity as template parameter in member functions.
   ************************************************************************************************/
  static constexpr unsigned int system_dim = system_dimension();

  /*!***********************************************************************************************
   * \brief   Step parameters of the local solve (see also hdg_gauss.pdf and HORK_GAUSS_PLAN.md).
   *
   * The only runtime time-discretisation inputs are:
   *   tau_      HDG stabilisation parameter (globally constant), tau > 0.
   *   delta_t_  time step Delta t.
   * The local system is assembled in the STAGE NORMALIZATION of hdg_gauss.pdf (eq 5): the step
   * enters the matrix only through sigma_l = 1/(theta_l*Delta t), which weights the (y,z)/(z,y)
   * couplings (see assemble_loc_matrix / assemble_schur). Deviation from the note: the z-rows keep
   * their C_u scaling -- (z,z) = M, (z,y) = -sigma*C_u*M -- instead of the note's C_u^{-1} d-form,
   * so massless welds (C_u = 0) stay regular (z = 0 there). This costs the literal symmetry of the
   * local saddle (row scaling), but leaves the condensed trace operator unchanged.
   ************************************************************************************************/
  const lSol_float_t tau_;
  const lSol_float_t delta_t_;
  /*!***********************************************************************************************
   * \brief   Gauss tableau data, one entry per stage (hdg_gauss.pdf; HORK_GAUSS_PLAN.md).
   *
   * Per stage l: Butcher eigenvalue theta_l (stage shift sigma_l = 1/(theta_l*dt)), collocation
   * node c_l, load weight omega_l = (T^{-1} 1)_l, recombination weight w_l = (d^T T)_l with
   * d^T = b^T A^{-1}; stage_affine_ = 1 - sum_j d_j multiplies the old state in the endpoint
   * update y+ = stage_affine_*y^n + sum_l w_l*y_l. Currently only the single-stage tableau
   * (implicit midpoint: theta = c = 1/2, omega = 1, w = 2, affine = -1), whose trajectory is
   * Crank-Nicolson's. s >= 2 has complex conjugate theta_l and lands with the complex LAPACK
   * plumbing (eigen-decomposition of the Butcher matrix via dgeev in the constructor).
   ************************************************************************************************/
  static_assert(n_stages == 1, "only the s = 1 (implicit midpoint) tableau is implemented; "
                               "s >= 2 requires the complex stage plumbing");
  static constexpr std::array<lSol_float_t, n_stages> stage_theta_{0.5};
  static constexpr std::array<lSol_float_t, n_stages> stage_c_{0.5};
  static constexpr std::array<lSol_float_t, n_stages> stage_omega_{1.};
  static constexpr std::array<lSol_float_t, n_stages> stage_w_{2.};
  static constexpr lSol_float_t stage_affine_ = -1.;
  /*!***********************************************************************************************
   * \brief   Mass-shift weight sigma_l = 1/(theta_l * Delta t) of the stage normalization.
   ************************************************************************************************/
  lSol_float_t sigma(const unsigned int stage = 0) const
  {
    return 1. / (stage_theta_[stage] * delta_t_);
  }
  /*!***********************************************************************************************
   * \brief   Solve local problems via a full-matrix LU instead of the displacement Schur path.
   *
   * The two solve routines are two views of the same local saddle system:
   *   assemble_loc_matrix   the full n_loc_dofs_ saddle form -- the authoritative definition;
   *   assemble_schur        its hand-eliminated displacement Schur complement (production path);
   * they must agree. -loc_lu_full selects the reference LU, kept for deep-dt convergence studies
   * where the Schur mass shift loses stiffness digits (see data_type below). The full-matrix path
   * is also the one that takes a complex sigma for the Gauss stages.
   ************************************************************************************************/
  const bool full_lu_ = false;

  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT>>,
    lSol_float_t>
    integrator;

  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
 public:
  /*!***********************************************************************************************
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef std::vector<double> constructor_value_type;

  struct data_type
  {
    SmallVec<space_dim*n_shape_fct_, lSol_float_t> u_old, v_old, r_old, s_old, n_old, m_old;
    // Cached Schur factorization of the local matrix (see assemble_schur). The matrix only depends
    // on geometry and (tau_, delta_t_, tableau), all constant across time steps / Krylov iterations.
    // Shared per-edge integral blocks (row-major n x n, M diagonal) used in reduction/back-sub:
    std::array<lSol_float_t, n_shape_fct_> M_inv;                       // 1 / diag(M)
    std::array<lSol_float_t, n_shape_fct_ * n_shape_fct_> Gmat;         // G   = int (grad phi) phi
    std::array<lSol_float_t, n_shape_fct_ * n_shape_fct_> BmG;          // B-G (boundary normal - G)
    SmallVec<4 * space_dim, lSol_float_t> extra;                        // C_n, C_m, C_u, C_r
    // LU factors (column-major) of the displacement Schur systems of the four components:
    std::array<lSol_float_t, nw_ * nw_> lu_A, lu_B;
    std::array<int, nw_> ipiv_A, ipiv_B;
    std::array<lSol_float_t, n2w_ * n2w_> lu_C, lu_D;
    std::array<int, n2w_> ipiv_C, ipiv_D;
    bool loc_mat_factorized = false;
    // Cached LU of the FULL local matrix (full_lu_ mode). Forming the displacement Schur
    // complement explicitly sums (C_u/(theta*dt^2))*M with the O(1) stiffness terms, so for very
    // small dt the stiffness drowns in round-off (e_trace floors ~1e-7 at dt~2e-5). The pivoted
    // full-matrix LU keeps the scales in separate entries and reaches ~1e-9; use it (option
    // -loc_lu_full) for deep convergence studies, the Schur path for production runs.
    SmallSquareMat<n_loc_dofs_, lSol_float_t> full_lu;
    std::array<int, n_loc_dofs_> full_ipiv;
    bool full_lu_factorized = false;
    // Per-stage local solution (q,y,z) stashed by set_data(zeta), consumed by finalize_step.
    // Real at s = 1; becomes complex per-representative with the multi-stage plumbing.
    std::array<SmallVec<n_loc_dofs_, lSol_float_t>, n_stages> stage_coeffs;
  };
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   tau           Penalty parameter of HDG scheme.
   ************************************************************************************************/
  // NOTE: tau, delta_t, [full_lu]
  TimoshenkoWave(const constructor_value_type& vals = std::vector(2, 1.)) : tau_(vals[0]),
    delta_t_(vals[1]), full_lu_(vals.size() > 2 && vals[2] != 0.) {}

  template <typename point_t, typename geom_t,
            lSol_float_t fun(const point_t&, const point_t&, const lSol_float_t),
            unsigned int n_comps = 3>
  std::array<lSol_float_t, n_comps> integrate_vol_phivecfunccomp_beam(const unsigned int i,
                                                                         std::array<int, n_comps> comps, geom_t& geom, const lSol_float_t time) const
  {
    std::array<lSol_float_t, n_comps> ret;
    for (unsigned int j = 0; j < n_comps; j++) {
      ret[j] = integrator::template integrate_vol_phivecfunccomp<
        point_t, geom_t, fun, Point<hyEdge_dimT, lSol_float_t>>(i, comps[j], geom,
                                                                         time);
    }
    return ret;
  }


  template <typename point_t,
            typename geom_t,
            lSol_float_t fun(const point_t&, const point_t&, const lSol_float_t),
            unsigned int n_comps = 3>
  std::array<lSol_float_t, n_comps> integrate_bdr_phivecfunccomp_beam(const unsigned int i,
                                               const unsigned int bdr,
                                               std::array<int, n_comps> comps,
                                               geom_t& geom,
                                               const lSol_float_t time = 0.) const
  {
   std::array<lSol_float_t, n_comps> ret;
   for (unsigned int j = 0; j < n_comps; j++) {
     ret[j] = integrator::template integrate_bdr_phivecfunccomp<
       point_t, geom_t, fun, Point<hyEdge_dimT, lSol_float_t>>(
                 i, bdr, comps[j], geom, time);
   }

    return ret;
  }



  /*!***********************************************************************************************
   * \brief   The local saddle system in stage normalization (hdg_gauss.pdf eq 5), one assembler.
   *
   * Blocks (per variable group dim, in units of n_shape_fct_; see the block decomposition above):
   *   (q,q) = M/C_sig,  (q,y) = -G,  (y,q) = B-G,  (y,y) = tau*F        [always]
   *   (y,z) = +sigma*M, (z,y) = -sigma*C_u*M, (z,z) = M                 [with_z only]
   * plus the t x n cross-product couplings. with_z = false yields the STATIC system (the sigma
   * couplings and z rows dropped) used by the initializers; sigma is ignored then. The z-rows keep
   * their C_u scaling (deviation from the note's C_u^{-1} d-form) so massless welds stay regular:
   * C_u = 0 gives z = 0 and a vanishing mass shift, not a division by zero.
   ************************************************************************************************/
  template <bool with_z, typename hyEdgeT>
  inline SmallSquareMat<(with_z ? 6 : 4) * space_dim * n_shape_fct_, lSol_float_t>
  assemble_loc_matrix(hyEdgeT& hyper_edge, const lSol_float_t sigma = 0.) const
  {
    SmallSquareMat<(with_z ? 6 : 4) * space_dim * n_shape_fct_, lSol_float_t> local_mat;
    SmallVec<4 * space_dim, lSol_float_t> extra_coeffs = get_extra_coeffs(hyper_edge);

    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        const lSol_float_t vol_integral =
          integrator::template integrate_vol_phiphi(i, j, hyper_edge.geometry);
        const auto grad_int_vec =
          integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dimT, lSol_float_t>,
                                                         decltype(hyEdgeT::geometry)>(
            i, j, hyper_edge.geometry);
        lSol_float_t face_integral = 0.;
        SmallVec<hyEdge_dimT, lSol_float_t> normal_int_vec(0.);
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        {
          const lSol_float_t helper =
            integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
              i, j, face, hyper_edge.geometry);
          face_integral += helper;
          normal_int_vec += helper * hyper_edge.geometry.local_normal(face);
        }

        for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
        {
          // compliance a(q,q) and equilibrium -b(q,w_y)
          local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
            vol_integral / extra_coeffs[dim];
          local_mat(dim * n_shape_fct_ + i, (2 * space_dim + dim) * n_shape_fct_ + j) -=
            grad_int_vec[0];
          // b^T and tau<y,w_y>
          local_mat((2 * space_dim + dim) * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
            (normal_int_vec[0] - grad_int_vec[0]);
          local_mat((2 * space_dim + dim) * n_shape_fct_ + i,
                    (2 * space_dim + dim) * n_shape_fct_ + j) += tau_ * face_integral;
          if constexpr (with_z)
          {
            // inertia (z,z) = M and the sigma-weighted couplings; z-rows carry the C_u scaling
            local_mat((4 * space_dim + dim) * n_shape_fct_ + i,
                      (4 * space_dim + dim) * n_shape_fct_ + j) += vol_integral;
            local_mat((4 * space_dim + dim) * n_shape_fct_ + i,
                      (2 * space_dim + dim) * n_shape_fct_ + j) -=
              sigma * extra_coeffs[2 * space_dim + dim] * vol_integral;
            local_mat((2 * space_dim + dim) * n_shape_fct_ + i,
                      (4 * space_dim + dim) * n_shape_fct_ + j) += sigma * vol_integral;
          }
        }

        // cross product (textbook orientation: + i x r and + i x n)
        local_mat(2 * n_shape_fct_ + i, (3 * space_dim + 1) * n_shape_fct_ + j) += vol_integral;
        local_mat(1 * n_shape_fct_ + i, (3 * space_dim + 2) * n_shape_fct_ + j) -= vol_integral;
        local_mat((3 * space_dim + 2) * n_shape_fct_ + i, 1 * n_shape_fct_ + j) += vol_integral;
        local_mat((3 * space_dim + 1) * n_shape_fct_ + i, 2 * n_shape_fct_ + j) -= vol_integral;
      }

    return local_mat;
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_lambda(
    const SmallMatT& lambda_values,
    hyEdgeT& hyper_edge) const;

  /*!***********************************************************************************************
   * \brief   Gauss stage loads (hdg_gauss.pdf eq 7): history enters through (y^n, z^n) only.
   *
   * g_y = (f(t_stage), w_y) + sigma*omega*(z^n, w_y);  g_z = -sigma*omega*C_u*(y^n, w_z) in our
   * C_u-scaled z-row normalization. Data (f, Dirichlet values) are sampled at the collocation node
   * t_stage = time - dt + c_l*dt, not endpoint-averaged: at s = 1 this is the O(dt^2) load
   * difference to the retired trapezoid scheme (exact match for f = 0 and time-independent BC).
   * No flux caches, no old trace: q and lambda are algebraic and never history.
   ************************************************************************************************/
  template <typename hyEdgeT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_stage(hyEdgeT& hyper_edge,
                                                                const lSol_float_t time,
                                                                const unsigned int stage = 0) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    SmallVec<n_loc_dofs_, lSol_float_t> rhs;
    const lSol_float_t t_stage = time - (1. - stage_c_[stage]) * delta_t_;
    const lSol_float_t sig_om = sigma(stage) * stage_omega_[stage];

    // body loads act on material only if the parameters opt in (massless welds unloaded)
    bool loaded = true;
    if constexpr (requires { parameters::massless_unloaded; })
      if (parameters::massless_unloaded && hyper_edge.geometry.has_extra_data())
        loaded = hyper_edge.geometry.extra_data()[0] > 0;

    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      if (loaded)
      {
        auto integrals = integrate_vol_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::right_hand_side_n>(i, {1, -1, -2},
                                                                      hyper_edge.geometry, t_stage);
        for (unsigned int comp = 0; comp < 3; ++comp)
          rhs[(2 * space_dim + comp) * n_shape_fct_ + i] += integrals[comp];
        integrals = integrate_vol_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::right_hand_side_m>(i, {1, -1, -2},
                                                                      hyper_edge.geometry, t_stage);
        for (unsigned int comp = 0; comp < 3; ++comp)
          rhs[(3 * space_dim + comp) * n_shape_fct_ + i] += integrals[comp];
      }

      // Dirichlet data at the collocation node (bit-6 static-only faces keep their trace)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        if (hyper_edge.node_descriptor[face] & (1 << 6)) continue;
        if (!hyper_edge.node_descriptor[face]) continue;
        auto integrals = integrate_bdr_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u>(i, face, {1, -1, -2},
                                                                      hyper_edge.geometry, t_stage);
        for (unsigned int comp = 0; comp < 3; ++comp)
        {
          rhs[(0 * space_dim + comp) * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * integrals[comp];
          rhs[(2 * space_dim + comp) * n_shape_fct_ + i] += tau_ * integrals[comp];
        }
        integrals = integrate_bdr_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi>(i, face, {1, -1, -2},
                                                                        hyper_edge.geometry,
                                                                        t_stage);
        for (unsigned int comp = 0; comp < 3; ++comp)
        {
          rhs[(1 * space_dim + comp) * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * integrals[comp];
          rhs[(3 * space_dim + comp) * n_shape_fct_ + i] += tau_ * integrals[comp];
        }
      }
    }

    // history: sigma*omega mass pairings of (y^n, z^n); M = area * Id (orthonormal basis)
    const auto extra = get_extra_coeffs(hyper_edge);
    const lSol_float_t mass = hyper_edge.geometry.area();
    for (unsigned int d = 0; d < space_dim; ++d)
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        rhs[(2 * space_dim + d) * n_shape_fct_ + i] +=
          sig_om * mass * hyper_edge.data.v_old[d * n_shape_fct_ + i];
        rhs[(3 * space_dim + d) * n_shape_fct_ + i] +=
          sig_om * mass * hyper_edge.data.s_old[d * n_shape_fct_ + i];
        rhs[(4 * space_dim + d) * n_shape_fct_ + i] -=
          sig_om * mass * extra[2 * space_dim + d] * hyper_edge.data.u_old[d * n_shape_fct_ + i];
        rhs[(5 * space_dim + d) * n_shape_fct_ + i] -=
          sig_om * mass * extra[3 * space_dim + d] * hyper_edge.data.r_old[d * n_shape_fct_ + i];
      }

    return rhs;
  }

  template <class hyEdgeT, typename SmallMatT>
  inline SmallMatT glob_dof_to_loc_dof(
    const SmallMatT& glob_dof,
    hyEdgeT& hyper_edge) const
  {
    SmallMatT loc_dof;

    Point<space_dim, lSol_float_t> normal_vec =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);

    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        loc_dof[i] += normal_vec[dim] * glob_dof[i+dim*n_shape_fct_];
      }

    for (unsigned int ind = 0; ind < space_dim - 1; ++ind)
    {
      normal_vec = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(ind);

      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          loc_dof[i+(1+ind)*n_shape_fct_] += normal_vec[dim] * glob_dof[i+n_shape_fct_*dim];
        }
    }

    return loc_dof;
  }

  template <class hyEdgeT, typename SmallMatT>
  inline SmallMatT loc_dof_to_glob_dof(
    const SmallMatT& loc_dof,
    hyEdgeT& hyper_edge) const
  {
    SmallMatT glob_dof;
    Point<space_dim, lSol_float_t> normal_vec =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        glob_dof[i + dim * n_shape_fct_] += normal_vec[dim] * loc_dof[i];

    for (unsigned int ind = 0; ind < space_dim - 1; ++ind)
    {
      normal_vec = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(ind);
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          glob_dof[i + dim * n_shape_fct_] += normal_vec[dim] * loc_dof[i + (1 + ind) * n_shape_fct_];
    }
    return glob_dof;
  }

  template <class hyEdgeT, typename SmallMatT>
  inline std::array<std::array<double, 2 * space_dim>, 2 * hyEdge_dimT> node_dof_to_edge_dof(
    const SmallMatT& glob_lambda,
    hyEdgeT& hyper_edge) const
  {
    std::array<std::array<double, 2 * space_dim>, 2 * hyEdge_dimT> loc_lambda;
    hy_assert(loc_lambda.size() == 2, "Only implemented in one dimension!");
    for (unsigned int i = 0; i < loc_lambda.size(); ++i)
    {
      hy_assert(loc_lambda[i].size() == 6, "Only implemented in one dimension!");
      loc_lambda[i].fill(0.);
    }

    Point<space_dim, lSol_float_t> normal_vec =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        loc_lambda[i][0] += normal_vec[dim] * glob_lambda[i][dim];
        loc_lambda[i][0 + space_dim] += normal_vec[dim] * glob_lambda[i][dim + space_dim];
      }

    for (unsigned int ind = 0; ind < space_dim - 1; ++ind)
    {
      normal_vec = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(ind);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          loc_lambda[i][1 + ind] += normal_vec[dim] * glob_lambda[i][dim];
          loc_lambda[i][1 + ind + space_dim] += normal_vec[dim] * glob_lambda[i][space_dim + dim];
        }
    }

    return loc_lambda;
  }

  template <class hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  inline SmallMatOutT& edge_dof_to_node_dof(const SmallMatInT& loc_lambda,
                                            SmallMatOutT& glob_lambda,
                                            hyEdgeT& hyper_edge) const
  {
    Point<space_dim, lSol_float_t> normal_vec =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        glob_lambda[i][dim] += normal_vec[dim] * loc_lambda[i][0];
        glob_lambda[i][dim + space_dim] += normal_vec[dim] * loc_lambda[i][space_dim];
      }

    for (unsigned int ind = 0; ind < space_dim - 1; ++ind)
    {
      normal_vec = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(ind);

      for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
        for (unsigned int dim = 0; dim < space_dim; ++dim)
        {
          glob_lambda[i][dim] += normal_vec[dim] * loc_lambda[i][1 + ind];
          glob_lambda[i][dim + space_dim] += normal_vec[dim] * loc_lambda[i][1 + ind + space_dim];
        }
    }

    return glob_lambda;
  }

  /*!***********************************************************************************************
   * \brief   Assemble and factorize the displacement Schur systems of the four components.
   *
   * Builds the shared per-edge integral blocks M (diagonal), G = int (grad phi) phi, B = boundary
   * normal flux and F = boundary mass; eliminates each triplet's sigma and wdot blocks (both pure
   * mass) explicitly to obtain the displacement Schur operator
   *   S(C_sig, C_u) = theta*tau*F + theta*C_sig*(B-G) M^{-1} G + (C_u/(theta*dt^2)) M.
   * Components A,B are single triplets (nw_ x nw_); C,D couple a force group (C_sig=C_sig_n) with a
   * rotation group via the cross product into an n2w_ x n2w_ system. M_inv, G and (B-G) are cached
   * for the reduction / back-substitution in solve_local_problem.
   ************************************************************************************************/
  template <typename hyEdgeT>
  inline void assemble_schur(hyEdgeT& hyper_edge) const
  {
    auto& data = hyper_edge.data;
    constexpr unsigned int n = n_shape_fct_;
    data.extra = get_extra_coeffs(hyper_edge);  // C_n, C_m, C_u, C_r

    std::array<lSol_float_t, nw_> M_diag;
    std::array<lSol_float_t, nw_ * nw_> G, Bm, F;  // row-major; Bm = B - G
    M_diag.fill(0.);
    G.fill(0.);
    Bm.fill(0.);
    F.fill(0.);
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j)
      {
        const lSol_float_t vol =
          integrator::template integrate_vol_phiphi(i, j, hyper_edge.geometry);
        const auto grad =
          integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dimT, lSol_float_t>,
                                                         decltype(hyEdgeT::geometry)>(
            i, j, hyper_edge.geometry);
        lSol_float_t face_integral = 0., normal_integral = 0.;
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        {
          const lSol_float_t h =
            integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
              i, j, face, hyper_edge.geometry);
          face_integral += h;
          normal_integral += h * hyper_edge.geometry.local_normal(face).operator[](0);
        }
        if (i == j)
          M_diag[i] = vol;  // M is diagonal: orthogonal Legendre basis on affine geometry
        G[i * n + j] = grad[0];
        Bm[i * n + j] = normal_integral - grad[0];
        F[i * n + j] = face_integral;
      }
    for (unsigned int i = 0; i < n; ++i)
      data.M_inv[i] = 1. / M_diag[i];
    data.Gmat = G;
    data.BmG = Bm;

    // K = (B-G) M^{-1} G  (row-major), the dense part of the Schur operator.
    std::array<lSol_float_t, nw_ * nw_> K;
    for (unsigned int r = 0; r < n; ++r)
      for (unsigned int c = 0; c < n; ++c)
      {
        lSol_float_t acc = 0.;
        for (unsigned int k = 0; k < n; ++k)
          acc += Bm[r * n + k] * data.M_inv[k] * G[k * n + c];
        K[r * n + c] = acc;
      }

    // Entry (r,c) of the single-triplet Schur operator S(C_sig, C_u) in stage normalization:
    // eliminating q and z from assemble_loc_matrix<true> leaves S = tau*F + C_sig*K + sigma^2*C_u*M.
    const lSol_float_t sig = sigma();
    auto S_entry = [&](lSol_float_t C_sig, lSol_float_t C_u, unsigned int r, unsigned int c) {
      lSol_float_t v = tau_ * F[r * n + c] + C_sig * K[r * n + c];
      if (r == c)
        v += sig * sig * C_u * M_diag[r];
      return v;
    };

    // Single component: factorize S directly (column-major for LAPACK).
    auto fill_single = [&](std::array<lSol_float_t, nw_ * nw_>& lu, unsigned int dim) {
      const lSol_float_t C_sig = data.extra[sig_blk(dim)];
      const lSol_float_t C_u = data.extra[w_blk(dim)];
      for (unsigned int c = 0; c < nw_; ++c)
        for (unsigned int r = 0; r < nw_; ++r)
          lu[c * nw_ + r] = S_entry(C_sig, C_u, r, c);
    };
    fill_single(data.lu_A, compA_dim);
    Wrapper::lapack_factorize(nw_, data.lu_A.data(), data.ipiv_A.data());
    fill_single(data.lu_B, compB_dim);
    Wrapper::lapack_factorize(nw_, data.lu_B.data(), data.ipiv_B.data());

    // Coupled component (force group dim_n, rotation group dim_r), cross sign s as it appears in
    // the full matrix (sigma_n row: +s*M*w_r; w_r row: +s*M*sigma_n up to sign). Eliminating
    // sigma_n = C_sn*(M^{-1}b + M^{-1}G w_n - s*w_r) puts -s on both off-diagonal Schur blocks:
    //   [ S_n            -s*C_sig_n*(B-G) ]
    //   [ -s*C_sig_n*G   S_r + C_sig_n*M  ]
    auto fill_coupled = [&](std::array<lSol_float_t, n2w_ * n2w_>& lu, lSol_float_t s,
                            unsigned int dim_n, unsigned int dim_r) {
      const lSol_float_t Csn = data.extra[sig_blk(dim_n)], Cun = data.extra[w_blk(dim_n)];
      const lSol_float_t Csr = data.extra[sig_blk(dim_r)], Cur = data.extra[w_blk(dim_r)];
      for (unsigned int c = 0; c < n2w_; ++c)
        for (unsigned int r = 0; r < n2w_; ++r)
        {
          lSol_float_t v;
          if (r < nw_ && c < nw_)
            v = S_entry(Csn, Cun, r, c);
          else if (r < nw_ && c >= nw_)
            v = -s * Csn * Bm[r * n + (c - nw_)];
          else if (r >= nw_ && c < nw_)
            v = -s * Csn * G[(r - nw_) * n + c];
          else
          {
            v = S_entry(Csr, Cur, r - nw_, c - nw_);
            if (r == c)
              v += Csn * M_diag[r - nw_];
          }
          lu[c * n2w_ + r] = v;
        }
    };
    fill_coupled(data.lu_C, +1., compC_dim_n, compC_dim_r);
    Wrapper::lapack_factorize(n2w_, data.lu_C.data(), data.ipiv_C.data());
    fill_coupled(data.lu_D, -1., compD_dim_n, compD_dim_r);
    Wrapper::lapack_factorize(n2w_, data.lu_D.data(), data.ipiv_D.data());
  }

  /*!***********************************************************************************************
   * \brief   Reduce triplet \c dim to its displacement right-hand side; stash mass-scaled rhs.
   *
   * rhs_w = b_w - C_sig*(B-G) M^{-1} b_sig - sigma*b_wdot (stage normalization). Outputs
   * M^{-1} b_sig and M^{-1} b_wdot for the back-substitution.
   ************************************************************************************************/
  template <typename DataT>
  inline std::array<lSol_float_t, nw_> triplet_reduce(
    const DataT& data, unsigned int dim, const SmallVec<n_loc_dofs_, lSol_float_t>& rhs,
    std::array<lSol_float_t, nw_>& m_inv_bsig, std::array<lSol_float_t, nw_>& m_inv_bwdot) const
  {
    constexpr unsigned int n = n_shape_fct_;
    const lSol_float_t C_sig = data.extra[sig_blk(dim)];
    std::array<lSol_float_t, nw_> rhs_w, bw;
    for (unsigned int k = 0; k < n; ++k)
    {
      m_inv_bsig[k] = data.M_inv[k] * rhs[sig_blk(dim) * n + k];
      m_inv_bwdot[k] = data.M_inv[k] * rhs[wdot_blk(dim) * n + k];
      bw[k] = rhs[w_blk(dim) * n + k];
    }
    for (unsigned int r = 0; r < n; ++r)
    {
      lSol_float_t bmg_x = 0.;
      for (unsigned int c = 0; c < n; ++c)
        bmg_x += data.BmG[r * n + c] * m_inv_bsig[c];
      rhs_w[r] = bw[r] - C_sig * bmg_x - sigma() * rhs[wdot_blk(dim) * n + r];
    }
    return rhs_w;
  }

  /*!***********************************************************************************************
   * \brief   Recover sigma, w, wdot of triplet \c dim from the solved displacement w and scatter.
   *
   * sigma_f = C_sig (M^{-1} b_sig + M^{-1} G w - s*w_partner),  wdot = M^{-1} b_wdot +
   * sigma*C_u*w (stage normalization; C_u = 0 on massless welds gives wdot from b_wdot alone).
   * \c s_cross is 0 for single components; for a coupled force group it is the cross sign of the
   * full matrix (+s*M*w_partner on the sigma_f row) and \c w_partner the rotation group's
   * displacement.
   ************************************************************************************************/
  template <typename DataT>
  inline void triplet_backsub(const DataT& data, unsigned int dim, lSol_float_t s_cross,
                              const std::array<lSol_float_t, nw_>& w,
                              const std::array<lSol_float_t, nw_>& w_partner,
                              const std::array<lSol_float_t, nw_>& m_inv_bsig,
                              const std::array<lSol_float_t, nw_>& m_inv_bwdot,
                              SmallVec<n_loc_dofs_, lSol_float_t>& result) const
  {
    constexpr unsigned int n = n_shape_fct_;
    const lSol_float_t C_sig = data.extra[sig_blk(dim)];
    const lSol_float_t C_u = data.extra[w_blk(dim)];
    for (unsigned int k = 0; k < n; ++k)
    {
      lSol_float_t Gw = 0.;
      for (unsigned int c = 0; c < n; ++c)
        Gw += data.Gmat[k * n + c] * w[c];
      result[sig_blk(dim) * n + k] =
        C_sig * (m_inv_bsig[k] + data.M_inv[k] * Gw - s_cross * w_partner[k]);
      result[w_blk(dim) * n + k] = w[k];
      result[wdot_blk(dim) * n + k] = m_inv_bwdot[k] + sigma() * C_u * w[k];
    }
  }

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> solve_local_problem(const SmallMatT& lambda_values,
                                                                 const unsigned int solution_type,
                                                                 hyEdgeT& hyper_edge,
                                                                 const lSol_float_t time,
                                                                 const unsigned int stage = 0) const
  {
    try
    {
      auto& data = hyper_edge.data;

      SmallVec<n_loc_dofs_, lSol_float_t> rhs;
      if (solution_type == 0)
        rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge);
      else if (solution_type == 1)
        rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge) +
              assemble_rhs_stage(hyper_edge, time, stage);
      else
        hy_assert(0 == 1, "This has not been implemented!");

      if (full_lu_)
      {
        if (!data.full_lu_factorized)
        {
          data.full_lu = assemble_loc_matrix<true>(hyper_edge, sigma());
          Wrapper::lapack_factorize(n_loc_dofs_, data.full_lu.data().data(),
                                    data.full_ipiv.data());
          data.full_lu_factorized = true;
        }
        Wrapper::lapack_solve_factored(n_loc_dofs_, 1, data.full_lu.data().data(),
                                       data.full_ipiv.data(), rhs.data().data());
        return rhs;
      }

      if (!data.loc_mat_factorized)
      {
        assemble_schur(hyper_edge);
        data.loc_mat_factorized = true;
      }

      constexpr unsigned int n = n_shape_fct_;
      SmallVec<n_loc_dofs_, lSol_float_t> result;
      const std::array<lSol_float_t, nw_> dummy{};  // unused w_partner for single components

      // Single components A, B.
      auto solve_single = [&](std::array<lSol_float_t, nw_ * nw_>& lu, std::array<int, nw_>& ipiv,
                              unsigned int dim) {
        std::array<lSol_float_t, nw_> mis, miw;
        auto rw = triplet_reduce(data, dim, rhs, mis, miw);
        Wrapper::lapack_solve_factored(nw_, 1, lu.data(), ipiv.data(), rw.data());
        triplet_backsub(data, dim, 0., rw, dummy, mis, miw, result);
      };
      solve_single(data.lu_A, data.ipiv_A, compA_dim);
      solve_single(data.lu_B, data.ipiv_B, compB_dim);

      // Coupled components C (s=+1), D (s=-1): force group dim_n, rotation group dim_r.
      auto solve_coupled = [&](std::array<lSol_float_t, n2w_ * n2w_>& lu,
                               std::array<int, n2w_>& ipiv, lSol_float_t s, unsigned int dim_n,
                               unsigned int dim_r) {
        std::array<lSol_float_t, nw_> mis_n, miw_n, mis_r, miw_r;
        auto rw_n = triplet_reduce(data, dim_n, rhs, mis_n, miw_n);
        auto rw_r = triplet_reduce(data, dim_r, rhs, mis_r, miw_r);
        const lSol_float_t Csn = data.extra[sig_blk(dim_n)];
        std::array<lSol_float_t, n2w_> rw;
        for (unsigned int k = 0; k < n; ++k)
        {
          rw[k] = rw_n[k];
          rw[nw_ + k] = rw_r[k] + s * Csn * rhs[sig_blk(dim_n) * n + k];
        }
        Wrapper::lapack_solve_factored(n2w_, 1, lu.data(), ipiv.data(), rw.data());
        std::array<lSol_float_t, nw_> w_n, w_r;
        for (unsigned int k = 0; k < n; ++k)
        {
          w_n[k] = rw[k];
          w_r[k] = rw[nw_ + k];
        }
        triplet_backsub(data, dim_n, s, w_n, w_r, mis_n, miw_n, result);  // force group: +s*w_r
        triplet_backsub(data, dim_r, 0., w_r, dummy, mis_r, miw_r, result);  // rotation group
      };
      solve_coupled(data.lu_C, data.ipiv_C, +1., compC_dim_n, compC_dim_r);
      solve_coupled(data.lu_D, data.ipiv_D, -1., compD_dim_n, compD_dim_r);

      return result;
    }
    catch (Wrapper::LAPACKexception& exc)
    {
      std::cout << hyper_edge.geometry.area() << std::endl;
      hy_assert(0 == 1, exc.what() << std::endl
                                   << "This can happen if quadrature is too inaccurate!");
      throw exc;
    }
  }

  template <typename hyEdgeT>
  inline SmallMat<2 * hyEdge_dimT, 2 * space_dim * n_shape_bdr_, lSol_float_t>
  extract_fluxes_from_coeffs(const SmallVec<n_loc_dofs_, lSol_float_t>& coeffs,
                             hyEdgeT& hyper_edge) const
  {
    SmallMat<2 * hyEdge_dimT, 2 * space_dim * n_shape_bdr_, lSol_float_t> bdr_values;
    lSol_float_t integral;

    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int j = 0; j < n_shape_bdr_; ++j)
        for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        {
          integral = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
            i, j, face, hyper_edge.geometry);
          for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
            bdr_values(face, j + dim) +=
              hyper_edge.geometry.local_normal(face).operator[](0) * integral *
                coeffs[dim * n_shape_fct_ + i] +
              tau_ * coeffs[(2 * space_dim + dim) * n_shape_fct_ + i] * integral;
        }

    return bdr_values;
  }

  /*!***********************************************************************************************
   * \brief   Zero every dynamic-Dirichlet trace dof of \c lambda in place (leaves bit-6 faces).
   *
   * A face flagged static-only (bit 6) keeps its trace; on every other face the components whose
   * Dirichlet bit is set are constrained, so their trace contribution is removed. Shared by the
   * two flux entry points below (and mirrors the masking done in the error evaluation).
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatT>
  void zero_dirichlet_trace(SmallMatT& lambda, hyEdgeT& hyper_edge) const
  {
    for (unsigned int node = 0; node < 2 * hyEdge_dimT; ++node)
      for (unsigned int dof = 0; dof < 2 * space_dim; ++dof)
        if (is_dirichlet(hyper_edge, node, dof))
          lambda[node][dof] = 0.;
  }

  /*!***********************************************************************************************
   * \brief   Shared core of trace_to_flux / residual_flux: apply the condensed local operator.
   *
   * \c solution_type selects the local-solve right-hand side (see solve_local_problem): 0 is the
   * homogeneous operator action A*lambda (used to assemble the time-constant system matrix); 1 is
   * the residual, i.e. operator action plus the old-time-step / body-load / Dirichlet data. The
   * Dirichlet trace dofs are removed on input and output so constrained faces stay pinned.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& apply_local_flux(const SmallMatInT& lambda_values_in,
                                 SmallMatOutT& lambda_values_out,
                                 const unsigned int solution_type,
                                 hyEdgeT& hyper_edge,
                                 const lSol_float_t time = 0.) const
  {
    hy_assert(lambda_values_in.size() == lambda_values_out.size() &&
                lambda_values_in.size() == 2 * hyEdge_dimT,
              "Both matrices must be of same size which corresponds to the number of faces!");
    for (unsigned int i = 0; i < lambda_values_in.size(); ++i)
      hy_assert(
        lambda_values_in[i].size() == lambda_values_out[i].size() &&
          lambda_values_in[i].size() == n_glob_dofs_per_node(),
        "Both matrices must be of same size which corresponds to the number of dofs per face!");

    SmallMatInT lambda_in = lambda_values_in;
    zero_dirichlet_trace(lambda_in, hyper_edge);

    SmallMatInT lambda_values_loc = node_dof_to_edge_dof(lambda_in, hyper_edge);

    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values_loc, solution_type, hyper_edge, time);

    auto result = extract_fluxes_from_coeffs(coeffs, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < 2 * space_dim; ++j)
        lambda_values_loc[i][j] = tau_ * lambda_values_loc[i][j] - result(i, j);

    lambda_values_out = edge_dof_to_node_dof(lambda_values_loc, lambda_values_out, hyper_edge);

    zero_dirichlet_trace(lambda_values_out, hyper_edge);

    return lambda_values_out;
  }

  /*!***********************************************************************************************
   * \brief   Homogeneous condensed operator action (assembles the time-constant system matrix).
   *
   * A distinct entry point from residual_flux (kept separate for historical reasons); both
   * delegate to apply_local_flux, differing only in the local-solve rhs (solution_type 0 vs 1).
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& trace_to_flux(const SmallMatInT& lambda_values_in,
                              SmallMatOutT& lambda_values_out,
                              hyEdgeT& hyper_edge,
                              const lSol_float_t time = 0.) const
  {
    return apply_local_flux(lambda_values_in, lambda_values_out, 0U, hyper_edge, time);
  }

  template <typename hyEdgeT>
  bool is_dirichlet(hyEdgeT& hyper_edge, unsigned int node, unsigned int dof) const
  {
    if (hyper_edge.node_descriptor[node] & (1ul << 6)) return false;  // static-only flag
    if (dof >= 2 * space_dim) return false;
    return hyper_edge.node_descriptor[node] & (1ul << dof);
  }

  /*!***********************************************************************************************
   * \brief   Residual of the condensed operator (operator action plus old-step / load data).
   *
   * A distinct entry point from trace_to_flux (kept separate for historical reasons); see
   * apply_local_flux for the shared body.
   ************************************************************************************************/
  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& residual_flux(const SmallMatInT& lambda_values_in,
                              SmallMatOutT& lambda_values_out,
                              hyEdgeT& hyper_edge,
                              const lSol_float_t time = 0.) const
  {
    return apply_local_flux(lambda_values_in, lambda_values_out, 1U, hyper_edge, time);
  }

  /*!***********************************************************************************************
   * \brief   Local squared contribution to the L2 error.
   *
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which analytic functions are evaluated.
   * \retval  vec_b             Local part of vector b.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<lSol_float_t, 2U> errors(
    const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>&
      lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    std::array<lSol_float_t,3> comps = {1,-1,-2};
    std::array<lSol_float_t, n_shape_fct_> coeffs;
    std::array<lSol_float_t, n_shape_bdr_> bcoeffs;
    lSol_float_t error = 0, trace = 0;
    // Weight each endpoint contribution by this edge's length: summed over all edges this
    // accumulates w_v = sum_{e per v} l_e at every node, i.e. the HDG skeleton norm
    // sum_e l_e ||.||^2_{de} -- a function-space (lumped L2) norm instead of a bare nodal
    // l2 sum whose value scales with the node count.
    const lSol_float_t len = hyper_edge.geometry.area();
    SmallVec<space_dim*n_shape_fct_, lSol_float_t> u_old = hyper_edge.data.u_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t> r_old = hyper_edge.data.r_old;

    // Input lambdas are in node-frame (global). Convert to edge-frame so we can compare against
    // analytic_result_u/phi which is evaluated against edge-local normals (comps = {1,-1,-2}).
    auto lambda_loc = node_dof_to_edge_dof(lambda_values, hyper_edge);

    for (unsigned int dim = 0; dim < space_dim; dim++) {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = u_old[i + dim * n_shape_fct_];
      error += integrator::template integrate_vol_diffsquare_discanacomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, comps[dim],
                                                                         hyper_edge.geometry, time);
    }

    for (unsigned int dim = 0; dim < space_dim; dim++) {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = r_old[i + dim * n_shape_fct_];
      error += integrator::template integrate_vol_diffsquare_discanacomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::analytic_result_phi, Point<hyEdge_dimT, lSol_float_t>>(coeffs, comps[dim],
                                                                         hyper_edge.geometry, time);
    }

    // Lambda is zeroed at dynamic-Dirichlet DOFs (see residual_flux), so for those DOFs we
    // substitute the projected Dirichlet value here -- otherwise the trace error picks up the
    // full analytic value at the boundary. Bit-6 faces (static-only dirichlet) keep their
    // lambda since make_initial_from_static fills it with the projected bulk solution.
    // For n_shape_bdr_ == 1 the projection coefficient equals dirichlet_value at the face point.
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dimT; ++bdr) {
      const bool dyn_dirichlet = hyper_edge.node_descriptor[bdr]
                                 && !(hyper_edge.node_descriptor[bdr] & (1u << 6));
      for (unsigned int dim = 0; dim < space_dim; ++dim) {
        for (unsigned int i = 0; i < bcoeffs.size(); ++i) {
          const unsigned int dof_j = dim * n_shape_bdr_ + i;
          if (dyn_dirichlet && (hyper_edge.node_descriptor[bdr] & (1u << dof_j)))
            bcoeffs[i] = integrator::template integrate_bdr_phivecfunccomp<
                Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
                decltype(hyEdgeT::geometry), parameters::dirichlet_value_u,
                Point<hyEdge_dimT, lSol_float_t>>(i, bdr, comps[dim], hyper_edge.geometry, time);
          else
            bcoeffs[i] = lambda_loc[bdr][i + dim * n_shape_bdr_];
        }
        lSol_float_t contrib = integrator::template integrate_bdr_diffsquare_discanacomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
          parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(
            bcoeffs, bdr, comps[dim], hyper_edge.geometry, time);
        trace += len * contrib;
      }
    }

    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dimT; ++bdr) {
      const bool dyn_dirichlet = hyper_edge.node_descriptor[bdr]
                                 && !(hyper_edge.node_descriptor[bdr] & (1u << 6));
      for (unsigned int dim = 0; dim < space_dim; ++dim) {
        for (unsigned int i = 0; i < bcoeffs.size(); ++i) {
          const unsigned int dof_j = (3 + dim) * n_shape_bdr_ + i;
          if (dyn_dirichlet && (hyper_edge.node_descriptor[bdr] & (1u << dof_j)))
            bcoeffs[i] = integrator::template integrate_bdr_phivecfunccomp<
                Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
                decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi,
                Point<hyEdge_dimT, lSol_float_t>>(i, bdr, comps[dim], hyper_edge.geometry, time);
          else
            bcoeffs[i] = lambda_loc[bdr][i + (3+dim) * n_shape_bdr_];
        }
        lSol_float_t contrib = integrator::template integrate_bdr_diffsquare_discanacomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
          parameters::analytic_result_phi, Point<hyEdge_dimT, lSol_float_t>>(
            bcoeffs, bdr, comps[dim], hyper_edge.geometry, time);
        trace += len * contrib;
      }
    }

    return std::array<lSol_float_t, 2U>({error, trace});
  }

  /*!***********************************************************************************************
   * \brief   Local squared contribution to the L2 norm of the analytic solution.
   *
   * Evaluates the same integrals as \c errors, but with a vanishing discrete solution, so that
   * the result is the squared L2 norm of the analytic solution. This allows to compute relative
   * errors.
   *
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time at which analytic functions are evaluated.
   * \retval  norm              Local squared L2 norm of the analytic solution.
   ************************************************************************************************/
  template <class hyEdgeT>
  std::array<lSol_float_t, 2U> norms(
    const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>&
      lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    (void)lambda_values;

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    std::array<lSol_float_t,3> comps = {1,-1,-2};
    std::array<lSol_float_t, n_shape_fct_> coeffs;
    coeffs.fill(0.);
    lSol_float_t norm = 0;

    for (unsigned int dim = 0; dim < 3; dim++) {
      norm += integrator::template integrate_vol_diffsquare_discanacomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, comps[dim],
                                                                         hyper_edge.geometry, time);
    }

    for (unsigned int dim = 0; dim < 3; dim++) {
      norm += integrator::template integrate_vol_diffsquare_discanacomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::analytic_result_phi, Point<hyEdge_dimT, lSol_float_t>>(coeffs, comps[dim],
                                                                         hyper_edge.geometry, time);
    }

    // trace norm of the analytic solution in the same edge-length-weighted skeleton norm as
    // the trace error in errors(): zero boundary coefficients turn the diff-square integral
    // into the square of the analytic value at the endpoints
    std::array<lSol_float_t, n_shape_bdr_> bcoeffs;
    bcoeffs.fill(0.);
    const lSol_float_t len = hyper_edge.geometry.area();
    lSol_float_t trace = 0;
    for (unsigned int bdr = 0; bdr < 2 * hyEdge_dimT; ++bdr)
      for (unsigned int dim = 0; dim < space_dim; ++dim) {
        trace += len * integrator::template integrate_bdr_diffsquare_discanacomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
          parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(
            bcoeffs, bdr, comps[dim], hyper_edge.geometry, time);
        trace += len * integrator::template integrate_bdr_diffsquare_discanacomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
          parameters::analytic_result_phi, Point<hyEdge_dimT, lSol_float_t>>(
            bcoeffs, bdr, comps[dim], hyper_edge.geometry, time);
      }

    return std::array<lSol_float_t, 2U>({norm, trace});
  }

  static constexpr unsigned int n_energy_components() { return 6 * space_dim; }

  // Per-edge energy split into 6*space_dim components, ordered (block of size space_dim each):
  //   0: ½ ∫ n²/C_n   1: ½ ∫ m²/C_m   2: ½ ∫ v²/C_u   3: ½ ∫ s²/C_r
  //   4: ½ τ Σ_bdr ∫ (u-λ_u)²            5: ½ τ Σ_bdr ∫ (r-λ_r)²
  template <class hyEdgeT>
  std::array<lSol_float_t, n_energy_components()> energy(
    const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>&
      lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t /*time*/ = 0.) const
  {
    std::array<lSol_float_t, n_energy_components()> result;
    result.fill(0.);

    auto& n_old = hyper_edge.data.n_old;
    auto& m_old = hyper_edge.data.m_old;
    auto& u_old = hyper_edge.data.u_old;
    auto& r_old = hyper_edge.data.r_old;
    auto& v_old = hyper_edge.data.v_old;
    auto& s_old = hyper_edge.data.s_old;

    auto extra = get_extra_coeffs(hyper_edge);
    auto lambda_loc = node_dof_to_edge_dof(lambda_values, hyper_edge);

    for (unsigned int d = 0; d < space_dim; ++d) {
      const lSol_float_t Cn = extra[0 * space_dim + d];
      const lSol_float_t Cm = extra[1 * space_dim + d];
      const lSol_float_t Cu = extra[2 * space_dim + d];
      const lSol_float_t Cr = extra[3 * space_dim + d];

      lSol_float_t strain_n = 0, strain_m = 0, kin_v = 0, kin_s = 0;
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        for (unsigned int j = 0; j < n_shape_fct_; ++j) {
          const lSol_float_t mij =
            integrator::template integrate_vol_phiphi<decltype(hyEdgeT::geometry)>(
              i, j, hyper_edge.geometry);
          strain_n += n_old[d * n_shape_fct_ + i] * n_old[d * n_shape_fct_ + j] * mij;
          strain_m += m_old[d * n_shape_fct_ + i] * m_old[d * n_shape_fct_ + j] * mij;
          kin_v    += v_old[d * n_shape_fct_ + i] * v_old[d * n_shape_fct_ + j] * mij;
          kin_s    += s_old[d * n_shape_fct_ + i] * s_old[d * n_shape_fct_ + j] * mij;
        }

      // massless welds have C == 0 -> 0/0 = NaN; zero mass carries zero energy
      result[0 * space_dim + d] = Cn > 0. ? 0.5 * strain_n / Cn : 0.;
      result[1 * space_dim + d] = Cm > 0. ? 0.5 * strain_m / Cm : 0.;
      result[2 * space_dim + d] = Cu > 0. ? 0.5 * kin_v / Cu : 0.;
      result[3 * space_dim + d] = Cr > 0. ? 0.5 * kin_s / Cr : 0.;

      lSol_float_t hyb_u = 0, hyb_r = 0;
      for (unsigned int bdr = 0; bdr < 2 * hyEdge_dimT; ++bdr) {
        // (y - λ)² = y² - 2 y λ + λ² on face bdr, with y = u or r.
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
          for (unsigned int j = 0; j < n_shape_fct_; ++j) {
            const lSol_float_t mij =
              integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            hyb_u += u_old[d * n_shape_fct_ + i] * u_old[d * n_shape_fct_ + j] * mij;
            hyb_r += r_old[d * n_shape_fct_ + i] * r_old[d * n_shape_fct_ + j] * mij;
          }
        for (unsigned int i = 0; i < n_shape_fct_; ++i)
          for (unsigned int j = 0; j < n_shape_bdr_; ++j) {
            const lSol_float_t mij =
              integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                i, j, bdr, hyper_edge.geometry);
            hyb_u -= 2 * u_old[d * n_shape_fct_ + i]
                       * lambda_loc[bdr][j + d * n_shape_bdr_] * mij;
            hyb_r -= 2 * r_old[d * n_shape_fct_ + i]
                       * lambda_loc[bdr][j + (space_dim + d) * n_shape_bdr_] * mij;
          }
        // hyEdge_dimT==1 ⇒ trace is 0-dimensional, ψ≡1 ⇒ ∫_∂e λ² = λ² directly.
        static_assert(hyEdge_dimT == 1, "trace-square shortcut only valid for hyEdge_dim==1");
        const lSol_float_t lam_u = lambda_loc[bdr][d * n_shape_bdr_];
        const lSol_float_t lam_r = lambda_loc[bdr][(space_dim + d) * n_shape_bdr_];
        hyb_u += lam_u * lam_u;
        hyb_r += lam_r * lam_r;
      }
      result[4 * space_dim + d] = 0.5 * tau_ * hyb_u;
      result[5 * space_dim + d] = 0.5 * tau_ * hyb_r;
    }

    return result;
  }

  // Edge-local frame vector: idx == 0 → inner_normal(0),
  // idx == -k (k>=1) → outer_normal(k-1).
  // Convention: nonneg → inner, neg → outer.
  template <class hyEdgeT, typename float_t>
  Point<space_dim, float_t>
  edge_frame_vector(hyEdgeT& hyper_edge, int idx) const
  {
    if (idx >= 0)
      return (Point<space_dim, float_t>)hyper_edge.geometry.inner_normal(idx);
    else
      return (Point<space_dim, float_t>)hyper_edge.geometry.outer_normal(-idx - 1);
  }

  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dimension()>
  bulk_values(const std::array<abscissa_float_t, sizeT>& abscissas,
              const input_array_t& lambda_values,
              hyEdgeT& hyper_edge,
              const lSol_float_t time = 0.) const
  {
    constexpr unsigned int n_pts    = Hypercube<hyEdge_dimT>::pow(sizeT);
    constexpr unsigned int n_fields = 6;  // n, m, u, r, v, s

   // Read stored fields directly (edge-local frame).
    // Order matches field index c: 0=n, 1=m, 2=u, 3=r, 4=v, 5=s.
    const SmallVec<space_dim*n_shape_fct_, lSol_float_t>* old_fields[n_fields] = {
      &hyper_edge.data.n_old,
      &hyper_edge.data.m_old,
      &hyper_edge.data.u_old,
      &hyper_edge.data.r_old,
      &hyper_edge.data.v_old,
      &hyper_edge.data.s_old,
    };

    SmallVec<n_shape_fct_, lSol_float_t> coeffs;
    SmallVec<static_cast<unsigned int>(sizeT), abscissa_float_t> helper(abscissas);

    // Edge-local-frame component values at each abscissa, per field.
    // shape: [field][local_component][pt]
    std::array<std::array<std::array<lSol_float_t, n_pts>, space_dim>, n_fields> point_vals{};

    for (unsigned int c = 0; c < n_fields; ++c)
      for (unsigned int dim = 0; dim < space_dim; ++dim) {
        for (unsigned int i = 0; i < coeffs.size(); ++i)
          coeffs[i] = (*old_fields[c])[dim * n_shape_fct_ + i];
        for (unsigned int pt = 0; pt < n_pts; ++pt)
          point_vals[c][dim][pt] = integrator::shape_fun_t::template lin_comb_fct_val<float>(
            coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT>>(pt, helper)
          );
      }

    std::array<std::array<lSol_float_t, n_pts>, system_dimension()> result{};

    static_assert(2 <= n_fields);

    // fields n m should be in local frame -> stay
    for (unsigned int c = 0; c < 2; ++c)
      for (unsigned int local = 0; local < space_dim; ++local)
        for (unsigned int q = 0; q < n_pts; ++q)
          result[c * space_dim + local][q] = point_vals[c][local][q];

    // fields u r should be in global frame -> need to transform
    for (unsigned int c = 2; c < n_fields; ++c)
      for (unsigned int local = 0; local < space_dim; ++local) {
        const int idx = -static_cast<int>(local);  // 0, -1, -2
        Point<space_dim, lSol_float_t> nv =
          edge_frame_vector<hyEdgeT, lSol_float_t>(hyper_edge, idx);
        for (unsigned int dim = 0; dim < space_dim; ++dim)
          for (unsigned int q = 0; q < n_pts; ++q)
            result[c * space_dim + dim][q] += point_vals[c][local][q] * nv[dim];
      }

    return result;
  }

  template <class hyEdgeT>
  SmallVec<4*space_dim, lSol_float_t>
  get_extra_coeffs(hyEdgeT& hyper_edge) const
  {
    SmallVec<4 * space_dim, lSol_float_t> extra_coeffs(1.); // C_n, C_m, C_u, C_r

    if (!hyper_edge.geometry.has_extra_data()) return extra_coeffs;

    // 17 quantities associated with each fiber
    // column      header                  description
    //  0          mass,                   mass
    //  1, 2, 3    EA,kG_1A,kG_2A,         displacement stiffness
    //  4, 5, 6    G_xI_x,E_1I_1,E_2I_2,   rotation stiffness
    //  7, 8, 9    n_11,n_12,n_13,         normal 1
    // 10,11,12    n_21,n_22,n_23,         normal 2
    // 13,14       width1,width2,          widths in direction of normals
    // 15          fiber_id                The fiber the edge is part of (-1 if a connection)
    // 16          fiber_edge_id           Where along the fiber the edge is (-1 if a connection)
    // material coefficients are given in the tangent,normal1,normal2 basis
    // so must be transformed into local basis chosen by HyperHDG, tangent coincides upto sign

    auto extra_data = hyper_edge.geometry.extra_data();
    hy_check(extra_data.size() == 17,
      "timowave expected 17 material coefficients, found " << extra_data.size());

    lSol_float_t mass = extra_data[0];
    SmallVec<space_dim, lSol_float_t> normal1 =
      std::array<lSol_float_t, space_dim>{{extra_data[7], extra_data[8], extra_data[9]}};
    SmallVec<space_dim, lSol_float_t> normal2 =
      std::array<lSol_float_t, space_dim>{{extra_data[10], extra_data[11], extra_data[12]}};
    SmallVec<space_dim, lSol_float_t> outer1 = hyper_edge.geometry.outer_normal(0);
    SmallVec<space_dim, lSol_float_t> outer2 = hyper_edge.geometry.outer_normal(1);

    auto length = hyper_edge.geometry.area();
    auto density = mass / (extra_data[13] * extra_data[14] * length);

    auto w1 = extra_data[13];
    auto w2 = extra_data[14];

    // approximation to second moment in direction of normal1, normal2
    // valid for rectangular cross section
    auto moment1 = w1*w1*w1 * w2 / 12;
    auto moment2 = w1 * w2*w2*w2 / 12;

    // displacement stiffness C_n
    extra_coeffs[0] = extra_data[1];
    extra_coeffs[1] = extra_data[2] * scalar_product(outer1, normal1) +
                      extra_data[3] * scalar_product(outer1, normal2);
    extra_coeffs[2] = extra_data[2] * scalar_product(outer2, normal1) +
                      extra_data[3] * scalar_product(outer2, normal2);

    // rotation stiffness C_m
    extra_coeffs[3] = extra_data[4];
    extra_coeffs[4] = extra_data[5] * scalar_product(outer1, normal1) +
                      extra_data[6] * scalar_product(outer1, normal2);
    extra_coeffs[5] = extra_data[5] * scalar_product(outer2, normal1) +
                      extra_data[6] * scalar_product(outer2, normal2);

    // displacement inertia C_u
    extra_coeffs[6] = mass / length;
    extra_coeffs[7] = mass / length;
    extra_coeffs[8] = mass / length;

    // rotation inertia
    extra_coeffs[9]  = density * (moment1 + moment2);
    extra_coeffs[10] = density * (moment1 * scalar_product(outer1, normal1) +
                                  moment2 * scalar_product(outer1, normal2));
    extra_coeffs[11] = density * (moment1 * scalar_product(outer2, normal1) +
                                  moment2 * scalar_product(outer2, normal2));

    for (unsigned int i = 0; i < extra_coeffs.size(); i++) {
      extra_coeffs[i] = std::abs(extra_coeffs[i]);
      // hy_check(std::isfinite(extra_coeffs[i]), "get_extra finite coeffs " << extra_coeffs);
      // hy_check(extra_coeffs[i] > 0, "get_extra zero coeffs " << extra_coeffs << " mass " << mass);
    }

    // hy_check(false, "extra_coeffs " << extra_coeffs);

    return extra_coeffs;
  }

  /*!***********************************************************************************************
   * \brief   Solve the stage-local problem from the stage trace zeta and stash its solution.
   *
   * \c lambda_values_in is the stage trace zeta of the Gauss step protocol; the stage-local
   * solution (q,y,z) is stashed in slot \c stage, the state advances in finalize_step (idempotent:
   * re-calling overwrites the slot).
   ************************************************************************************************/
  template <class hyEdgeT>
  void set_data(
    const std::array<std::array<lSol_float_t, 2*n_shape_bdr_*space_dim>, 2 * hyEdge_dimT>& lambda_values_in,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.,
    const unsigned int stage = 0) const
  {
    auto lambda_values = node_dof_to_edge_dof(lambda_values_in, hyper_edge);
    hyper_edge.data.stage_coeffs[stage] =
      solve_local_problem(lambda_values, 1U, hyper_edge, time, stage);
  }

  /*!***********************************************************************************************
   * \brief   Gauss protocol: combine the stashed stage solutions into the new state.
   *
   * Endpoint update y+ = stage_affine_ * y^n + sum_l Re(w_l y_l); at s = 1 that is
   * y+ = 2 y_1 - y^n. q = (n, m) is algebraic and linear in (y, lambda), so the same extrapolation
   * is exact for it at s = 1 (multi-stage recovers q via recover_dual instead). Explicitly driver-
   * called once all stage set_data calls are in -- no per-edge completion tally.
   ************************************************************************************************/
  template <class hyEdgeT>
  void finalize_step(hyEdgeT& hyper_edge) const
  {
    auto& data = hyper_edge.data;
    constexpr unsigned int nf = space_dim * n_shape_fct_;
    SmallVec<nf, lSol_float_t>* fields[6] = {&data.n_old, &data.m_old, &data.u_old,
                                             &data.r_old, &data.v_old, &data.s_old};
    for (unsigned int f = 0; f < 6; ++f)
      for (unsigned int i = 0; i < nf; ++i)
      {
        lSol_float_t val = stage_affine_ * (*fields[f])[i];
        for (unsigned int l = 0; l < n_stages; ++l)
          val += stage_w_[l] * data.stage_coeffs[l][f * nf + i];
        (*fields[f])[i] = val;
      }
  }

  /*!***********************************************************************************************
   * \brief   Dirichlet boundary data of the static problem, added to a local right-hand side.
   *
   * Any face with a nonzero node descriptor is treated as Dirichlet (the initializers' static
   * problems consider every constrained face, including bit-6 static-only ones): q-rows receive
   * -<g, w_q nu>, y-rows tau<g, w_y>, for g = dirichlet_value_u / _phi. Shared by make_initial and
   * make_initial_from_static. HACK (inherited): assumes dirichlet_value_* returns zero in
   * non-constrained directions.
   ************************************************************************************************/
  template <class hyEdgeT>
  void add_dirichlet_rhs_static(SmallVec<n_loc_dofs_, lSol_float_t>& rhs,
                                hyEdgeT& hyper_edge,
                                const lSol_float_t time = 0.) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        if (!hyper_edge.node_descriptor[face])
          continue;
        auto integrals = integrate_bdr_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u>(i, face, {1, -1, -2},
                                                                      hyper_edge.geometry, time);
        for (unsigned int comp = 0; comp < 3; ++comp)
        {
          rhs[(0 * space_dim + comp) * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * integrals[comp];
          rhs[(2 * space_dim + comp) * n_shape_fct_ + i] += tau_ * integrals[comp];
        }
        integrals = integrate_bdr_phivecfunccomp_beam<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi>(i, face, {1, -1, -2},
                                                                        hyper_edge.geometry, time);
        for (unsigned int comp = 0; comp < 3; ++comp)
        {
          rhs[(1 * space_dim + comp) * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * integrals[comp];
          rhs[(3 * space_dim + comp) * n_shape_fct_ + i] += tau_ * integrals[comp];
        }
      }
  }

  /*!***********************************************************************************************
   * \brief   Recover the algebraic dual (n, m) from the stored (u, r) and static boundary data.
   *
   * The constitutive rows of the static system read (M/C) q - G w (+ cross terms) = b_q, so
   * q = C * M^{-1} * (b_q - (A_static w)_q). \c rhs must hold the q-row boundary data (lambda +
   * Dirichlet contributions); only its first 2*space_dim blocks are read. The dual is algebraic --
   * slaved to (y, lambda) at the same instant -- so the same recovery serves the initializers now
   * and the Gauss output-time diagnostics later. (This also fixes the formerly missing material
   * factor C: for unit coefficients, i.e. the wave4-type tests, nothing changes.)
   ************************************************************************************************/
  template <class hyEdgeT>
  void recover_dual(hyEdgeT& hyper_edge, const SmallVec<n_loc_dofs_, lSol_float_t>& rhs) const
  {
    auto& data = hyper_edge.data;
    SmallVec<4 * space_dim * n_shape_fct_, lSol_float_t> coeffs_w;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
      {
        coeffs_w[(2 * space_dim + dim) * n_shape_fct_ + i] = data.u_old[i + dim * n_shape_fct_];
        coeffs_w[(3 * space_dim + dim) * n_shape_fct_ + i] = data.r_old[i + dim * n_shape_fct_];
      }
    const auto res = assemble_loc_matrix<false>(hyper_edge) * coeffs_w;
    const auto extra = get_extra_coeffs(hyper_edge);  // C_n, C_m in components 0 .. 2*space_dim-1

    SmallSquareMat<n_shape_fct_> mmat;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
        mmat(i, j) = integrator::template integrate_vol_phiphi(i, j, hyper_edge.geometry);

    SmallVec<n_shape_fct_> temp;
    for (unsigned int d = 0; d < 2 * space_dim; ++d)
    {
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        temp[i] = rhs[d * n_shape_fct_ + i] - res[d * n_shape_fct_ + i];
      temp = temp / mmat;
      auto& target = d < space_dim ? data.n_old : data.m_old;
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        target[(d % space_dim) * n_shape_fct_ + i] = extra[d] * temp[i];
    }
  }

  template <class hyEdgeT, typename SmallMatT>
  SmallMatT& make_initial(SmallMatT& lambda_values,
                          hyEdgeT& hyper_edge,
                          const lSol_float_t time = 0.) const
  {
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& u_old = hyper_edge.data.u_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& r_old = hyper_edge.data.r_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& v_old = hyper_edge.data.v_old;
    SmallVec<space_dim*n_shape_fct_, lSol_float_t>& s_old = hyper_edge.data.s_old;

    // first u then r in skeletal variables
    SmallVec<space_dim, lSol_float_t> helper;

    static_assert(lambda_values[0].size() == 2*space_dim);

    // Set skeltal variable!
    for (unsigned int i = 0; i < lambda_values.size(); ++i)
    {
        for (unsigned int j = 0; j < n_shape_bdr_; ++j) {
          helper = integrator::template integrate_bdrUni_psivecfunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
            decltype(hyEdgeT::geometry), parametersT<space_dim, lSol_float_t>::initial_u,
            Point<hyEdge_dimT, lSol_float_t>>(j, i, hyper_edge.geometry, time);
          for (unsigned int dim = 0; dim < space_dim; dim++)
            lambda_values[i][j+dim] = helper[dim];
           helper = integrator::template integrate_bdrUni_psivecfunc<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
            decltype(hyEdgeT::geometry), parametersT<space_dim, lSol_float_t>::initial_r,
            Point<hyEdge_dimT, lSol_float_t>>(j, i, hyper_edge.geometry, time);
          for (unsigned int dim = 0; dim < space_dim; dim++)
            lambda_values[i][j+space_dim+dim] = helper[dim];
        }
    }

    for (unsigned int i = 0; i < lambda_values.size(); ++i)
      for (unsigned int j = 0; j < 2*space_dim; ++j)
        if (hyper_edge.node_descriptor[i] & (1<<j))
            lambda_values[i][j] = 0.;

    auto lambda_values_loc = node_dof_to_edge_dof(lambda_values, hyper_edge);

    // set u, r, n, m old
    // NOTE: computing in global dofs
    // Define primary as L^2 projection!
    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      auto res = integrator::template integrate_volUni_phivecfunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parametersT<space_dim, lSol_float_t>::initial_u, Point<hyEdge_dimT, lSol_float_t>>(i, hyper_edge.geometry, time);
      for (unsigned int dim = 0; dim < space_dim; dim++)
        u_old[i+dim*n_shape_fct_] = res[dim];
    }

    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      auto res = integrator::template integrate_volUni_phivecfunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parametersT<space_dim, lSol_float_t>::initial_v, Point<hyEdge_dimT, lSol_float_t>>(i, hyper_edge.geometry, time);
      for (unsigned int dim = 0; dim < space_dim; dim++)
        v_old[i+dim*n_shape_fct_] = res[dim];
    }

    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      auto res = integrator::template integrate_volUni_phivecfunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parametersT<space_dim, lSol_float_t>::initial_s, Point<hyEdge_dimT, lSol_float_t>>(i, hyper_edge.geometry, time);
      for (unsigned int dim = 0; dim < space_dim; dim++)
        s_old[i+dim*n_shape_fct_] = res[dim];
    }

    for (unsigned int i = 0; i < n_shape_fct_; ++i) {
      auto res = integrator::template integrate_volUni_phivecfunc<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parametersT<space_dim, lSol_float_t>::initial_r, Point<hyEdge_dimT, lSol_float_t>>(i, hyper_edge.geometry, time);
      for (unsigned int dim = 0; dim < space_dim; dim++)
        r_old[i+dim*n_shape_fct_] = res[dim];
    }

    // transform global dofs to edge dofs
    u_old = glob_dof_to_loc_dof(u_old, hyper_edge);
    r_old = glob_dof_to_loc_dof(r_old, hyper_edge);
    v_old = glob_dof_to_loc_dof(v_old, hyper_edge);
    s_old = glob_dof_to_loc_dof(s_old, hyper_edge);

    // n, m are algebraic: recover them from (u, r) and the boundary data of the static problem.
    SmallVec<n_loc_dofs_, lSol_float_t> rhs =
      assemble_rhs_from_lambda(lambda_values_loc, hyper_edge);
    add_dirichlet_rhs_static(rhs, hyper_edge, time);
    recover_dual(hyper_edge, rhs);

    return lambda_values;
  }

  template <class hyEdgeT, typename SmallMatT>
  SmallMatT& make_initial_from_static(SmallMatT& lambda_values_in,
                          hyEdgeT& hyper_edge,
                          const lSol_float_t time = 0.) const
  {
    // HACK
    // solves static problem
    // with rhs from lambdas and global rhs but where we consider all dirichlet

    auto lambda_values = node_dof_to_edge_dof(lambda_values_in, hyper_edge);

    // At bit-6 (static-only Dirichlet) faces the static problem assumes the trace
    // lambda is zero. The global loop's set_dof_values write-back from a previous
    // edge may have left a non-zero projection in x_vec at shared bit-6 junctions;
    // zero those entries so the rhs assembly stays invariant to edge ordering.
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      if (hyper_edge.node_descriptor[face] & (1u << 6))
        for (unsigned int d = 0; d < 2 * space_dim; ++d)
          lambda_values[face][d] = 0;

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    int comps[] = {1, -1, -2};
    static_assert(space_dim <= 3);

    // same massless_unloaded opt-in as assemble_rhs_stage: the static problem
    // this reconstructs must match the network solver's (unloaded weld) RHS
    bool loaded = true;
    if constexpr (requires { parameters::massless_unloaded; })
      if (parameters::massless_unloaded && hyper_edge.geometry.has_extra_data())
        loaded = hyper_edge.geometry.extra_data()[0] > 0;

    // static rhs = lambda traces + body loads + Dirichlet data (helpers shared with make_initial)
    SmallVec<n_loc_dofs_, lSol_float_t> rhs_full =
      assemble_rhs_from_lambda(lambda_values, hyper_edge);
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int c = 0; loaded && c < space_dim; c++) {
        rhs_full[(2 * space_dim + c)* n_shape_fct_ + i] +=
          integrator::template integrate_vol_phivecfunccomp<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::right_hand_side_n, Point<hyEdge_dimT, lSol_float_t>
          >(i, comps[c], hyper_edge.geometry, 0.);

        rhs_full[(3 * space_dim + c)* n_shape_fct_ + i] +=
          integrator::template integrate_vol_phivecfunccomp<
            Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
            parameters::right_hand_side_m, Point<hyEdge_dimT, lSol_float_t>
          >(i, comps[c], hyper_edge.geometry, 0.);
      }
    add_dirichlet_rhs_static(rhs_full, hyper_edge, 0.);

    auto mat = assemble_loc_matrix<false>(hyper_edge);
    SmallVec<4 * space_dim * n_shape_fct_, lSol_float_t> rhs, coeffs;
    for (unsigned int i = 0; i < 4 * space_dim * n_shape_fct_; ++i)
      rhs[i] = rhs_full[i];

    coeffs = rhs / mat;

    for (unsigned int i = 0; i < space_dim * n_shape_fct_; i++) {
      hyper_edge.data.n_old[i] = coeffs[0*space_dim*n_shape_fct_+i];
      hyper_edge.data.m_old[i] = coeffs[1*space_dim*n_shape_fct_+i];
      hyper_edge.data.u_old[i] = coeffs[2*space_dim*n_shape_fct_+i];
      hyper_edge.data.r_old[i] = coeffs[3*space_dim*n_shape_fct_+i];
    }

    // rest is zero initialized

    // at the static-only dirichlet nodes the trace lambda is zero
    // at these dirichlet nodes, the trace lambda should be non zero for the wave problem,
    // hence we need to project the static bulk solution to the static-only dirichlet
    // trace lambda

    // Project bulk u_old, r_old onto trace lambda on bit-6 faces,
    // only in components where bit-j (j=0..2*space_dim-1) is unset.
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if (!(hyper_edge.node_descriptor[face] & (1u << 6))) continue;

      for (unsigned int k = 0; k < n_shape_bdr_; ++k)
      {
        for (unsigned int d = 0; d < space_dim; ++d)
        {
            lSol_float_t num = 0;
            for (unsigned int i = 0; i < n_shape_fct_; ++i)
              num += integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                       i, k, face, hyper_edge.geometry)
                     * hyper_edge.data.u_old[d * n_shape_fct_ + i];
            lambda_values[face][k + d] = num;

            num = 0;
            for (unsigned int i = 0; i < n_shape_fct_; ++i)
              num += integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
                       i, k, face, hyper_edge.geometry)
                     * hyper_edge.data.r_old[d * n_shape_fct_ + i];
            lambda_values[face][k + space_dim + d] = num;
        }
      }
    }

    // write edge-local lambda back to global frame; edge_dof_to_node_dof accumulates,
    // so zero the destination first
    for (unsigned int i = 0; i < lambda_values_in.size(); ++i)
      for (unsigned int j = 0; j < lambda_values_in[i].size(); ++j)
        lambda_values_in[i][j] = 0.;
    edge_dof_to_node_dof(lambda_values, lambda_values_in, hyper_edge);

    return lambda_values_in;
  }
};  // end of class TimoshenkoWave


// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_lambda
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT,
          typename lSol_float_t,
          unsigned int n_stages>
template <typename hyEdgeT, typename SmallMatT>
inline SmallVec<
  TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t,
                 n_stages>::n_loc_dofs_,
  lSol_float_t>
TimoshenkoWave<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t, n_stages>::
  assemble_rhs_from_lambda(const SmallMatT& lambda_values, hyEdgeT& hyper_edge) const
{
  static_assert(std::is_same<typename SmallMatT::value_type::value_type, lSol_float_t>::value,
                "Lambda values should have same floating point arithmetics as local solver!");
  hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
            "The size of the lambda values should be twice the dimension of a hyperedge.");
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert(lambda_values[i].size() == 2 * space_dim * n_shape_bdr_,
              "The size of lambda should be the amount of ansatz functions at boundary.");

  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
    for (unsigned int j = 0; j < n_shape_bdr_; ++j)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        integral = integrator::template integrate_bdr_phipsi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
        {
          right_hand_side[(2 * space_dim + dim) * n_shape_fct_ + i] +=
            tau_ * lambda_values[face][j + dim] * integral;
          right_hand_side[dim * n_shape_fct_ + i] -=
            hyper_edge.geometry.local_normal(face).operator[](0) * lambda_values[face][j + dim] *
            integral;
        }
      }

  return right_hand_side;
}  // end of TimoshenkoWave::assemble_rhs_from_lambda


}  // namespace LocalSolver
