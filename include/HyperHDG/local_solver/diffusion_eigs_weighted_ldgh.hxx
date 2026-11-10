#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hypercube.hxx>
#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

#include <algorithm>
#include <array>
#include <tuple>

namespace LocalSolver
{
/*!*************************************************************************************************
 * \brief   Default parameters of the weighted diffusion eigenvalue problem, cf. below.
 *
 * Faces of node type 1 are Dirichlet faces, and faces of node type 3 carry hypernode mass. Every
 * other face that is not shared by several hyperedges (e.g. node type 2) is a natural boundary.
 *
 * \authors   Andreas Rupp, Saarland University, 2026.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct DiffusionEigsWeightedParametersDefault
{
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types corresponding to Dirichlet boundary.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 1U> dirichlet_nodes{1};
  /*!***********************************************************************************************
   * \brief   Array containing hypernode types of hypernodes that carry mass.
   ************************************************************************************************/
  static constexpr std::array<unsigned int, 1U> mass_nodes{3};
};  // end of struct DiffusionEigsWeightedParametersDefault

/*!*************************************************************************************************
 * \brief   LDG-H method for weighted diffusion eigenvalue problems with hypernode mass.
 *
 * The local solver discretizes the eigenvalue problem
 * \f[
 *  - \nabla \cdot ( d \kappa \nabla u ) + d c u = \lambda d \rho u \quad \text{ in } E,
 *  \qquad \sum_{E \ni N} d \kappa \nabla u \cdot \nu_E = \lambda \mu_N u \quad \text{ on } N,
 * \f]
 * on every hyperedge \f$E\f$ and every hypernode \f$N\f$ that carries mass, where \f$\nu_E\f$ is
 * the outward unit normal (conormal) of \f$E\f$, and the hypernodes of Dirichlet type have
 * homogeneous Dirichlet values. The plate weight \f$d\f$, the diffusion coefficient \f$\kappa\f$,
 * the mass density \f$\rho\f$, the reaction coefficient \f$c\f$, and the contribution \f$\mu\f$ of
 * the hyperedge to the mass density \f$\mu_N\f$ of its mass hypernodes are constant on every
 * hyperedge of the initial (unrefined) hypergraph. They are read from the hyperedge properties
 * (keyword HYPEREDGE_PROPERTIES of a .geo file) in the order \f$(d,\kappa,\rho,c,\mu)\f$; missing
 * entries default to \f$(1,1,1,0,0)\f$. Thus, \f$\mu_N\f$ is the sum of \f$\mu\f$ over all
 * hyperedges that contain \f$N\f$.
 *
 * In contrast to the other eigenvalue local solvers, the bulk unknowns are not eliminated. Only the
 * flux \f$q = - \kappa \nabla u\f$ is eliminated, i.e., for \f$(w, \hat w)\f$ in the space of
 * bulk and skeletal unknowns, \f$Q(w,\hat w)\f$ is the local solution of
 * \f[
 *  (\kappa^{-1} Q(w,\hat w), s)_E = (w, \nabla \cdot s)_E - \langle \hat w, s \cdot \nu \rangle_{
 *  \partial E}
 * \f]
 * for all discrete fluxes \f$s\f$. Then, the eigenvalue problem is the linear and generalized
 * eigenvalue problem \f$a_h(w, v) = \lambda m_h(w, v)\f$ with the symmetric bilinear forms
 * \f[
 *  a_h(w, v) = \sum_E d \left[ (\kappa^{-1} Q w, Q v)_E + c (w, v)_E + \tau \langle w - \hat w,
 *  v - \hat v \rangle_{\partial E} \right], \qquad m_h(w, v) = \sum_E d (\rho w, v)_E + \sum_N
 *  \mu_N \langle \hat w, \hat v \rangle_N.
 * \f]
 * This local solver provides the local matrices of these forms, which are assembled by
 * GlobalLoop::GeneralizedEigenvalue. The bulk unknowns and the flux are tensor-product polynomials
 * of degree \c poly_deg in each direction, and the skeletal unknowns are tensor-product
 * polynomials of the same degree on every face.
 *
 * \tparam  hyEdge_dimT   Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  poly_deg      The polynomial degree of test and trial functions.
 * \tparam  quad_deg      The order of the quadrature rule.
 * \tparam  parametersT   Struct depending on templates \c space_dimTP and \c lSol_float_TP that
 *                        contains the arrays \c dirichlet_nodes and (optionally) \c mass_nodes.
 * \tparam  lSol_float_t  The floating point type calculations are executed in. Defaults to double.
 *
 * \authors   Andreas Rupp, Saarland University, 2026.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT =
            DiffusionEigsWeightedParametersDefault,
          typename lSol_float_t = double>
class DiffusionEigsWeighted
{
 public:
  /*!***********************************************************************************************
   *  \brief  Define type of (hyperedge related) data that is stored in HyDataContainer.
   ************************************************************************************************/
  struct data_type
  {
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

  // -----------------------------------------------------------------------------------------------
  // Public, static constexpr functions
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Dimension of hyper edge type that this object solves on.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim() { return hyEdge_dimT; }
  /*!***********************************************************************************************
   * \brief   Evaluate amount of global degrees of freedom per hypernode.
   ************************************************************************************************/
  static constexpr unsigned int n_glob_dofs_per_node()
  {
    return Hypercube<hyEdge_dimT - 1>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief   Number of bulk degrees of freedom (of the primal unknown) per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_bulk_dofs()
  {
    return Hypercube<hyEdge_dimT>::pow(poly_deg + 1);
  }
  /*!***********************************************************************************************
   * \brief   Number of rows and columns of the local matrices.
   *
   * The local unknowns are ordered as follows: first the bulk unknowns, then the skeletal unknowns
   * of the faces in the order of the faces of the hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_dofs()
  {
    return n_loc_bulk_dofs() + 2 * hyEdge_dimT * n_glob_dofs_per_node();
  }
  /*!***********************************************************************************************
   * \brief   Dimension of the solution evaluated with respect to a hyperedge (flux and primal).
   ************************************************************************************************/
  static constexpr unsigned int system_dimension() { return hyEdge_dimT + 1; }
  /*!***********************************************************************************************
   * \brief   Dimension of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return 1; }
  /*!***********************************************************************************************
   * \brief   Number of hyperedge coefficients, i.e., weight, diffusion, density, reaction, mass.
   ************************************************************************************************/
  static constexpr unsigned int n_coefficients() { return 5U; }

 private:
  // -----------------------------------------------------------------------------------------------
  // Private, static constexpr functions and types
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Number of local shape functions (with respect to all spatial dimensions).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_fct_ = n_loc_bulk_dofs();
  /*!***********************************************************************************************
   * \brief   Number of local shape functions (with respect to a face / hypernode).
   ************************************************************************************************/
  static constexpr unsigned int n_shape_bdr_ = n_glob_dofs_per_node();
  /*!***********************************************************************************************
   * \brief   Number of local unknowns of the flux.
   ************************************************************************************************/
  static constexpr unsigned int n_flux_dofs_ = hyEdge_dimT * n_shape_fct_;
  /*!***********************************************************************************************
   * \brief   Number of rows and columns of the local matrices.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_dofs_ = n_loc_dofs();
  /*!***********************************************************************************************
   * \brief   Find out whether a node is of Dirichlet type.
   ************************************************************************************************/
  template <typename parameters>
  static constexpr bool is_dirichlet(const unsigned int node_type)
  {
    return std::find(parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                     node_type) != parameters::dirichlet_nodes.end();
  }
  /*!***********************************************************************************************
   * \brief   Find out whether a node carries mass.
   ************************************************************************************************/
  template <typename parameters>
  static constexpr bool is_mass_node(const unsigned int node_type)
  {
    if constexpr (requires { parameters::mass_nodes.begin(); })
      return std::find(parameters::mass_nodes.begin(), parameters::mass_nodes.end(), node_type) !=
             parameters::mass_nodes.end();
    else
      return false;
  }
  /*!***********************************************************************************************
   * \brief   The integrator (quadrature and shape functions of the bulk).
   ************************************************************************************************/
  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT> >,
    lSol_float_t>
    integrator;
  /*!***********************************************************************************************
   * \brief   The shape functions of the skeletal unknowns on a face.
   ************************************************************************************************/
  typedef TPP::ShapeFunction<
    TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT - 1> >
    face_shape_t;

  // -----------------------------------------------------------------------------------------------
  // Private, const members
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_;

  // -----------------------------------------------------------------------------------------------
  // Private, internal functions for the local solver
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Assemble the matrices that define the local flux reconstruction.
   *
   * The coefficients of \f$Q(w, \hat w)\f$ are \f$M_q^{-1} G (w, \hat w)\f$, where \f$M_q\f$ is
   * the (weighted) mass matrix of the flux, and \f$G\f$ is returned by this function.
   *
   * \param   hyper_edge    The hyperedge.
   * \param   kappa         The diffusion coefficient of the hyperedge.
   * \param   mass_flux     The matrix \f$M_q\f$ (output).
   * \retval  G             The matrix \f$G\f$.
   ************************************************************************************************/
  template <typename hyEdgeT>
  SmallMat<n_flux_dofs_, n_loc_dofs_, lSol_float_t> flux_matrices(
    hyEdgeT& hyper_edge,
    const lSol_float_t kappa,
    SmallSquareMat<n_flux_dofs_, lSol_float_t>& mass_flux) const
  {
    using geom_t = decltype(hyEdgeT::geometry);
    SmallMat<n_flux_dofs_, n_loc_dofs_, lSol_float_t> G;
    mass_flux = SmallSquareMat<n_flux_dofs_, lSol_float_t>();
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
      {
        const lSol_float_t mass =
          integrator::template integrate_vol_phiphi<geom_t>(i, j, hyper_edge.geometry);
        // (w, div s) with w = phi_i and s = phi_j e_dim: integral of d_dim phi_j * phi_i
        const SmallVec<hyEdge_dimT, lSol_float_t> grad =
          integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dimT, lSol_float_t>,
                                                         geom_t>(j, i, hyper_edge.geometry);
        for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
        {
          mass_flux(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) = mass / kappa;
          G(dim * n_shape_fct_ + j, i) += grad[dim];
        }
      }
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      const auto normal = hyper_edge.geometry.local_normal(face);  // outward unit normal
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
        for (unsigned int l = 0; l < n_shape_bdr_; ++l)
        {
          const lSol_float_t integral =
            integrator::template integrate_bdr_phipsi<geom_t>(j, l, face, hyper_edge.geometry);
          for (unsigned int dim = 0; dim < hyEdge_dimT; ++dim)
            G(dim * n_shape_fct_ + j, n_shape_fct_ + face * n_shape_bdr_ + l) -=
              normal[dim] * integral;
        }
    }
    return G;
  }

 public:
  // -----------------------------------------------------------------------------------------------
  // Public functions (and one typedef) to be utilized by external functions.
  // -----------------------------------------------------------------------------------------------

  /*!***********************************************************************************************
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef lSol_float_t constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   tau           Penalty parameter of HDG scheme.
   ************************************************************************************************/
  DiffusionEigsWeighted(const constructor_value_type& tau = 1.) : tau_(tau) {}
  /*!***********************************************************************************************
   * \brief   Coefficients \f$(d,\kappa,\rho,c,\mu)\f$ of a hyperedge, cf. class description.
   ************************************************************************************************/
  template <typename hyEdgeT>
  static std::array<lSol_float_t, 5U> coefficients(hyEdgeT& hyper_edge)
  {
    std::array<lSol_float_t, 5U> coeffs{1., 1., 1., 0., 0.};
    if constexpr (requires {
                    hyper_edge.geometry.has_extra_data();
                    hyper_edge.geometry.extra_data();
                  })
      if (hyper_edge.geometry.has_extra_data())
      {
        const auto& properties = hyper_edge.geometry.extra_data();
        for (unsigned int i = 0; i < std::min<std::size_t>(properties.size(), coeffs.size()); ++i)
          coeffs[i] = properties[i];
      }
    return coeffs;
  }
  /*!***********************************************************************************************
   * \brief   Local matrix of the bilinear form \f$a_h\f$ (stiffness matrix).
   *
   * \param   hyper_edge    The hyperedge.
   * \retval  stiff_mat     The local matrix with respect to bulk and skeletal unknowns.
   ************************************************************************************************/
  template <typename hyEdgeT>
  SmallSquareMat<n_loc_dofs_, lSol_float_t> stiffness_matrix(hyEdgeT& hyper_edge) const
  {
    using geom_t = decltype(hyEdgeT::geometry);
    const std::array<lSol_float_t, 5U> coeffs = coefficients(hyper_edge);
    const lSol_float_t weight = coeffs[0], kappa = coeffs[1], reaction = coeffs[3];

    // Energy of the reconstructed flux: G^T M_q^{-1} G.
    SmallSquareMat<n_flux_dofs_, lSol_float_t> mass_flux;
    const SmallMat<n_flux_dofs_, n_loc_dofs_, lSol_float_t> G =
      flux_matrices(hyper_edge, kappa, mass_flux);
    SmallMat<n_flux_dofs_, n_loc_dofs_, lSol_float_t> rhs(G);
    const SmallMat<n_flux_dofs_, n_loc_dofs_, lSol_float_t> minv_G = rhs / mass_flux;
    SmallSquareMat<n_loc_dofs_, lSol_float_t> stiff_mat = transposed_mat_times_mat(G, minv_G);

    // Reaction term.
    if (reaction != 0.)
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
          stiff_mat(i, j) +=
            reaction * integrator::template integrate_vol_phiphi<geom_t>(i, j, hyper_edge.geometry);

    // Stabilization tau < w - w_hat, v - v_hat > on all faces.
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      const unsigned int offset = n_shape_fct_ + face * n_shape_bdr_;
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
      {
        for (unsigned int j = 0; j < n_shape_fct_; ++j)
          stiff_mat(i, j) += tau_ * integrator::template integrate_bdr_phiphi<geom_t>(
                                      i, j, face, hyper_edge.geometry);
        for (unsigned int l = 0; l < n_shape_bdr_; ++l)
        {
          const lSol_float_t integral = tau_ * integrator::template integrate_bdr_phipsi<geom_t>(
                                                 i, l, face, hyper_edge.geometry);
          stiff_mat(i, offset + l) -= integral;
          stiff_mat(offset + l, i) -= integral;
        }
      }
      for (unsigned int l = 0; l < n_shape_bdr_; ++l)
        for (unsigned int k = 0; k < n_shape_bdr_; ++k)
          stiff_mat(offset + l, offset + k) +=
            tau_ *
            integrator::template integrate_bdr_psipsi<geom_t>(l, k, face, hyper_edge.geometry);
    }

    return weight * stiff_mat;
  }
  /*!***********************************************************************************************
   * \brief   Local matrix of the bilinear form \f$m_h\f$ (mass matrix).
   *
   * \param   hyper_edge    The hyperedge.
   * \retval  mass_mat      The local matrix with respect to bulk and skeletal unknowns.
   ************************************************************************************************/
  template <typename hyEdgeT>
  SmallSquareMat<n_loc_dofs_, lSol_float_t> mass_matrix(hyEdgeT& hyper_edge) const
  {
    using geom_t = decltype(hyEdgeT::geometry);
    using parameters = parametersT<geom_t::space_dim(), lSol_float_t>;
    const std::array<lSol_float_t, 5U> coeffs = coefficients(hyper_edge);
    const lSol_float_t bulk_density = coeffs[0] * coeffs[2], node_mass = coeffs[4];

    SmallSquareMat<n_loc_dofs_, lSol_float_t> mass_mat;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
        mass_mat(i, j) = bulk_density * integrator::template integrate_vol_phiphi<geom_t>(
                                          i, j, hyper_edge.geometry);

    if (node_mass != 0.)
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        if (is_mass_node<parameters>(hyper_edge.node_descriptor[face]))
        {
          const unsigned int offset = n_shape_fct_ + face * n_shape_bdr_;
          for (unsigned int l = 0; l < n_shape_bdr_; ++l)
            for (unsigned int k = 0; k < n_shape_bdr_; ++k)
              mass_mat(offset + l, offset + k) +=
                node_mass *
                integrator::template integrate_bdr_psipsi<geom_t>(l, k, face, hyper_edge.geometry);
        }

    return mass_mat;
  }
  /*!***********************************************************************************************
   * \brief   Fill an array with 1 if the node is Dirichlet and 0 otherwise.
   ************************************************************************************************/
  template <typename hyEdgeT>
  std::array<unsigned int, 2 * hyEdge_dimT> node_types(hyEdgeT& hyper_edge) const
  {
    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    std::array<unsigned int, 2 * hyEdge_dimT> result;
    result.fill(0);
    for (unsigned int i = 0; i < result.size(); ++i)
      result[i] = is_dirichlet<parameters>(hyper_edge.node_descriptor[i]) ? 1 : 0;
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate a bulk function (given by its coefficients) at a reference point.
   ************************************************************************************************/
  template <typename coeffsT, typename pointT>
  static lSol_float_t bulk_value(const coeffsT& coeffs, const pointT& ref_point)
  {
    return integrator::shape_fun_t::template lin_comb_fct_val<lSol_float_t>(coeffs, ref_point);
  }
  /*!***********************************************************************************************
   * \brief   Evaluate a skeletal function (given by its coefficients) at a reference point.
   *
   * The reference point is a point of the reference face, whose coordinates are those of the
   * reference hyperedge without the coordinate that is constant on the face.
   ************************************************************************************************/
  template <typename coeffsT, typename pointT>
  static lSol_float_t face_value(const coeffsT& coeffs, const pointT& ref_point)
  {
    return face_shape_t::template lin_comb_fct_val<lSol_float_t>(coeffs, ref_point);
  }
  /*!***********************************************************************************************
   * \brief   Evaluate the local reconstruction of an eigenfunction at tensorial abscissas.
   *
   * This function is used for plotting: The skeletal values of an eigenfunction and the eigenvalue
   * determine its bulk values by the local equations \f$(A_{uu} - \lambda M_{uu}) u = - A_{u\hat u}
   * \hat u\f$ (if the eigenvalue is not an eigenvalue of the local problem). The flux is
   * reconstructed from the bulk and skeletal values.
   *
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The hyperedge.
   * \param   eigenvalue        The eigenvalue (passed as time by the plot functions).
   * \retval  function_values   Function values (flux components, then primal) at the points.
   ************************************************************************************************/
  template <typename abscissa_float_t,
            std::size_t abscissas_sizeT,
            class input_array_t,
            class hyEdgeT>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>,
             hyEdge_dimT + 1>
  bulk_values(const std::array<abscissa_float_t, abscissas_sizeT>& abscissas,
              const input_array_t& lambda_values,
              hyEdgeT& hyper_edge,
              const lSol_float_t eigenvalue = 0.) const
  {
    const SmallSquareMat<n_loc_dofs_, lSol_float_t> stiff_mat = stiffness_matrix(hyper_edge),
                                                    mass_mat = mass_matrix(hyper_edge);
    SmallVec<n_loc_dofs_, lSol_float_t> loc_dofs;
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      for (unsigned int l = 0; l < n_shape_bdr_; ++l)
        loc_dofs[n_shape_fct_ + face * n_shape_bdr_ + l] = lambda_values[face][l];

    SmallSquareMat<n_shape_fct_, lSol_float_t> loc_mat;
    SmallVec<n_shape_fct_, lSol_float_t> rhs;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
    {
      for (unsigned int j = 0; j < n_shape_fct_; ++j)
        loc_mat(i, j) = stiff_mat(i, j) - eigenvalue * mass_mat(i, j);
      for (unsigned int j = n_shape_fct_; j < n_loc_dofs_; ++j)
        rhs[i] -= stiff_mat(i, j) * loc_dofs[j];
    }
    const SmallVec<n_shape_fct_, lSol_float_t> bulk = rhs / loc_mat;
    for (unsigned int i = 0; i < n_shape_fct_; ++i)
      loc_dofs[i] = bulk[i];

    SmallSquareMat<n_flux_dofs_, lSol_float_t> mass_flux;
    const SmallMat<n_flux_dofs_, n_loc_dofs_, lSol_float_t> G =
      flux_matrices(hyper_edge, coefficients(hyper_edge)[1], mass_flux);
    const SmallVec<n_flux_dofs_, lSol_float_t> flux = (G * loc_dofs) / mass_flux;

    SmallVec<static_cast<unsigned int>(abscissas_sizeT), abscissa_float_t> helper(abscissas);
    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(abscissas_sizeT)>,
               hyEdge_dimT + 1>
      point_vals;
    SmallVec<n_shape_fct_, lSol_float_t> coeffs;
    for (unsigned int d = 0; d < hyEdge_dimT + 1; ++d)
    {
      for (unsigned int i = 0; i < n_shape_fct_; ++i)
        coeffs[i] = (d < hyEdge_dimT) ? flux[d * n_shape_fct_ + i] : loc_dofs[i];
      for (unsigned int pt = 0; pt < Hypercube<hyEdge_dimT>::pow(abscissas_sizeT); ++pt)
        point_vals[d][pt] = bulk_value(
          coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT, lSol_float_t> >(
                    pt, helper));
    }
    return point_vals;
  }
};  // end of class DiffusionEigsWeighted

}  // namespace LocalSolver
