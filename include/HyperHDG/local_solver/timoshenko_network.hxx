#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/global_loop/kinds.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/local_solver/local_operator.hxx>
#include <tpp/quadrature/tensorial.hxx>
#include <tpp/shape_function/shape_function.hxx>

#include <array>
#include <cmath>
#include <limits>
#include <tuple>

namespace LocalSolver
{

/*!*************************************************************************************************
 * \brief   Default parameters of the Timoshenko beam solver: a uniform force -e_z on every beam, no
 *          distributed moment, zero displacement and rotation at all Dirichlet nodes.
 *
 * There is no closed-form solution to this problem, so analytic_result_u and analytic_result_phi
 * are zero and errors() returns the norm of the discrete solution. Convergence studies use their
 * own parameters, cf. tests/parameters/timoshenko.hxx. The functions receive a
 * point and a direction and return the component of the vector field along that direction.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TimoschenkoBeamParametersDefault
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
   * \brief   Distributed force: uniform, -e_z.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>&,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return -normal[2];
  }
  /*!***********************************************************************************************
   * \brief   Distributed moment: none.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>&,
                                         const Point<space_dimT, param_float_t>&,
                                         const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of the displacement: zero.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>&,
                                         const Point<space_dimT, param_float_t>&,
                                         const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of the rotation: zero.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>&,
                                           const Point<space_dimT, param_float_t>&,
                                           const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   No closed-form displacement: zero, so errors() is the norm of the discrete solution.
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>&,
                                         const Point<space_dimT, param_float_t>&,
                                         const param_float_t = 0.)
  {
    return 0.;
  }
  /*!***********************************************************************************************
   * \brief   No closed-form rotation: zero.
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>&,
                                           const Point<space_dimT, param_float_t>&,
                                           const param_float_t = 0.)
  {
    return 0.;
  }
};  // end of struct TimoschenkoBeamParametersDefault

/*!*************************************************************************************************
 * \brief   Default parameters for the diffusion equation, cf. below.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct TimoschenkoBeamParametersClamped
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
   * \brief   Inverse diffusionbeam_network_bilaplacian.hxx coefficient in PDE as analytic function.
   ************************************************************************************************/
  // static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
  //                                              const param_float_t = 0.)
  // {
  //   return 1.;
  //   // return 1. / M_PI / M_PI;
  // }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_n(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return 0;
    // return -M_PI * cos(M_PI * point[0]) * normal[2] * (point[1] == 0. && point[2] == 0.);
    // return M_PI * M_PI * cos(M_PI * point[0]) * normal[2] * (point[1] == 0. && point[2] == 0.);
    // return M_PI * M_PI * sin(M_PI * point[0]) * normal[0];
    // return M_PI * M_PI * sin(M_PI * point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Right-hand side in PDE as analytic function.
   ************************************************************************************************/
  static param_float_t right_hand_side_m(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return 0.;
    // return  (M_PI * M_PI + 1.) * sin(M_PI * point[0]) * normal[1] * (point[1] == 0. && point[2]
    // == 0.); return  -M_PI * sin(M_PI * point[0]) * normal[1] * (point[1] == 0. && point[2] ==
    // 0.); return M_PI * M_PI * sin(M_PI * point[0]) * normal[0]; return M_PI * M_PI * sin(M_PI *
    // point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    return analytic_result_u(point, normal);
  }
  /*!***********************************************************************************************
   * \brief   Dirichlet values of solution as analytic function.
   ************************************************************************************************/
  static param_float_t dirichlet_value_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t = 0.)
  {
    return analytic_result_phi(point, normal);
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_u(const Point<space_dimT, param_float_t>& point,
                                         const Point<space_dimT, param_float_t>& normal,
                                         const param_float_t = 0.)
  {
    // return 0.;
    return 5e-4 * ((point[0] > 1e-3) * normal[0] + (point[1] > .5e-3) * normal[1] + normal[2]);
    // return point[0] * normal[0];
    // return sin(M_PI * point[0]) * normal[0];
  }
  /*!***********************************************************************************************
   * \brief   Analytic result of PDE (for convergence tests).
   ************************************************************************************************/
  static param_float_t analytic_result_phi(const Point<space_dimT, param_float_t>& point,
                                           const Point<space_dimT, param_float_t>& normal,
                                           const param_float_t = 0.)
  {
    return 0.;
    // return point[0] * normal[0];
    // return sin(M_PI * point[0]) * normal[0];
  }
};  // end of struct TimoschenkoBeamParametersClamped

/*!*************************************************************************************************
 * \brief   Layout of the Timoshenko beam: the frame of the beam and its six stiffnesses.
 *
 * In the frame of the beam (tangent, two normals) the local problem falls into four independent
 * ones: stretching (u_0), twisting (r_0), and bending in the two planes, which couples the
 * transverse displacement with the rotation about the other normal, (u_2, r_1) and (u_1, r_2).
 * Without a load, the fluxes of each are in balance and vanish for rigid motions, so the condensed
 * operator of a problem is fixed by few numbers: with the differences d(.) between the two nodes,
 *   stretching:  flux = -+ k_a d(u_0),            twisting likewise with k_t and r_0,
 *   bending:     shear force  = -+ k_s e_s,       e_s = d(u) - h (r at node 0 + r at node 1),
 *                moment       = -+ k_b d(r) + h k_s e_s,
 * (upper sign at node 1) with h = +- half the length of the beam. The operator starts with the
 * frame, h and the six stiffnesses, which the solver computes from the four small local problems
 * (TimoshenkoBeam::beam_operator): 16 numbers, two cache lines in double, whatever the degree,
 * and all that apply reads. apply rotates the traces into the frame, applies the formulas and
 * rotates back.
 *
 * The residual also needs the fluxes of the local load, K l with K = C A^-1. K falls into the same
 * four parts, 40 n numbers for n shape functions, which the operator keeps behind the 16 (about a
 * quarter of the size of Condensed in total), so that the residual factorizes nothing.
 * local_solution, for errors, solves the four problems of its beam on the fly.
 **************************************************************************************************/
struct BeamStiffness
{
  // The components (in the frame of the beam) of the four parts: stretching and twisting, and the
  // force and the moment components of the two bending problems.
  static constexpr std::array<unsigned int, 2> single_dim{0, 3};
  static constexpr std::array<unsigned int, 2> force_dim{2, 1}, moment_dim{4, 5};

  template <typename T, unsigned int n>
  struct alignas(64) Operator
  {
    using Trace = std::array<std::array<T, 6>, 2>;  // per node: u and r in the frame of the beam

    // ---- read by apply
    std::array<T, 9> frame;  // rows: the tangent and the two normals of the beam
    T half_length;           // h for (u_2, r_1), -h for (u_1, r_2); signed with the tangent
    T axial, torsion;
    std::array<T, 2> shear, bending;  // of (u_2, r_1) and of (u_1, r_2)
    // ---- read by residual only: K of the four parts. Rows: the flux at node 0 and node 1 (for
    // bending: of the force, then of the moment); columns: the unknowns of the part, as
    // TimoshenkoBeam::BeamProblem orders them.
    std::array<SmallMat<2, 2 * n, T>, 2> load_to_flux_single;
    std::array<SmallMat<4, 4 * n, T>, 2> load_to_flux_coupled;

    template <unsigned int k, class SolverT, class Edge, class Args>
    void apply(const SolverT&,
               Edge&,
               const SmallMat<12, k, T>& in,
               const Args&,
               SmallMat<12, k, T>& out) const
    {
      for (unsigned int col = 0; col < k; ++col)
        to_nodes(flux(to_beam(in, col)), out, col);
    }
    template <class SolverT, class Edge, class Args>
    void residual(const SolverT& solver,
                  Edge& edge,
                  const SmallVec<12, T>& lambda,
                  const Args& args,
                  SmallVec<12, T>& r) const
    {
      const auto load = solver.local_load(edge, args);
      Trace result = flux(to_beam(lambda, 0));
      for (unsigned int c = 0; c < 2; ++c)
      {
        const std::array<unsigned int, 2> first{single_dim[c] * n, (6 + single_dim[c]) * n};
        for (unsigned int row = 0; row < 2; ++row)
          result[row][single_dim[c]] += flux_of_load(load_to_flux_single[c], row, first, load);
        const std::array<unsigned int, 4> first_coupled{
          force_dim[c] * n, moment_dim[c] * n, (6 + force_dim[c]) * n, (6 + moment_dim[c]) * n};
        for (unsigned int row = 0; row < 4; ++row)
          result[row % 2][row < 2 ? force_dim[c] : moment_dim[c]] +=
            flux_of_load(load_to_flux_coupled[c], row, first_coupled, load);
      }
      to_nodes(result, r, 0);
    }
    template <class SolverT, class Edge, class Args>
    auto local_solution(const SolverT& solver,
                        Edge& edge,
                        const SmallVec<12, T>& lambda,
                        const Args& args) const
    {
      return solved(solver, edge, solver.beam_problem(edge), to_beam(lambda, 0), args);
    }

   private:
    // Row row of K of one part times the load of its unknowns: blocks of n starting at first[.].
    template <unsigned int n_rows, std::size_t n_blocks, class LoadT>
    static T flux_of_load(const SmallMat<n_rows, n_blocks * n, T>& load_to_flux,
                          const unsigned int row,
                          const std::array<unsigned int, n_blocks>& first,
                          const LoadT& load)
    {
      T value = 0.;
      for (unsigned int block = 0; block < n_blocks; ++block)
        for (unsigned int i = 0; i < n; ++i)
          value += load_to_flux(row, block * n + i) * load[first[block] + i];
      return value;
    }
    // The local solution of the beam for the trace (in its frame) and the local load.
    template <class SolverT, class Edge, class ProblemT, class Args>
    static auto solved(const SolverT& solver,
                       Edge& edge,
                       const ProblemT& problem,
                       const Trace& trace,
                       const Args& args)
    {
      auto x = problem.rhs(trace);
      const auto load = solver.local_load(edge, args);
      for (unsigned int i = 0; i < x.size(); ++i)
        x[i] += load[i];
      problem.solve(x);
      return x;
    }
    // Column col of the traces at the nodes (u and r in space), rotated into the frame.
    template <unsigned int k>
    Trace to_beam(const SmallMat<12, k, T>& values, const unsigned int col) const
    {
      Trace trace;
      for (unsigned int block = 0; block < 4; ++block)  // u and r at node 0, u and r at node 1
        for (unsigned int c = 0; c < 3; ++c)
        {
          T value = 0.;
          for (unsigned int d = 0; d < 3; ++d)
            value += frame[3 * c + d] * values(3 * block + d, col);
          trace[block / 2][3 * (block % 2) + c] = value;
        }
      return trace;
    }
    template <unsigned int k>
    void to_nodes(const Trace& trace, SmallMat<12, k, T>& values, const unsigned int col) const
    {
      for (unsigned int block = 0; block < 4; ++block)
        for (unsigned int d = 0; d < 3; ++d)
        {
          T value = 0.;
          for (unsigned int c = 0; c < 3; ++c)
            value += frame[3 * c + d] * trace[block / 2][3 * (block % 2) + c];
          values(3 * block + d, col) = value;
        }
    }
    Trace flux(const Trace& trace) const
    {
      Trace result;
      const T stretch = axial * (trace[1][0] - trace[0][0]);
      result[0][0] = stretch, result[1][0] = -stretch;
      const T twist = torsion * (trace[1][3] - trace[0][3]);
      result[0][3] = twist, result[1][3] = -twist;
      bend(trace, result, 2, 4, half_length, shear[0], bending[0]);
      bend(trace, result, 1, 5, -half_length, shear[1], bending[1]);
      return result;
    }
    static void bend(const Trace& trace,
                     Trace& result,
                     const unsigned int u,
                     const unsigned int r,
                     const T h,
                     const T k_s,
                     const T k_b)
    {
      const T force = k_s * (trace[1][u] - trace[0][u] - h * (trace[0][r] + trace[1][r]));
      const T moment = k_b * (trace[1][r] - trace[0][r]);
      result[0][u] = force, result[1][u] = -force;
      result[0][r] = moment + h * force, result[1][r] = -moment + h * force;
    }
  };

  template <class SolverT, class Edge, class Kind, class Args>
  static auto make(const SolverT& solver, Edge& edge, Kind, const Args&)
  {
    return solver.beam_operator(edge);
  }
};

/*!*************************************************************************************************
 * \brief   Local solver for the equation that governs the bending and change of length of an
 *          elastic Bernoulli beam.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename> typename parametersT = TimoschenkoBeamParametersDefault,
          typename lSol_float_t = double,
          class LayoutT = BeamStiffness>
class TimoshenkoBeam
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
  static constexpr unsigned int system_dimension() { return space_dim; }
  /*!***********************************************************************************************
   * \brief   Dimension of of the solution evaluated with respect to a hypernode.
   ************************************************************************************************/
  static constexpr unsigned int node_system_dimension() { return space_dim; }

  template <typename parameters>
  static constexpr bool is_dirichlet(const unsigned int node_type)
  {
    return std::find(parameters::dirichlet_nodes.begin(), parameters::dirichlet_nodes.end(),
                     node_type) != parameters::dirichlet_nodes.end();
  }

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
  static constexpr unsigned int n_loc_dofs_ = 4 * space_dim * n_shape_fct_;
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
   * \brief   (Globally constant) penalty parameter for HDG scheme.
   ************************************************************************************************/
  const lSol_float_t tau_;

  typedef TPP::Quadrature::Tensorial<
    TPP::Quadrature::GaussLegendre<quad_deg>,
    TPP::ShapeFunction<TPP::ShapeType::Tensorial<TPP::ShapeType::Legendre<poly_deg>, hyEdge_dimT>>,
    lSol_float_t>
    integrator;

 public:
  /*!***********************************************************************************************
   * \brief   Class is constructed using a single double indicating the penalty parameter.
   ************************************************************************************************/
  typedef lSol_float_t constructor_value_type;
  /*!***********************************************************************************************
   * \brief   Constructor for local solver.
   *
   * \param   tau           Penalty parameter of HDG scheme.
   ************************************************************************************************/
  TimoshenkoBeam(const constructor_value_type& tau = 1.) : tau_(tau) {}

  /*!***********************************************************************************************
   * \brief   The six stiffness coefficients of a beam in its frame (EA, the two shear stiffnesses,
   *          the torsional and the two bending stiffnesses); one without data of the geometry.
   ************************************************************************************************/
  template <typename hyEdgeT>
  SmallVec<2 * space_dim, lSol_float_t> beam_coefficients(hyEdgeT& hyper_edge) const
  {
    SmallVec<2 * space_dim, lSol_float_t> extra_coeffs(1.);

    if (hyper_edge.geometry.has_extra_data())
    {
      auto extra_data = hyper_edge.geometry.extra_data();
      // (EA, kG_1A, kG_2A, G_xI_x, E_1I_1, E_2I_2):   6 structural constants
      // (n_11,n_12,n_13) : normal 1
      // (n_21,n_22,n_23) : normal 2
      SmallVec<space_dim, lSol_float_t> normal1 =
        std::array<lSol_float_t, space_dim>{{extra_data[6], extra_data[7], extra_data[8]}};
      SmallVec<space_dim, lSol_float_t> normal2 =
        std::array<lSol_float_t, space_dim>{{extra_data[9], extra_data[10], extra_data[11]}};
      SmallVec<space_dim, lSol_float_t> outer1 = hyper_edge.geometry.outer_normal(0);
      SmallVec<space_dim, lSol_float_t> outer2 = hyper_edge.geometry.outer_normal(1);

      extra_coeffs[0] = extra_data[0];
      extra_coeffs[1] = extra_data[1] * scalar_product(outer1, normal1) +
                        extra_data[2] * scalar_product(outer1, normal2);
      extra_coeffs[2] = extra_data[1] * scalar_product(outer2, normal1) +
                        extra_data[2] * scalar_product(outer2, normal2);
      extra_coeffs[3] = extra_data[3];
      extra_coeffs[4] = extra_data[4] * scalar_product(outer1, normal1) +
                        extra_data[5] * scalar_product(outer1, normal2);
      extra_coeffs[5] = extra_data[4] * scalar_product(outer2, normal1) +
                        extra_data[5] * scalar_product(outer2, normal2);

      extra_coeffs[0] *= 1e4;
      extra_coeffs[1] *= 1e4;
      extra_coeffs[2] *= 1e4;
      extra_coeffs[3] *= 1e12;
      extra_coeffs[4] *= 1e12;
      extra_coeffs[5] *= 1e12;
    }

    for (unsigned int i = 0; i < extra_coeffs.size(); ++i)
      extra_coeffs[i] = std::abs(extra_coeffs[i]);
    // extra_coeffs[i] = 1.;
    // hy_assert(false, extra_coeffs);
    return extra_coeffs;
  }

  template <typename hyEdgeT>
  inline SmallSquareMat<n_loc_dofs_, lSol_float_t> assemble_loc_matrix(
    hyEdgeT& hyper_edge,
    const lSol_float_t time) const;

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_lambda(
    const SmallMatT& lambda_values,
    hyEdgeT& hyper_edge) const;

  template <typename hyEdgeT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> assemble_rhs_from_global_rhs(
    hyEdgeT& hyper_edge,
    const unsigned int dim) const;

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

  template <typename hyEdgeT, typename SmallMatT>
  inline SmallVec<n_loc_dofs_, lSol_float_t> solve_local_problem(const SmallMatT& lambda_values,
                                                                 const unsigned int solution_type,
                                                                 hyEdgeT& hyper_edge,
                                                                 const lSol_float_t time) const
  {
    try
    {
      SmallVec<n_loc_dofs_, lSol_float_t> rhs;
      if (solution_type == 0)
        rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge);
      else if (solution_type == 1)
        rhs = assemble_rhs_from_lambda(lambda_values, hyper_edge) +
              assemble_rhs_from_global_rhs(hyper_edge, time);
      else
        hy_assert(0 == 1, "This has not been implemented!");
      return rhs / assemble_loc_matrix(hyper_edge, time);
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

  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& trace_to_flux(const SmallMatInT& lambda_values_in,
                              SmallMatOutT& lambda_values_out,
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

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    SmallMatInT lambda_values_loc = node_dof_to_edge_dof(lambda_values_in, hyper_edge);

    // for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    //   for (unsigned int j = 0; j < 2 * space_dim; ++j)
    //     std::cout << lambda_values_loc[i][j] << " ";
    // std::cout << std::endl;

    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values_loc, 0U, hyper_edge, time);

    auto result = extract_fluxes_from_coeffs(coeffs, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < 2 * space_dim; ++j)
        lambda_values_loc[i][j] = result(i, j) - tau_ * lambda_values_loc[i][j];

    // for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    //   for (unsigned int j = 0; j < 2 * space_dim; ++j)
    //     std::cout << result(i,j) << " ";
    // std::cout << std::endl;

    // for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    //   for (unsigned int j = 0; j < 2 * space_dim; ++j)
    //     std::cout << lambda_values_loc[i][j] << " ";
    // std::cout << std::endl << std::endl;

    lambda_values_out = edge_dof_to_node_dof(lambda_values_loc, lambda_values_out, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        lambda_values_out[i].fill(0.);

    return lambda_values_out;
  }

  template <typename hyEdgeT, typename SmallMatInT, typename SmallMatOutT>
  SmallMatOutT& residual_flux(const SmallMatInT& lambda_values_in,
                              SmallMatOutT& lambda_values_out,
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

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;

    SmallMatInT lambda_values_loc = node_dof_to_edge_dof(lambda_values_in, hyper_edge);

    SmallVec<n_loc_dofs_, lSol_float_t> coeffs =
      solve_local_problem(lambda_values_loc, 1U, hyper_edge, time);

    auto result = extract_fluxes_from_coeffs(coeffs, hyper_edge);
    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      for (unsigned int j = 0; j < 2 * space_dim; ++j)
        lambda_values_loc[i][j] = result(i, j) - tau_ * lambda_values_loc[i][j];
    lambda_values_out = edge_dof_to_node_dof(lambda_values_loc, lambda_values_out, hyper_edge);

    for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[i]))
        lambda_values_out[i].fill(0.);

    return lambda_values_out;
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
  std::array<lSol_float_t, 1U> errors(
    const std::array<std::array<lSol_float_t, n_glob_dofs_per_node()>, 2 * hyEdge_dimT>&
      lambda_values,
    hyEdgeT& hyper_edge,
    const lSol_float_t time = 0.) const
  {
    hy_assert(lambda_values.size() == 2 * hyEdge_dimT, "Matrix must have appropriate size!");
    for (unsigned int i = 0; i < lambda_values.size(); ++i)
      hy_assert(lambda_values[i].size() == n_glob_dofs_per_node(),
                "Matrix must have appropriate size!");

    using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    auto lambda_values_loc = node_dof_to_edge_dof(lambda_values, hyper_edge);

    SmallVec<n_loc_dofs_, lSol_float_t> coefficients =
      solve_local_problem(lambda_values_loc, 1U, hyper_edge, time);
    lSol_float_t error = 0.;
    std::array<lSol_float_t, n_shape_fct_> coeffs;

    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[i + (2 * space_dim + 0) * n_shape_fct_];
    error += integrator::template integrate_vol_diffsquare_discanacomp<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, 1,
                                                                       hyper_edge.geometry, time);

    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[i + (2 * space_dim + 1) * n_shape_fct_];
    error += integrator::template integrate_vol_diffsquare_discanacomp<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, -1,
                                                                       hyper_edge.geometry, time);

    for (unsigned int i = 0; i < coeffs.size(); ++i)
      coeffs[i] = coefficients[i + (2 * space_dim + 2) * n_shape_fct_];
    error += integrator::template integrate_vol_diffsquare_discanacomp<
      Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
      parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, -2,
                                                                       hyper_edge.geometry, time);

    // error += len_beam.errors(lambda_values, hyper_edge, time)[0];
    // error += twist_beam.errors(lambda_values, hyper_edge, time)[0];
    // error += ben_beam.errors(lambda_values, hyper_edge, time)[0];
    return std::array<lSol_float_t, 1U>({error});
  }
  /*!***********************************************************************************************
   * \brief   Evaluate local local reconstruction at tensorial products of abscissas.
   *
   * \tparam  absc_float_t      Floating type for the abscissa values.
   * \tparam  abscissas_sizeT   Size of the array of array of abscissas.
   * \tparam  input_array_t     Type of input array.
   * \tparam  hyEdgeT           The geometry type / typename of the considered hyEdge's geometry.
   * \param   abscissas         Abscissas of the supporting points.
   * \param   lambda_values     The values of the skeletal variable's coefficients.
   * \param   hyper_edge        The geometry of the considered hyperedge (of typename GeomT).
   * \param   time              Time.
   * \retval  func_values       Function values at tensorial points.
   ************************************************************************************************/
  template <typename abscissa_float_t, std::size_t sizeT, class input_array_t, class hyEdgeT>
  std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dimension()>
  bulk_values(const std::array<abscissa_float_t, sizeT>& abscissas,
              const input_array_t& lambda_values,
              hyEdgeT& hyper_edge,
              const lSol_float_t time = 0.) const
  {
    // using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
    input_array_t lambda_values_loc = node_dof_to_edge_dof(lambda_values, hyper_edge);
    SmallVec<n_loc_dofs_, lSol_float_t> coefficients =
      solve_local_problem(lambda_values_loc, 1U, hyper_edge, time);
    SmallVec<n_shape_fct_, lSol_float_t> coeffs;
    SmallVec<static_cast<unsigned int>(sizeT), abscissa_float_t> helper(abscissas);

    // for(unsigned int i = 0; i < 2; ++i)
    //   for(unsigned int j = 0; j < 3; ++j)
    //     std::cout << lambda_values_loc[i][j] << " ";
    // std::cout << std::endl;

    // std::cout << coefficients;

    std::array<std::array<lSol_float_t, Hypercube<hyEdge_dimT>::pow(sizeT)>, system_dimension()>
      point_vals, result;

    for (unsigned int dim = 0; dim < space_dim; ++dim)
    {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[(2 * space_dim + dim) * n_shape_fct_ + i];
      for (unsigned int pt = 0; pt < Hypercube<hyEdge_dimT>::pow(sizeT); ++pt)
        point_vals[dim][pt] = integrator::shape_fun_t::template lin_comb_fct_val<float>(
          coeffs, Hypercube<hyEdge_dimT>::template tensorial_pt<Point<hyEdge_dimT>>(pt, helper));
    }

    Point<space_dim, lSol_float_t> normal_vector =
      (Point<space_dim, lSol_float_t>)hyper_edge.geometry.inner_normal(0);
    for (unsigned int dim = 0; dim < result.size(); ++dim)
      for (unsigned int q = 0; q < result[dim].size(); ++q)
        result[dim][q] = point_vals[0][q] * normal_vector[dim];

    normal_vector = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(0);
    for (unsigned int dim = 0; dim < result.size(); ++dim)
      for (unsigned int q = 0; q < result[dim].size(); ++q)
        result[dim][q] += point_vals[1][q] * normal_vector[dim];

    normal_vector = (Point<space_dim, lSol_float_t>)hyper_edge.geometry.outer_normal(1);
    for (unsigned int dim = 0; dim < result.size(); ++dim)
      for (unsigned int q = 0; q < result[dim].size(); ++q)
        result[dim][q] += point_vals[2][q] * normal_vector[dim];

    return result;
  }

  // -----------------------------------------------------------------------------------------------
  // Local-operator protocol of GlobalLoop::Generic
  // -----------------------------------------------------------------------------------------------
  //
  // The blocks of the local problem, A x = B lambda + l, flux = C x + D lambda, with lambda the
  // traces in the frame of the nodes (as the global loop stores them): A is assemble_loc_matrix;
  // B, C and D are the legacy maps applied to unit vectors -- the rotation into the frame of the
  // beam with assemble_rhs_from_lambda (B), extract_fluxes_from_coeffs with the rotation back (C),
  // and -tau lambda in the frame of the beam, rotated back (D) -- so the operator is that of
  // trace_to_flux / residual_flux by construction. The leading 2 space_dim n_shape_fct_ unknowns,
  // the stress resultants n and m, have a diagonal mass block (Legendre basis, coefficients
  // constant on a beam) and are eliminated by the factorization. Dirichlet values enter through
  // the local load; constraining the trace on Dirichlet nodes is the loop's business. The blocks
  // serve the layouts of local_operator.hxx (Condensed, Factored, Streaming); the default layout,
  // BeamStiffness, does not use them (see below).

  struct Args
  {
    lSol_float_t time = 0.;
  };
  static constexpr unsigned int n_loc = 2 * hyEdge_dimT * n_glob_dofs_per_node();
  using Load = SmallVec<n_loc_dofs_, lSol_float_t>;
  using Blocks =
    LocalOperator::Blocks<n_loc_dofs_, n_loc, lSol_float_t, 2 * space_dim * n_shape_fct_>;

  template <unsigned int space_dimT>
  static constexpr bool constrained_node(const unsigned int node_type)
  {
    return is_dirichlet<parametersT<space_dimT, lSol_float_t>>(node_type);
  }
  template <class Edge>
  Blocks local_blocks(Edge& edge, GlobalLoop::Stiffness, const Args& args) const
  {
    static_assert(hyEdge_dimT == 1, "TimoshenkoBeam is implemented for beams (hyEdge_dimT = 1).");
    constexpr unsigned int n_node = n_glob_dofs_per_node();
    using Trace = std::array<std::array<lSol_float_t, n_node>, 2 * hyEdge_dimT>;
    Blocks blocks;
    blocks.A = assemble_loc_matrix(edge, args.time);
    for (unsigned int col = 0; col < n_loc; ++col)
    {
      Trace unit{};
      unit[col / n_node][col % n_node] = 1.;
      const auto local = node_dof_to_edge_dof(unit, edge);
      const Load rhs = assemble_rhs_from_lambda(local, edge);
      for (unsigned int row = 0; row < n_loc_dofs_; ++row)
        blocks.B(row, col) = rhs[row];
      Trace local_flux{}, flux{};
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        for (unsigned int k = 0; k < n_node; ++k)
          local_flux[face][k] = -tau_ * local[face][k];
      edge_dof_to_node_dof(local_flux, flux, edge);
      for (unsigned int row = 0; row < n_loc; ++row)
        blocks.D(row, col) = flux[row / n_node][row % n_node];
    }
    for (unsigned int i = 0; i < n_loc_dofs_; ++i)
    {
      Load unit;
      unit[i] = 1.;
      const auto bdr_values = extract_fluxes_from_coeffs(unit, edge);
      Trace local_flux{}, flux{};
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
        for (unsigned int k = 0; k < n_node; ++k)
          local_flux[face][k] = bdr_values(face, k);
      edge_dof_to_node_dof(local_flux, flux, edge);
      for (unsigned int row = 0; row < n_loc; ++row)
        blocks.C(row, i) = flux[row / n_node][row % n_node];
    }
    return blocks;
  }
  template <class Edge>
  auto local_operator(Edge& edge, GlobalLoop::Stiffness kind, const Args& args) const
  {
    return LayoutT::make(*this, edge, kind, args);
  }
  template <class Edge>
  Load local_load(Edge& edge, const Args& args) const
  {
    return assemble_rhs_from_global_rhs(edge, args.time);
  }
  // -----------------------------------------------------------------------------------------------
  // The local problem by its four independent parts (the layout BeamStiffness)
  // -----------------------------------------------------------------------------------------------
  //
  // In the frame of the beam the local matrix couples a stress resultant only with its own
  // displacement or rotation, except for the cross product: the shear force n_2 with the rotation
  // r_1 and n_1 with r_2. So A falls into the problems of (n_0, u_0), of (m_0, r_0), of
  // (n_2, m_1, u_2, r_1) and of (n_1, m_2, u_1, r_2), each with a diagonal mass block of the
  // stress resultants, which is eliminated: LU factors of size n_shape_fct_ twice and
  // 2 n_shape_fct_ twice, instead of 6 n_shape_fct_ for the whole matrix.

  using BeamTrace = std::array<std::array<lSol_float_t, 2 * space_dim>, 2>;  // in the beam's frame
  struct BeamProblem
  {
    static constexpr unsigned int n = n_shape_fct_;
    static constexpr std::array<unsigned int, 2> single_dim = BeamStiffness::single_dim,
                                                 force_dim = BeamStiffness::force_dim,
                                                 moment_dim = BeamStiffness::moment_dim;

    std::array<LocalOperator::detail::Factors<2 * n, n, lSol_float_t>, 2> single;
    std::array<LocalOperator::detail::Factors<4 * n, 2 * n, lSol_float_t>, 2> coupled;
    std::array<std::array<lSol_float_t, n>, 2> shape_at_node;  // of the shape functions
    std::array<lSol_float_t, 2> normal;                        // local normal at the nodes
    lSol_float_t tau;

    // B lambda, as assemble_rhs_from_lambda.
    Load rhs(const BeamTrace& lambda) const
    {
      Load result;
      for (unsigned int node = 0; node < 2; ++node)
        for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
          for (unsigned int i = 0; i < n; ++i)
          {
            const lSol_float_t value = lambda[node][dim] * shape_at_node[node][i];
            result[dim * n + i] -= normal[node] * value;
            result[(2 * space_dim + dim) * n + i] += tau * value;
          }
      return result;
    }
    // x = A^-1 x.
    void solve(Load& x) const
    {
      for (unsigned int c = 0; c < 2; ++c)
      {
        const std::array<unsigned int, 2> first{single_dim[c] * n,
                                                (2 * space_dim + single_dim[c]) * n};
        solve_part(single[c], first, x);
      }
      for (unsigned int c = 0; c < 2; ++c)
      {
        const std::array<unsigned int, 4> first{force_dim[c] * n, moment_dim[c] * n,
                                                (2 * space_dim + force_dim[c]) * n,
                                                (2 * space_dim + moment_dim[c]) * n};
        solve_part(coupled[c], first, x);
      }
    }
    // C x + D lambda, as extract_fluxes_from_coeffs and -tau lambda.
    BeamTrace flux(const Load& x, const BeamTrace& lambda) const
    {
      BeamTrace result;
      for (unsigned int node = 0; node < 2; ++node)
        for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
        {
          lSol_float_t value = -tau * lambda[node][dim];
          for (unsigned int i = 0; i < n; ++i)
            value += shape_at_node[node][i] *
                     (normal[node] * x[dim * n + i] + tau * x[(2 * space_dim + dim) * n + i]);
          result[node][dim] = value;
        }
      return result;
    }

    // K = C A^-1 of part c, stretching and twisting: the fluxes at the two nodes of a load on its
    // unknowns (n, u). One solve with the transposed factors per flux.
    SmallMat<2, 2 * n, lSol_float_t> load_to_flux_single(const unsigned int c) const
    {
      SmallMat<2 * n, 2, lSol_float_t> transposed_k;
      for (unsigned int node = 0; node < 2; ++node)
        for (unsigned int i = 0; i < n; ++i)
        {
          transposed_k(i, node) = normal[node] * shape_at_node[node][i];
          transposed_k(n + i, node) = tau * shape_at_node[node][i];
        }
      single[c].solve(transposed_k, 'T');
      return transposed(transposed_k);
    }
    // The same for bending, part c: the fluxes of the force at the two nodes, then those of the
    // moment, of a load on its unknowns (n_a, m_b, u_a, r_b).
    SmallMat<4, 4 * n, lSol_float_t> load_to_flux_coupled(const unsigned int c) const
    {
      SmallMat<4 * n, 4, lSol_float_t> transposed_k;
      for (unsigned int field = 0; field < 2; ++field)
        for (unsigned int node = 0; node < 2; ++node)
          for (unsigned int i = 0; i < n; ++i)
          {
            transposed_k(field * n + i, 2 * field + node) = normal[node] * shape_at_node[node][i];
            transposed_k((2 + field) * n + i, 2 * field + node) = tau * shape_at_node[node][i];
          }
      coupled[c].solve(transposed_k, 'T');
      return transposed(transposed_k);
    }

   private:
    // Solve with the unknowns of one part: blocks of n unknowns starting at first[.].
    template <class FactorsT, std::size_t n_blocks>
    static void solve_part(const FactorsT& factors,
                           const std::array<unsigned int, n_blocks>& first,
                           Load& x)
    {
      SmallVec<n_blocks * n, lSol_float_t> part;
      for (unsigned int block = 0; block < n_blocks; ++block)
        for (unsigned int i = 0; i < n; ++i)
          part[block * n + i] = x[first[block] + i];
      factors.solve(part);
      for (unsigned int block = 0; block < n_blocks; ++block)
        for (unsigned int i = 0; i < n; ++i)
          x[first[block] + i] = part[block * n + i];
    }
  };
  /*!***********************************************************************************************
   * \brief   Set up and factorize the four parts of the local problem of a beam (the entries are
   *          those of assemble_loc_matrix).
   ************************************************************************************************/
  template <class Edge>
  BeamProblem beam_problem(Edge& edge) const
  {
    static_assert(hyEdge_dimT == 1 && space_dim == 3,
                  "The parts of the local problem are those of a beam in three dimensions.");
    constexpr unsigned int n = n_shape_fct_;
    const SmallVec<2 * space_dim, lSol_float_t> coeff = beam_coefficients(edge);
    BeamProblem problem;
    problem.tau = tau_;
    for (unsigned int node = 0; node < 2; ++node)
    {
      problem.normal[node] = edge.geometry.local_normal(node)[0];
      for (unsigned int i = 0; i < n; ++i)
        problem.shape_at_node[node][i] =
          integrator::template integrate_bdr_phipsi<decltype(Edge::geometry)>(i, 0, node,
                                                                              edge.geometry);
    }
    // mass, gradient, its partner (normal - gradient) and the penalty, as in assemble_loc_matrix
    SmallSquareMat<n, lSol_float_t> mass, grad, grad_t, face;
    for (unsigned int i = 0; i < n; ++i)
      for (unsigned int j = 0; j < n; ++j)
      {
        mass(i, j) = integrator::integrate_vol_phiphi(i, j, edge.geometry);
        grad(i, j) =
          integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dimT, lSol_float_t>,
                                                         decltype(Edge::geometry)>(
            i, j, edge.geometry)[0];
        grad_t(i, j) = -grad(i, j);
        for (unsigned int node = 0; node < 2; ++node)
        {
          const lSol_float_t value =
            integrator::template integrate_bdr_phiphi<decltype(Edge::geometry)>(i, j, node,
                                                                                edge.geometry);
          face(i, j) += tau_ * value;
          grad_t(i, j) += problem.normal[node] * value;
        }
      }
    for (unsigned int c = 0; c < 2; ++c)
    {
      SmallSquareMat<2 * n, lSol_float_t> part;
      for (unsigned int i = 0; i < n; ++i)
        for (unsigned int j = 0; j < n; ++j)
        {
          part(i, j) = mass(i, j) / coeff[BeamProblem::single_dim[c]];
          part(i, n + j) = -grad(i, j);
          part(n + i, j) = grad_t(i, j);
          part(n + i, n + j) = face(i, j);
        }
      problem.single[c] = LocalOperator::detail::Factors<2 * n, n, lSol_float_t>(part);
    }
    for (unsigned int c = 0; c < 2; ++c)
    {
      // unknowns (n_a, m_b, u_a, r_b); the cross product couples n_a with r_b, with this sign
      const unsigned int a = BeamProblem::force_dim[c], b = BeamProblem::moment_dim[c];
      const lSol_float_t sign = a == 1 ? 1. : -1.;
      SmallSquareMat<4 * n, lSol_float_t> part;
      for (unsigned int i = 0; i < n; ++i)
        for (unsigned int j = 0; j < n; ++j)
        {
          part(i, j) = mass(i, j) / coeff[a];
          part(i, 2 * n + j) = -grad(i, j);
          part(i, 3 * n + j) = sign * mass(i, j);
          part(n + i, n + j) = mass(i, j) / coeff[b];
          part(n + i, 3 * n + j) = -grad(i, j);
          part(2 * n + i, j) = grad_t(i, j);
          part(2 * n + i, 2 * n + j) = face(i, j);
          part(3 * n + i, j) = -sign * mass(i, j);
          part(3 * n + i, n + j) = grad_t(i, j);
          part(3 * n + i, 3 * n + j) = face(i, j);
        }
      problem.coupled[c] = LocalOperator::detail::Factors<4 * n, 2 * n, lSol_float_t>(part);
    }
    return problem;
  }
  /*!***********************************************************************************************
   * \brief   The operator of the layout BeamStiffness: the frame of the beam and its stiffnesses.
   *
   * The stiffnesses are the responses to unit differences of the traces between the two nodes
   * (two solves of the four parts); the fluxes of loads take twelve transposed solves. Debug
   * builds check what the formulas of the operator rest on:
   * no flux for rigid motions, and no coupling of shear and curvature.
   ************************************************************************************************/
  template <class Edge>
  BeamStiffness::Operator<lSol_float_t, n_shape_fct_> beam_operator(Edge& edge) const
  {
    static_assert(poly_deg >= 1, "Rigid rotations need linear displacements.");
    const BeamProblem problem = beam_problem(edge);
    auto response = [&problem](const BeamTrace& lambda)
    {
      Load x = problem.rhs(lambda);
      problem.solve(x);
      return problem.flux(x, lambda);
    };
    // the stiffness of a component: flux = -+ k (difference of the traces)
    auto stiffness = [](const BeamTrace& flux, const unsigned int dim)
    { return (flux[0][dim] - flux[1][dim]) / 2; };

    BeamStiffness::Operator<lSol_float_t, n_shape_fct_> op;
    const Point<space_dim, lSol_float_t> tangent = edge.geometry.inner_normal(0),
                                         normal_1 = edge.geometry.outer_normal(0),
                                         normal_2 = edge.geometry.outer_normal(1);
    for (unsigned int d = 0; d < space_dim; ++d)
    {
      op.frame[d] = tangent[d];
      op.frame[space_dim + d] = normal_1[d];
      op.frame[2 * space_dim + d] = normal_2[d];
    }
    op.half_length = problem.normal[1] * edge.geometry.area() / 2;

    BeamTrace stretch_and_rotate{}, shear{};
    for (const unsigned int dim : {0, 3, 4, 5})
      stretch_and_rotate[0][dim] = -0.5, stretch_and_rotate[1][dim] = 0.5;
    for (const unsigned int dim : {1, 2})
      shear[0][dim] = -0.5, shear[1][dim] = 0.5;
    const BeamTrace flux_rotate = response(stretch_and_rotate), flux_shear = response(shear);
    op.axial = stiffness(flux_rotate, 0);
    op.torsion = stiffness(flux_rotate, 3);
    op.bending = {stiffness(flux_rotate, 4), stiffness(flux_rotate, 5)};
    op.shear = {stiffness(flux_shear, 2), stiffness(flux_shear, 1)};
    for (unsigned int c = 0; c < 2; ++c)
    {
      op.load_to_flux_single[c] = problem.load_to_flux_single(c);
      op.load_to_flux_coupled[c] = problem.load_to_flux_coupled(c);
    }

    if constexpr (hy_debug_mode)
    {
      const lSol_float_t h = op.half_length;
      const lSol_float_t scale =
        (op.axial + op.torsion + op.bending[0] + op.bending[1] + op.shear[0] + op.shear[1]) *
        (1. + std::abs(h));
      const lSol_float_t tol = 1e4 * std::numeric_limits<lSol_float_t>::epsilon() * scale;
      // rigid motions: a translation, and a rotation about each axis of the frame
      BeamTrace rigid;
      rigid[0] = {1., 1. + h, 1. - h, 1., 1., 1.};
      rigid[1] = {1., 1. - h, 1. + h, 1., 1., 1.};
      const BeamTrace flux_rigid = response(rigid);
      for (unsigned int node = 0; node < 2; ++node)
        for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
          hy_assert(std::abs(flux_rigid[node][dim]) <= tol,
                    "A rigid motion of the beam must not cause fluxes, but causes "
                      << (double)flux_rigid[node][dim] << " (node " << node << ", component " << dim
                      << ").");
      for (unsigned int c = 0; c < 2; ++c)
        hy_assert(std::abs(stiffness(flux_rotate, BeamProblem::force_dim[c])) <= tol &&
                    std::abs(stiffness(flux_shear, BeamProblem::moment_dim[c])) <= tol,
                  "Shear and curvature of the beam must not be coupled.");
    }
    return op;
  }
  /*!***********************************************************************************************
   * \brief   Squared L2 error of the displacement on one beam, as errors() computes it.
   ************************************************************************************************/
  template <class Edge, class OperatorT>
  SmallVec<1, lSol_float_t> error_sq(Edge& edge,
                                     const OperatorT& op,
                                     const SmallVec<n_loc, lSol_float_t>& lambda,
                                     const Args& args) const
  {
    using parameters = parametersT<decltype(Edge::geometry)::space_dim(), lSol_float_t>;
    Load coefficients;
    if constexpr (requires { op.local_solution(*this, edge, lambda, args); })
      coefficients = op.local_solution(*this, edge, lambda, args);
    else
      coefficients = LocalOperator::Factored::make(*this, edge, GlobalLoop::Stiffness{}, args)
                       .local_solution(*this, edge, lambda, args);
    SmallVec<1, lSol_float_t> error;
    std::array<lSol_float_t, n_shape_fct_> coeffs;
    const int components[3] = {1, -1, -2};  // as errors(): u along the tangent and the normals
    for (unsigned int c = 0; c < 3; ++c)
    {
      for (unsigned int i = 0; i < coeffs.size(); ++i)
        coeffs[i] = coefficients[i + (2 * space_dim + c) * n_shape_fct_];
      error[0] += integrator::template integrate_vol_diffsquare_discanacomp<
        Point<decltype(Edge::geometry)::space_dim(), lSol_float_t>, decltype(Edge::geometry),
        parameters::analytic_result_u, Point<hyEdge_dimT, lSol_float_t>>(coeffs, components[c],
                                                                         edge.geometry, args.time);
    }
    return error;
  }
};  // end of class LengtheningBernoulliBendingBeam

// -------------------------------------------------------------------------------------------------
// assemble_loc_matrix
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t,
          class LayoutT>
template <typename hyEdgeT>
inline SmallSquareMat<
  TimoshenkoBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t, LayoutT>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t, LayoutT>::
  assemble_loc_matrix(hyEdgeT& hyper_edge, const lSol_float_t time) const
{
  // using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  SmallSquareMat<n_loc_dofs_, lSol_float_t> local_mat;
  lSol_float_t vol_integral, face_integral, helper;
  SmallVec<hyEdge_dimT, lSol_float_t> grad_int_vec, normal_int_vec;
  const SmallVec<2 * space_dim, lSol_float_t> extra_coeffs = beam_coefficients(hyper_edge);

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    for (unsigned int j = 0; j < n_shape_fct_; ++j)
    {
      // Integral_element phi_i phi_j dx in diagonal blocks
      vol_integral = integrator::integrate_vol_phiphi(i, j, hyper_edge.geometry);

      grad_int_vec =
        integrator::template integrate_vol_nablaphiphi<SmallVec<hyEdge_dimT, lSol_float_t>,
                                                       decltype(hyEdgeT::geometry)>(
          i, j, hyper_edge.geometry);

      face_integral = 0.;
      normal_int_vec = 0.;
      for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
      {
        helper = integrator::template integrate_bdr_phiphi<decltype(hyEdgeT::geometry)>(
          i, j, face, hyper_edge.geometry);
        face_integral += helper;
        normal_int_vec += helper * hyper_edge.geometry.local_normal(face);
      }

      for (unsigned int dim = 0; dim < 2 * space_dim; ++dim)
      {
        local_mat(dim * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          vol_integral / extra_coeffs[dim];
        local_mat(dim * n_shape_fct_ + i, (2 * space_dim + dim) * n_shape_fct_ + j) -=
          grad_int_vec[0];
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i, dim * n_shape_fct_ + j) +=
          normal_int_vec[0] - grad_int_vec[0];
        local_mat((2 * space_dim + dim) * n_shape_fct_ + i,
                  (2 * space_dim + dim) * n_shape_fct_ + j) += tau_ * face_integral;
      }

      // Consider the cross product
      local_mat(2 * n_shape_fct_ + i, (3 * space_dim + 1) * n_shape_fct_ + j) -= vol_integral;
      local_mat(1 * n_shape_fct_ + i, (3 * space_dim + 2) * n_shape_fct_ + j) += vol_integral;
      local_mat((3 * space_dim + 2) * n_shape_fct_ + i, 1 * n_shape_fct_ + j) -= vol_integral;
      local_mat((3 * space_dim + 1) * n_shape_fct_ + i, 2 * n_shape_fct_ + j) += vol_integral;
    }
  }

  return local_mat;
}  // end of Diffusion::assemble_loc_matrix

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_lambda
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t,
          class LayoutT>
template <typename hyEdgeT, typename SmallMatT>
inline SmallVec<
  TimoshenkoBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t, LayoutT>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t, LayoutT>::
  assemble_rhs_from_lambda(const SmallMatT& lambda_values, hyEdgeT& hyper_edge) const
{
  static_assert(std::is_same<typename SmallMatT::value_type::value_type, lSol_float_t>::value,
                "Lambda values should have same floating point arithmetics as local solver!");
  hy_assert(lambda_values.size() == 2 * hyEdge_dimT,
            "The size of the lambda values should be twice the dimension of a hyperedge.");
  for (unsigned int i = 0; i < 2 * hyEdge_dimT; ++i)
    hy_assert(lambda_values[i].size() == 2 * space_dim * n_shape_bdr_,
              "The size of lambda should be the amount of ansatz functions at boundary.");

  // for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
  //   for (unsigned int j = 0; j < 1 * space_dim; ++j)
  //     std::cout << lambda_values[face][j] << " ";
  // std::cout << std::endl;

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

  // std::cout << right_hand_side << std::endl;
  return right_hand_side;
}  // end of Diffusion::assemble_rhs_from_lambda

// -------------------------------------------------------------------------------------------------
// assemble_rhs_from_global_rhs
// -------------------------------------------------------------------------------------------------

template <unsigned int hyEdge_dimT,
          unsigned int space_dim,
          unsigned int poly_deg,
          unsigned int quad_deg,
          template <unsigned int, typename>
          typename parametersT,
          typename lSol_float_t,
          class LayoutT>
template <typename hyEdgeT>
inline SmallVec<
  TimoshenkoBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t, LayoutT>::
    n_loc_dofs_,
  lSol_float_t>
TimoshenkoBeam<hyEdge_dimT, space_dim, poly_deg, quad_deg, parametersT, lSol_float_t, LayoutT>::
  assemble_rhs_from_global_rhs(hyEdgeT& hyper_edge, const unsigned int dim) const
{
  using parameters = parametersT<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>;
  // constexpr unsigned int n_dofs_lap = n_loc_dofs_ / 2;
  SmallVec<n_loc_dofs_, lSol_float_t> right_hand_side;
  lSol_float_t integral;

  for (unsigned int i = 0; i < n_shape_fct_; ++i)
  {
    right_hand_side[2 * space_dim * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_n, Point<hyEdge_dimT, lSol_float_t>>(i, 1, hyper_edge.geometry,
                                                                         0.);
    right_hand_side[(2 * space_dim + 1) * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_n, Point<hyEdge_dimT, lSol_float_t>>(i, -1, hyper_edge.geometry,
                                                                         0.);
    right_hand_side[(2 * space_dim + 2) * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_n, Point<hyEdge_dimT, lSol_float_t>>(i, -2, hyper_edge.geometry,
                                                                         0.);

    right_hand_side[3 * space_dim * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_m, Point<hyEdge_dimT, lSol_float_t>>(i, 1, hyper_edge.geometry,
                                                                         0.);
    right_hand_side[(3 * space_dim + 1) * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_m, Point<hyEdge_dimT, lSol_float_t>>(i, -1, hyper_edge.geometry,
                                                                         0.);
    right_hand_side[(3 * space_dim + 2) * n_shape_fct_ + i] =
      integrator::template integrate_vol_phivecfunccomp<
        Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>, decltype(hyEdgeT::geometry),
        parameters::right_hand_side_m, Point<hyEdge_dimT, lSol_float_t>>(i, -2, hyper_edge.geometry,
                                                                         0.);
    for (unsigned int face = 0; face < 2 * hyEdge_dimT; ++face)
    {
      if (is_dirichlet<parameters>(hyper_edge.node_descriptor[face]))
      {
        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, 1, hyper_edge.geometry, 0.);
        right_hand_side[(0 * space_dim + 0) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(2 * space_dim + 0) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, -1, hyper_edge.geometry, 0.);
        right_hand_side[(0 * space_dim + 1) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(2 * space_dim + 1) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_u,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, -2, hyper_edge.geometry, 0.);
        right_hand_side[(0 * space_dim + 2) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(2 * space_dim + 2) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, 1, hyper_edge.geometry, 0.);
        right_hand_side[(1 * space_dim + 0) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(3 * space_dim + 0) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, -1, hyper_edge.geometry, 0.);
        right_hand_side[(1 * space_dim + 1) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(3 * space_dim + 1) * n_shape_fct_ + i] += tau_ * integral;

        integral = integrator::template integrate_bdr_phivecfunccomp<
          Point<decltype(hyEdgeT::geometry)::space_dim(), lSol_float_t>,
          decltype(hyEdgeT::geometry), parameters::dirichlet_value_phi,
          Point<hyEdge_dimT, lSol_float_t>>(i, face, -2, hyper_edge.geometry, 0.);
        right_hand_side[(1 * space_dim + 2) * n_shape_fct_ + i] -=
          hyper_edge.geometry.local_normal(face).operator[](0) * integral;
        right_hand_side[(3 * space_dim + 2) * n_shape_fct_ + i] += tau_ * integral;
      }
    }
  }

  return right_hand_side;
}  // end of Bilaplacian::assemble_rhs_from_global_rhs

}  // namespace LocalSolver
