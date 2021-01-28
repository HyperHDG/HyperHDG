#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/plot.hxx>
#include <algorithm>
#include <array>
#include <cmath>

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Combine local solver and global information for nonlinear eigenvalue problems.
 *
 * \tparam  TopologyT       Class type containing topological information.
 * \tparam  GeometryT       Class type containing geometrical information.
 * \tparam  NodeDescriptorT Class type containing the information of nodes of hyperedges.
 * \tparam  LocalSolverT    Class type of the local solver.
 * \tparam  LargeVecT       Clas type of large, global vector.
 * \tparam  dof_index_t     Index type of hyperedges. Default is \c unsigned \c int.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template <class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename LargeVecT = std::vector<double>,
          typename dof_index_t = unsigned int>
class NonlinearEigenvalue
{
  /*!***********************************************************************************************
   * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
   ************************************************************************************************/
  HAS_MEMBER_FUNCTION(trace_to_flux, has_trace_to_flux);
  /*!***********************************************************************************************
   * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
   ************************************************************************************************/
  HAS_MEMBER_FUNCTION(jacobian_of_trace_to_flux, has_jacobian_of_trace_to_flux);
  /*!***********************************************************************************************
   * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
   ************************************************************************************************/
  HAS_MEMBER_FUNCTION(make_initial, has_make_initial);

  /*!***********************************************************************************************
   * \brief   Floating type is determined by floating type of large vector's entries.
   ************************************************************************************************/
  using dof_value_t = typename LargeVecT::value_type;

 private:
  /*!***********************************************************************************************
   * \brief   Instantiation of a hypergraph.
   ************************************************************************************************/
  HDGHyperGraph<LocalSolverT::n_glob_dofs_per_node(),
                TopologyT,
                GeometryT,
                NodeDescriptorT,
                typename LocalSolverT::data_type>
    hyper_graph_;
  /*!***********************************************************************************************
   * \brief   Instantiation of a local solver.
   ************************************************************************************************/
  const LocalSolverT local_solver_;
  /*!***********************************************************************************************
   * \brief   Struct encoding the options for plotting.
   ************************************************************************************************/
  PlotOptions plot_options;

 public:
  /*!***********************************************************************************************
   * \brief   Abstract problem constructor.
   *
   * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
   * hyperedge.
   *
   * \param   construct_topo    Information to construct a topology.
   * \param   construct_geom    Information to construct a geometry.
   * \param   construct_loc_sol Information to construct a local solver.
   ************************************************************************************************/
  NonlinearEigenvalue(const typename TopologyT::constructor_value_type& construct_topo,
                      const typename GeometryT::constructor_value_type& construct_geom,
                      const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo, construct_geom), local_solver_(construct_loc_sol)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Abstract problem constructor.
   *
   * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
   * hyperedge.
   *
   * \param   construct_topo    Information to construct a topology.
   * \param   construct_loc_sol Information to construct a local solver.
   ************************************************************************************************/
  NonlinearEigenvalue(const typename TopologyT::constructor_value_type& construct_topo,
                      const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo), local_solver_(construct_loc_sol)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Abstract problem constructor.
   *
   * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
   * hyperedge.
   *
   * \param   construct_topo    Information to construct a topology.
   ************************************************************************************************/
  NonlinearEigenvalue(const typename TopologyT::constructor_value_type& construct_topo)
  : hyper_graph_(construct_topo)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product.
   *
   * This function corresponds to the evaluation of the residual. Here, the vector contains the
   * eigenfunction representation, while the floating point is the eigenvalue.
   *
   * \param   x_vec         A vector containing the input vector.
   * \param   eig           Eigenvalue.
   * \retval  y_vec         A vector containing the residual.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT trace_to_flux(const LargeVecT& x_vec, const dof_value_t eig = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    LargeVecT vec_Ax(x_vec.size(), 0.);
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs_old,
      hyEdge_dofs_new;

    // Do matrix--vector multiplication by iterating over all hyperedges.
    std::for_each(hyper_graph_.begin(), hyper_graph_.end(), [&](auto hyper_edge) {
      // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
      hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
      {
        hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec,
                                                     hyEdge_dofs_old[hyNode]);
        hyEdge_dofs_new[hyNode].fill(0.);
      }

      // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
      if constexpr (
        has_trace_to_flux<
          LocalSolverT,
          std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&(
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            dof_value_t)>::value)
        local_solver_.trace_to_flux(hyEdge_dofs_old, hyEdge_dofs_new, eig);
      else if constexpr (
        has_trace_to_flux<
          LocalSolverT,
          std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&(
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            decltype(hyper_edge)&, dof_value_t)>::value)
        local_solver_.trace_to_flux(hyEdge_dofs_old, hyEdge_dofs_new, hyper_edge, eig);
      else
        hy_assert(false, "Function seems not to be implemented!");

      // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().add_to_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                        hyEdge_dofs_new[hyNode]);
    });

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product.
   *
   * This function corresponds to evaluating the Jacobian at point \c x_val, \c eig_val in direction
   * \c x_vec, \c eig.
   *
   * \param   x_vec         Direction in which Jacobian is evaluated.
   * \param   eig           Direction in which Jacobian is evaluated.
   * \param   x_val         Point at which Jacobian is evaluated.
   * \param   eig_val       Point at which Jacobian is evaluated.
   * \retval  y_vec         Corresponds to directional derivative.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT jacobian_of_trace_to_flux(const LargeVecT& x_vec,
                                      const dof_value_t eig,
                                      const LargeVecT& x_val,
                                      const dof_value_t eig_val)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    LargeVecT vec_Ax(x_vec.size(), 0.);
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs_old,
      hyEdge_dofs_new, hyEdge_vals;

    // Do matrix--vector multiplication by iterating over all hyperedges.
    std::for_each(hyper_graph_.begin(), hyper_graph_.end(), [&](auto hyper_edge) {
      // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
      hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
      {
        hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec,
                                                     hyEdge_dofs_old[hyNode]);
        hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_val,
                                                     hyEdge_vals[hyNode]);
        hyEdge_dofs_new[hyNode].fill(0.);
      }

      // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
      if constexpr (
        has_jacobian_of_trace_to_flux<
          LocalSolverT,
          std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&(
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            dof_value_t,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            dof_value_t)>::value)
      {
        local_solver_.jacobian_of_trace_to_flux(hyEdge_dofs_old, hyEdge_dofs_new, eig, hyEdge_vals,
                                                eig_val);
      }
      else if constexpr (
        has_jacobian_of_trace_to_flux<
          LocalSolverT,
          std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&(
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            dof_value_t,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>&,
            dof_value_t, decltype(hyper_edge)&)>::value)
      {
        local_solver_.jacobian_of_trace_to_flux(hyEdge_dofs_old, hyEdge_dofs_new, eig, hyEdge_vals,
                                                eig_val, hyper_edge);
      }
      else
        hy_assert(false, "Function seems not to be implemented!");

      // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().add_to_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                        hyEdge_dofs_new[hyNode]);
    });

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Create initial starting value for Newton's methods.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$ which should be zero.
   * \param   eig           Eigenvalue (approximation).
   * \retval  y_vec         A vector containing initial flux vector.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT make_initial(const LargeVecT& x_vec, const dof_value_t eig = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    LargeVecT vec_Ax(x_vec.size(), 0.);
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

    // Do matrix--vector multiplication by iterating over all hyperedges.
    std::for_each(hyper_graph_.begin(), hyper_graph_.end(), [&](auto hyper_edge) {
      // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
      hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec,
                                                     hyEdge_dofs[hyNode]);

      // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
      if constexpr (has_make_initial<LocalSolverT,
                                     std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                2 * TopologyT::hyEdge_dim()>&(
                                       std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                  2 * TopologyT::hyEdge_dim()>&,
                                       dof_value_t)>::value)
      {
        local_solver_.make_initial(hyEdge_dofs, eig);
      }
      else if constexpr (has_make_initial<LocalSolverT,
                                          std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                     2 * TopologyT::hyEdge_dim()>&(
                                            std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                       2 * TopologyT::hyEdge_dim()>&,
                                            decltype(hyper_edge)&, dof_value_t)>::value)
      {
        local_solver_.make_initial(hyEdge_dofs, hyper_edge, eig);
      }
      else
        hy_assert(false, "Function seems not to be implemented!");

      // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().set_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                     hyEdge_dofs[hyNode]);
    });

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Determine size of condensed system for the skeletal unknowns.
   *
   * Function that returns the size \f$n\f$ of the \f$n \times n\f$ linear, sparse system
   * \f$Ax = b\f$ that is solved by the program in a matrix-free fashion.
   *
   * This function is needed to define a \c LinearOperator from Python's \c scipy.sparse.linalg
   * package which can be used to define iterative solvers for sparse systems.
   *
   * \retval  n             An \c int which Python needs and actually is a parsed \c unsigned
   *                        \c int.
   ************************************************************************************************/
  dof_index_t size_of_system() const { return hyper_graph_.n_global_dofs(); }
  /*!***********************************************************************************************
   * \brief   Set plot option and return old plot option.
   *
   * Function to set and / or read the current plot option.
   *
   * \param   option        A \c std::string containing the plot option to be considered.
   * \param   value         A \c std::string containing the new value of the considered option.
   *                        If empty, the old value is kept.
   * \retval  opt_value     A \c std::string containing the value of the plot option.
   ************************************************************************************************/
  std::string plot_option(const std::string& option, std::string value = "")
  {
    return set_plot_option(plot_options, option, value);
  }
  /*!***********************************************************************************************
   * \brief   Plot solution in vtu format.
   *
   * Function that plots the solution of the problem to a predefined file.
   *
   * \param   lambda        A vector of unknowns containing the data vector.
   * \param   time          Time.
   * \retval  file          A file in the output directory.
   ************************************************************************************************/
  void plot_solution(const std::vector<dof_value_t>& lambda, const dof_value_t time = 0.)
  {
    plot(hyper_graph_, local_solver_, lambda, plot_options, time);
  }
  /*!***********************************************************************************************
   * \brief   Return refinement level.
   ************************************************************************************************/
  unsigned int get_refinement() { return hyper_graph_.get_refinement(); }
  /*!***********************************************************************************************
   * \brief   Set refinement level.
   ************************************************************************************************/
  void set_refinement(unsigned int level) { hyper_graph_.set_refinement(level); }
};  // end of class NonlinearEigenvalue

}  // end of namespace GlobalLoop
