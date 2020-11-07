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
 * \brief   Combine local solver and global information for parabolic problems.
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
class Parabolic
{
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
  Parabolic(const typename TopologyT::constructor_value_type& construct_topo,
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
  Parabolic(const typename TopologyT::constructor_value_type& construct_topo,
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
  Parabolic(const typename TopologyT::constructor_value_type& construct_topo)
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
   * \brief   Returns vector of appropriate size for the predefined problem.
   *
   * Returns a vector containing only the value zero, but of the size \f$n\f$ which is also the
   * number which is returned if \c size_of_system() is evaluated.
   *
   * \retval  zero          A vector of the correct size for the unknowns of the given problem.
   ************************************************************************************************/
  LargeVecT return_zero_vector() const { return LargeVecT(hyper_graph_.n_global_dofs(), 0.); }
  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product.
   *
   * Function that evaluates the condensed, matrix-free version of the matrix-vector product
   * \f$A x = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method that needs to be
   * inverted for a time step, \f$x\f$ is the vector of parameters to define the skeletal variable
   * \f$\lambda\f$, and \f$y\f$ is the resulting vector, which has the same size as the input vector
   * \f$x\f$.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time at which the new time step will end.
   * \retval  y_vec         A vector containing the product \f$y = Ax\f$.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT matrix_vector_multiply(const LargeVecT& x_vec, const dof_value_t time = 0.)
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
      if constexpr (not_uses_geometry<LocalSolverT,
                                      std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                 2 * TopologyT::hyEdge_dim()>(
                                        std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                   2 * TopologyT::hyEdge_dim()>&)>::value)
        hyEdge_dofs = local_solver_.numerical_flux_from_lambda(hyEdge_dofs, time);
      else
        hyEdge_dofs = local_solver_.numerical_flux_from_lambda(hyEdge_dofs, hyper_edge, time);

      // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().add_to_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                        hyEdge_dofs[hyNode]);
    });

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product.
   *
   * Function that evaluates the condensed, matrix-free version of the matrix-vector product
   * \f$A x = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method that needs to be
   * inverted to do a time step, \f$x\f$ is the vector of parameters to define the skeletal variable
   * \f$\lambda\f$, and \f$y\f$ is the resulting vector, which has the same size as the input vector
   * \f$x\f$.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time at which the time step ends.
   * \retval  y_vec         A vector containing the product \f$y = Ax\f$.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  std::vector<dof_value_t> total_flux_vector(const std::vector<dof_value_t>& x_vec,
                                             const dof_value_t time = 0.)
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
      if constexpr (not_uses_geometry<LocalSolverT,
                                      std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                 2 * TopologyT::hyEdge_dim()>(
                                        std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                   2 * TopologyT::hyEdge_dim()>&)>::value)
      {
        hyEdge_dofs = local_solver_.numerical_flux_total(hyEdge_dofs, time);
      }
      else
      {
        hyEdge_dofs = local_solver_.numerical_flux_total(hyEdge_dofs, hyper_edge, time);
      }
      // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().add_to_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                        hyEdge_dofs[hyNode]);
    });

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Set data using the result of the old time step.
   *
   * \param   x_vec         A \c std::vector containing the input vector \f$x\f$.
   * \param   time          Time at which the old time step ended.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  void set_data(const LargeVecT& x_vec, const dof_value_t time = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

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
      if constexpr (not_uses_geometry<LocalSolverT,
                                      std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                 2 * TopologyT::hyEdge_dim()>(
                                        std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                   2 * TopologyT::hyEdge_dim()>&)>::value)
      {
        local_solver_.set_data(hyEdge_dofs, time);
      }
      else
      {
        local_solver_.set_data(hyEdge_dofs, hyper_edge, time);
      }
    });
  }
  /*!***********************************************************************************************
   * \brief   Evaluate the initial flux of the problem.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time for initial data.
   * \retval  y_vec         A vector containing the initial fluxes.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT initial_flux_vector(const LargeVecT& x_vec, const dof_index_t time = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    std::vector<dof_value_t> vec_Ax(x_vec.size(), 0.);
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
      if constexpr (not_uses_geometry<LocalSolverT,
                                      std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                 2 * TopologyT::hyEdge_dim()>(
                                        std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                   2 * TopologyT::hyEdge_dim()>&)>::value)
      {
        hyEdge_dofs = local_solver_.numerical_flux_initial(hyEdge_dofs, time);
      }
      else
      {
        hyEdge_dofs = local_solver_.numerical_flux_initial(hyEdge_dofs, hyper_edge, time);
      }
      // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().set_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                     hyEdge_dofs[hyNode]);
    });

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Calculate L2 error.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time at which error is evaluated.
   * \retval  error         L2 error.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  dof_value_t calculate_L2_error(const LargeVecT& x_vec, const dof_value_t time = 0.)
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
    constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

    dof_value_t result = 0.;

    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

    // Calculate errors by iteration over all hyperedges.
    std::for_each(hyper_graph_.begin(), hyper_graph_.end(), [&](auto hyper_edge) {
      // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
      hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec,
                                                     hyEdge_dofs[hyNode]);

      // Turn degrees of freedom of x_vec that have been stored locally into local errors.
      // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
      if constexpr (not_uses_geometry<LocalSolverT,
                                      std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                 2 * TopologyT::hyEdge_dim()>(
                                        std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                   2 * TopologyT::hyEdge_dim()>&)>::value)
      {
        result += local_solver_.calc_L2_error_squared(hyEdge_dofs, time);
      }
      else
      {
        result += local_solver_.calc_L2_error_squared(hyEdge_dofs, hyper_edge, time);
      }
    });

    hy_assert(result >= 0., "The squared error must be non-negative, but was " << result);
    return std::sqrt(result);
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
   * \retval  n             Size of condensed system of equations.
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
   * \param   time          Time at which analytic functions are evaluated.
   * \retval  file          A file in the output directory.
   ************************************************************************************************/
  void plot_solution(const std::vector<dof_value_t>& lambda, const dof_value_t time = 0.)
  {
    plot(hyper_graph_, local_solver_, lambda, plot_options, time);
  }
};  // end of class Parabolic

}  // end of namespace GlobalLoop
