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
 * \brief   Combine local solver and global information for eigenvalue problems with shifted inverse
 *          approach.
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
class ShiftedEigenvalue
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
   * \brief   Vector containing the indices of Dirichlet type nodes.
   ************************************************************************************************/
  std::vector<dof_index_t> dirichlet_indices_;
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
  ShiftedEigenvalue(const typename TopologyT::constructor_value_type& construct_topo,
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
  ShiftedEigenvalue(const typename TopologyT::constructor_value_type& construct_topo,
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
  ShiftedEigenvalue(const typename TopologyT::constructor_value_type& construct_topo)
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
   * \brief   Return indices of Dirichlet degrees of freedom.
   *
   * \tparam  hyNode_index_t      Typename of the index type for hypernodes.
   * \retval  dirichlet_indices_  Vector containing all dirichlet dof indices.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  std::vector<unsigned int> dirichlet_indices()
  {
    constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();

    std::array<dof_index_t, LocalSolverT::n_glob_dofs_per_node()> dof_indices;
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<unsigned int, 2 * hyEdge_dim> hyNode_types;

    std::for_each(hyper_graph_.begin(), hyper_graph_.end(), [&](auto hyper_edge) {
      hyNode_types = local_solver_.node_types(hyper_edge);
      hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();

      for (unsigned int i = 0; i < hyNode_types.size(); ++i)
        if (hyNode_types[i] == 1)
        {
          hyper_graph_.hyNode_factory().get_dof_indices(hyEdge_hyNodes[i], dof_indices);
          for (unsigned int j = 0; j < dof_indices.size(); ++j)
            dirichlet_indices_.push_back(dof_indices[j]);
        }
    });

    std::sort(dirichlet_indices_.begin(), dirichlet_indices_.end());
    auto last = std::unique(dirichlet_indices_.begin(), dirichlet_indices_.end());
    dirichlet_indices_.erase(last, dirichlet_indices_.end());

    return dirichlet_indices_;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product.
   *
   * Evaluate \f$ (A - \sigma M) x = y\f$.
   *
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   sigma         Approximation to eigenvalue / shifting parameter.
   * \retval  y_vec         A vector containing the product \f$y = Ax\f$.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT matrix_vector_multiply(const LargeVecT& x_vec, const dof_value_t sigma = 0.)
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
      if constexpr (not_uses_geometry<LocalSolverT,
                                      std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                 2 * TopologyT::hyEdge_dim()>(
                                        std::array<std::array<dof_value_t, n_dofs_per_node>,
                                                   2 * TopologyT::hyEdge_dim()>&)>::value)
        local_solver_.numerical_flux_from_lambda(hyEdge_dofs_old, hyEdge_dofs_new, sigma);
      else
        local_solver_.numerical_flux_from_lambda(hyEdge_dofs_old, hyEdge_dofs_new, hyper_edge,
                                                 sigma);

      // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode)
        hyper_graph_.hyNode_factory().add_to_dof_values(hyEdge_hyNodes[hyNode], vec_Ax,
                                                        hyEdge_dofs_new[hyNode]);
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
   * \retval  n             Size of condensed global system of equations.
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
};  // end of class ShiftedEigenvalue

}  // end of namespace GlobalLoop
