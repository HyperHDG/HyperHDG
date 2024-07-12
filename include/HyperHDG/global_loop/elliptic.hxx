#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/global_loop/prototype.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/plot.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Combine local solver and global information for elliptic problems.
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
class Elliptic
{
 private:
  /*!***********************************************************************************************
   * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
   ************************************************************************************************/
  HAS_MEMBER_FUNCTION(trace_to_flux, has_trace_to_flux);
  /*!***********************************************************************************************
   * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
   ************************************************************************************************/
  HAS_MEMBER_FUNCTION(residual_flux, has_residual_flux);
  /*!***********************************************************************************************
   * \brief   Prepare struct to check for function to exist (cf. compile_time_tricks.hxx).
   ************************************************************************************************/
  HAS_MEMBER_FUNCTION(errors, has_errors);
  /*!***********************************************************************************************
   * \brief   Some constant variable that might be helpful.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
  /*!***********************************************************************************************
   * \brief   Some constant variable that might be helpful.
   ************************************************************************************************/
  static constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();

  /*!***********************************************************************************************
   * \brief   Floating type is determined by floating type of large vector's entries.
   ************************************************************************************************/
  using dof_value_t = typename LargeVecT::value_type;
  /*!***********************************************************************************************
   * \brief   Instantiation of a hypergraph.
   ************************************************************************************************/
  HDGHyperGraph<LocalSolverT::n_glob_dofs_per_node(),
                TopologyT,
                GeometryT,
                NodeDescriptorT,
                typename LocalSolverT::data_type>
    hyper_graph_, hyper_graph_ref_;
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
  Elliptic(const typename TopologyT::constructor_value_type& construct_topo,
           const typename GeometryT::constructor_value_type& construct_geom,
           const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo, construct_geom), hyper_graph_ref_(construct_topo, construct_geom),
            local_solver_(construct_loc_sol)
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
  Elliptic(const typename TopologyT::constructor_value_type& construct_topo,
           const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo), hyper_graph_ref_(construct_topo), local_solver_(construct_loc_sol)
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
  Elliptic(const typename TopologyT::constructor_value_type& construct_topo)
  : hyper_graph_(construct_topo), hyper_graph_ref_(construct_topo)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Read indices of Dirichlet type hypernodes/faces.
   *
   * Read the indices ot the hypernodes/faces that are of Dirichlet type and therefore do not
   * contain degrees of freedom that are allowed to change within iterations of the iterative
   * solver and other processes. In contrast, these degrees of freedom are set by the user.
   *
   * The user creates a vector that contains the coefficients of the corresponding degrees of
   * freedom (read by this function) and defines the Dirichlet values by this choice. The
   * remaining elements of the global vector of unknowns (which is \b not the vector \c indices
   * are supposed to be zero).
   *
   * \param   indices         A vector containing the (global) indices of Dirichlet type hypernodes.
   ************************************************************************************************/
  void read_dirichlet_indices(const std::vector<unsigned int>& indices)
  {
    dirichlet_indices_.resize(indices.size());
    for (unsigned int i = 0; i < indices.size(); ++i)
    {
      hy_assert(
        (dof_index_t)indices[i] >= 0 && (dof_index_t)indices[i] < hyper_graph_.n_global_dofs(),
        "All indices of Dirichlet nodes need to be larger than or equal to zero and "
          << "smaller than the total amount of degrees of freedom." << std::endl
          << "In this case, the index is " << indices[i] << " and the total amount of "
          << "hypernodes is " << hyper_graph_.n_global_dofs() << ".");
      dirichlet_indices_[i] = (dof_index_t)indices[i];
    }
    std::sort(dirichlet_indices_.begin(), dirichlet_indices_.end());
    auto last = std::unique(dirichlet_indices_.begin(), dirichlet_indices_.end());
    dirichlet_indices_.erase(last, dirichlet_indices_.end());
  }
  /*!***********************************************************************************************
   * \brief   Returns vector of appropriate size for the predefined problem.
   *
   * Returns a vector containing only the value zero, but of the size \f$n\f$ which is also the
   * number which is returned if \c size_of_system() is evaluated.
   *
   * \retval  zero            A vector of the correct size for the unknowns of the given problem.
   ************************************************************************************************/
  LargeVecT zero_vector() const { return LargeVecT(hyper_graph_.n_global_dofs(), 0.); }
  /*!***********************************************************************************************
   * \brief   VP: testing out the function to find the interior coefficients for u
   *
   * \param   coeffs_in     A vector containing the skeletal lambda coefficients.
   * \param   time          Time at which the new time step will end.
   * \retval  all_coeffs    A vector containing the interior u coefficients.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT refine_solution(const LargeVecT& coeffs_in, const dof_value_t time = 0.)
  {
    /* The function currently maps the previous values to the corresponding locations.
    - That is: interior to new nodes, old nodes to refined nodes.
    NOTE:
    - The implementation should work for cases:
        'coarse.refinement * 2 = refined.refinement' splits each edge into 2 ** hyEdge_dim edges.
    AND
        'coarse.refinement == 1 -> refined.refinement = K'
    *
    - Refining from coarse.refinement != 1 to refined.refinement != 2 * coarse.refinement, hasn't been tested.
    But, the following has been built upon an assumption that each refined HyEdge belongs to only one coarse HyEdge.
    */
    const unsigned int target_refinement = 2;

    unsigned int n_splits = 2;
    // Refine the auxilary Hyper Graph
    if (hyper_graph_.get_refinement() == 1) {
      n_splits = target_refinement;
      hyper_graph_ref_.set_refinement(target_refinement);
    } else {
      hyper_graph_ref_.set_refinement(hyper_graph_.get_refinement() * 2);
    }

    // Check the edge counts (just in case)
    // printf("\nn_edges coarse: %d\n", hyper_graph_.n_hyEdges());
    // printf("n_edges refined: %d\n", hyper_graph_ref_.n_hyEdges());

    // Initialize the output vector for the new dof values:
    LargeVecT coeffs_out(hyper_graph_ref_.n_global_dofs(), 0.);

    // Auxilary vectors / matrices for the evaluations:
    std::array<dof_value_t, n_dofs_per_node> hyNode_dofs_coarse;
    SmallVec<2 * hyEdge_dim, hyNode_index_t> hyEdge_hyNodes;
    std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs_old;

    // Number of nodes on the refined graph for each coarse node (hyNode_dim???):
    const unsigned int n_subNodes = pow(n_splits, hyEdge_dim - 1);

    // Constant values of the coarse hyper graph:
    const unsigned int n_hyNodes_coarse = hyper_graph_.n_hyNodes();
    const unsigned int n_hyEdges_coarse = hyper_graph_.n_hyEdges();

    // Number of coefficients on the node: n_dofs_per_node

    // Indexing 'unsigned int' variables: (for clarity)
    unsigned int hyNode_ref = 0;

    // Go through the Hyper Nodes of the coarse graph:
    for (unsigned int hyNode_coarse = 0; hyNode_coarse < n_hyNodes_coarse; ++hyNode_coarse) {
      // Take the dof-values from the input vector corresponding to the hyper node:
      hyper_graph_.hyNode_factory().get_dof_values(hyNode_coarse, coeffs_in, hyNode_dofs_coarse);

      // Add the dof-values to the output:
      // (For each subdivision of a coarse node one after another:)
      for (unsigned int sub_node = 0; sub_node < n_subNodes; ++sub_node) {
        // Evaluate the index of the new node:
        /* hyNode_ref = hyNode_coarse * n_subNodes + sub_node; */

        // Perform the projection to get the nodal coefficients:
        // Mapping : node of coarse edge -> node of refined edge
        auto hyNode_dofs_ref = local_solver_.node_to_node(hyNode_dofs_coarse);

        // Add the nodal dofs to the global dofs (old values atm)
        hyper_graph_ref_.hyNode_factory().set_dof_values(hyNode_ref, coeffs_out, hyNode_dofs_coarse);

        ++hyNode_ref;
      }
    }

    // Go over the hyper edges of the coarse graph:
    for (unsigned int hyEdge_coarse = 0; hyEdge_coarse < n_hyEdges_coarse; ++hyEdge_coarse) {
      // Get the corresponding hyper edge:
      auto hyper_edge = hyper_graph_[hyEdge_coarse];

      // Get the node indices corresponding to the considered hyper edge:
      hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();

      // Go over the hyper nodes of the edge to get the dofs:
      for (unsigned int hyNode = 0; hyNode < hyEdge_hyNodes.size(); ++hyNode) {
        hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], coeffs_in, hyEdge_dofs_old[hyNode]);
      }

      // Evaluate the local interior coefficients:
      auto interior_dofs = local_solver_.local_interior(hyEdge_dofs_old, hyper_edge);

      // Map the interior coefficients to the new nodes:
      for (unsigned int dim = 0; dim < hyEdge_dim; ++dim) {
        for (unsigned int split = 0; split < n_splits - 1; ++split) {
          // Go over the nodes which need to be added:
          for (unsigned int sub_node = 0; sub_node < n_subNodes; ++sub_node) {
            // Evaluate the node index:
            /* hyNode_ref = n_hyNodes_coarse * n_subNodes
                    + hyEdge_coarse * (hyEdge_dim * n_subNodes * (n_splits - 1))
                    + dim * n_subNodes * (n_splits - 1)
                    + split * (n_splits - 1)
                    + sub_node; */

            // Perform the projection to get the nodal coefficients:
            // Mapping : interior of coarse edge -> node of refined edge
            // auto hyNode_dofs_ref = local_solver_.interior_to_node(interior_dofs, hyper_edge_ref, 0);

            // Add the nodal dofs to the global dofs: (old dofs atm)
            hyper_graph_ref_.hyNode_factory().set_dof_values(hyNode_ref, coeffs_out, interior_dofs);

            // Increment the refined hyper node
            ++hyNode_ref;
          }
        }
      }
    }

    return coeffs_out;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product.
   *
   * Function that evaluates the condensed, matrix-free version of the matrix-vector product
   * \f$A x = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method, \f$x\f$ is the
   * vector of parameters to define the skeletal variable \f$\lambda\f$, and \f$y\f$ is the
   * resulting vector, which has the same size as the input vector \f$x\f$.
   *
   * \tparam  hyNode_index_t  Typename of the hypernode index. Defaults to \c unsigned \c int.
   * \param   x_vec           A vector containing the input vector \f$x\f$.
   * \param   time            Time at which given analyitic functions will be evaluated.
   * \retval  y_vec           A vector containing the product \f$y = Ax\f$.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT trace_to_flux(const LargeVecT& x_vec, const dof_value_t time = 0.)
  {
    auto vec_Ax = prototype_mat_vec_multiply(trace_to_flux, has_trace_to_flux);

    // Set all Dirichlet values to zero.
    for (dof_index_t i = 0; i < dirichlet_indices_.size(); ++i)
    {
      hy_assert(dirichlet_indices_[i] >= 0 && dirichlet_indices_[i] < hyper_graph_.n_global_dofs(),
                "All indices of Dirichlet nodes need to be larger than or equal to zero and "
                  << "smaller than the total amount of degrees of freedom." << std::endl
                  << "In this case, the index is " << dirichlet_indices_[i] << " and the total "
                  << "amount of hypernodes is " << hyper_graph_.n_global_dofs() << ".");
      vec_Ax[dirichlet_indices_[i]] = 0.;
    }

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Evaluate condensed matrix-vector product containing data.
   *
   * Function that evaluates the condensed, matrix-free version of the matrix-vector product
   * \f$A x - f = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method, \f$x\f$ is the
   * vector of parameters to define the skeletal variable \f$\lambda\f$, \f$f\f$ is the right-hand
   * side of the problem and \f$y\f$ is the resulting vector, which has the same size as the input
   * vector \f$x\f$.
   *
   * \tparam  hyNode_index_t  Typename of the hypernode index. Defaults to \c unsigned \c int.
   * \param   x_vec         A vector containing the input vector \f$x\f$.
   * \param   time          Time at which analytical functions will be evaluated.
   * \retval  y_vec         A vector containing the product \f$y = Ax\f$.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT residual_flux(const LargeVecT& x_vec, const dof_value_t time = 0.)
  {
    auto vec_Ax = prototype_mat_vec_multiply(residual_flux, has_residual_flux);

    // Set all Dirichlet values to zero.
    for (dof_index_t i = 0; i < dirichlet_indices_.size(); ++i)
    {
      hy_assert(dirichlet_indices_[i] >= 0 && dirichlet_indices_[i] < hyper_graph_.n_global_dofs(),
                "All indices of Dirichlet nodes need to be larger than or equal to zero and "
                  << "smaller than the total amount of degrees of freedom." << std::endl
                  << "In this case, the index is " << dirichlet_indices_[i] << " and the total "
                  << "amount of hypernodes is " << hyper_graph_.n_global_dofs() << ".");
      vec_Ax[dirichlet_indices_[i]] = 0.;
    }

    return vec_Ax;
  }
  /*!***********************************************************************************************
   * \brief   Calculate L2 error of approximated function.
   *
   * \tparam  hyNode_index_t  Typename of the hypernode index. Defaults to \c unsigned \c int.
   * \param   x_vec           A vector containing the input vector \f$x\f$.
   * \param   time            Time at which analytical functions will be evaluated.
   * \retval  error           A vector containing the errors.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  std::vector<dof_value_t> errors(const LargeVecT& x_vec, const dof_value_t time = 0.)
  {
    auto result = prototype_errors(errors, has_errors);
    return std::vector<dof_value_t>(result.begin(), result.end());
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
   * \retval  n               Size of the condensed (square) system of equations.
   ************************************************************************************************/
  dof_index_t size_of_system() const { return hyper_graph_.n_global_dofs(); }
  /*!***********************************************************************************************
   * \brief   Set plot option and return old plot option.
   *
   * Function to set and / or read the current plot option.
   *
   * \param   option          A \c std::string containing the plot option to be considered.
   * \param   value           A \c std::string containing the new value of the considered option.
   *                          If empty, the old value is kept.
   * \retval  opt_value       A \c std::string containing the value of the plot option.
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
   * \param   time          Time at which analytical functions are evaluated.
   * \retval  file          A file in the output directory.
   ************************************************************************************************/
  void plot_solution(const LargeVecT& lambda, const dof_value_t time = 0.)
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
};  // end of class Elliptic

}  // end of namespace GlobalLoop
