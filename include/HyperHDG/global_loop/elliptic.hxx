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
          typename dof_index_t = unsigned int,
          typename param_time_t = double>
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
  Elliptic(const typename TopologyT::constructor_value_type& construct_topo,
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
  Elliptic(const typename TopologyT::constructor_value_type& construct_topo,
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
  Elliptic(const typename TopologyT::constructor_value_type& construct_topo)
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
  LargeVecT trace_to_flux(const LargeVecT& x_vec, const param_time_t time)
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
  
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT trace_to_flux(const LargeVecT& x_vec)
  {
    param_time_t time = 0;
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
  sparse_mat<LargeVecT> trace_to_flux_mat(const param_time_t time = 0)
  {
    return prototype_mat_generate(trace_to_flux, has_trace_to_flux);
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
  LargeVecT residual_flux(const LargeVecT& x_vec, const param_time_t time)
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
  
  template <typename hyNode_index_t = dof_index_t>
  LargeVecT residual_flux(const LargeVecT& x_vec)
  {
    param_time_t time = 0;
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
  std::vector<dof_value_t> errors(const LargeVecT& x_vec, const param_time_t time = 0)
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
  void plot_solution(const LargeVecT& lambda, const param_time_t time = 0)
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
