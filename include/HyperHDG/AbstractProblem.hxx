#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/CompileTimeTricks.hxx>
#include <HyperHDG/HDGHyperGraph.hxx>
#include <HyperHDG/Plot.hxx>
#include <HyperHDG/HyAssert.hxx>
#include <algorithm>
#include <array>


/*!*************************************************************************************************
 * \bief  Check in reduced complexity whether an element is contained in a \c std::vector.
 **************************************************************************************************/
// template<class T>
// bool contains(const std::vector<T>& container, const T& v)
// {
//   auto it = std::lower_bound( container.begin(), container.end(), v);
//   return (it != container.end() && *it == v);
// }
/*!*************************************************************************************************
 * \brief   This is an abstract example problem class.
 * 
 * This file provides an AbstractProblem class defining an abstract problem. This abstract problem
 * serves as one possible, simple interface to Python. At the moment, it can be used to quickly
 * prototype testcases and others.
 *
 * \todo  The loop in matrix_vector_multiply() only combines properties of HyperGraph with local
 *        solvers, right? Dirichlet boundary conditions? Post filtering!
 *        -> A: I believe that we have to discuss, how to do this best. Note that the now, there is
 *        a for_each loop (cf. HDGHyperGraph.hxx)!
 * 
 * \todo  We should discuss, whether or not it makes sense to turn this class into an abstract class
 *        that receives a HyperGraph Topology, Geometry, and a LocalSolver as template parameters.
 *        -> A: This is the case already. I do not really see the difference.
 * 
 * \todo  We should rewrite this explanation appropriately and think whether this is general enough.
 *        (With explanation, I mean this definition and the following function explanations, etc.)
 *        -> A: Agreed.
 *        
 * \tparam  TopologyT     Class type containing topological information.
 * \tparam  GeometryT     Class type containing geometrical information.
 * \tparam  LocalSolverT  Class type of the local solver.
 * \tparam  dof_index_t   Index type of hyperedges. Default is \c unsigned \c int.
 *
 * \authors   Guido Kanschat, Heidelberg University, 2019--2020.
 * \authors   Andreas Rupp, Heidelberg University, 2019--2020.
 **************************************************************************************************/
template
< 
  class TopologyT, class GeometryT, class NodeDescriptorT, class LocalSolverT, 
  typename dof_index_t = unsigned int
>
class AbstractProblem
{
  private:
    /*!*********************************************************************************************
     * \brief   Instantiation of a hypergraph.
     **********************************************************************************************/
    HDGHyperGraph < LocalSolverT::n_glob_dofs_per_node(), TopologyT, GeometryT, NodeDescriptorT >
    hyper_graph_;
    /*!*********************************************************************************************
     * \brief   Vector containing the indices of Dirichlet type nodes.
     **********************************************************************************************/
    std::vector<dof_index_t> dirichlet_indices_, neumann_indices_;
    /*!*********************************************************************************************
     * \brief   Instantiation of a local solver.
     **********************************************************************************************/
    const LocalSolverT local_solver_;
    /*!*********************************************************************************************
     * \brief   Struct encoding the options for plotting.
     **********************************************************************************************/
    PlotOptions plot_options;
  public:
    /*!*********************************************************************************************
     * \brief   Abstract problem constructor.
     *
     * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
     * hyperedge.
     *
     * \param   construct_topo    Information to construct a topology.
     * \param   construct_geom    Information to construct a geometry.
     * \param   construct_loc_sol Information to construct a local solver.
     **********************************************************************************************/
    AbstractProblem
    ( 
      const typename TopologyT::constructor_value_type& construct_topo,
      const typename GeometryT::constructor_value_type& construct_geom,
      const typename LocalSolverT::constructor_value_type& construct_loc_sol
    )
    : hyper_graph_  ( construct_topo, construct_geom ), local_solver_ ( construct_loc_sol )
    {
      static_assert( TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim() ,
                     "Hyperedge dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::space_dim() == GeometryT::space_dim() ,
                     "Space dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim() ,
                     "Hyperedge dimension of hypergraph and local solver must be equal!" );
    }
    /*!*********************************************************************************************
     * \brief   Abstract problem constructor.
     *
     * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
     * hyperedge.
     *
     * \param   construct_topo    Information to construct a topology.
     * \param   construct_loc_sol Information to construct a local solver.
     **********************************************************************************************/
    AbstractProblem
    ( 
      const typename TopologyT::constructor_value_type& construct_topo,
      const typename LocalSolverT::constructor_value_type& construct_loc_sol
    )
    : hyper_graph_  ( construct_topo ), local_solver_ ( construct_loc_sol )
    {
      static_assert( TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim() ,
                     "Hyperedge dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::space_dim() == GeometryT::space_dim() ,
                     "Space dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim() ,
                     "Hyperedge dimension of hypergraph and local solver must be equal!" );
    }
    /*!*********************************************************************************************
     * \brief   Abstract problem constructor.
     *
     * Constructor for class containing a HyperGraph and a local solver that solve a PDE on a
     * hyperedge.
     *
     * \param   construct_topo    Information to construct a topology.
     **********************************************************************************************/
    AbstractProblem ( const typename TopologyT::constructor_value_type& construct_topo )
    : hyper_graph_  ( construct_topo )
    {
      static_assert( TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim() ,
                     "Hyperedge dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::space_dim() == GeometryT::space_dim() ,
                     "Space dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim() ,
                     "Hyperedge dimension of hypergraph and local solver must be equal!" );
    }
    /*!*********************************************************************************************
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
     * \param   indices       A \c std::vector containing the (global) indices of Dirichlet type
     *                        hypernodes/faces.
     **********************************************************************************************/
    void read_dirichlet_indices( const std::vector<unsigned int>& indices )
    {
      dirichlet_indices_.resize(indices.size());
      for (unsigned int i = 0; i < indices.size(); ++i)
      {
        hy_assert( (dof_index_t) indices[i] >= 0
                      && (dof_index_t) indices[i] < hyper_graph_.n_global_dofs(), 
                   "All indices of Dirichlet nodes need to be larger than or equal to zero and "
                   << "smaller than the total amount of degrees of freedom." << std::endl
                   << "In this case, the index is " << indices[i] << " and the total amount of " <<
                   "hypernodes is " << hyper_graph_.n_global_dofs() << "." );
        dirichlet_indices_[i] = (dof_index_t) indices[i];
      }
      std::sort(dirichlet_indices_.begin(),dirichlet_indices_.end());
    }
    /*!*********************************************************************************************
     * \brief   Read indices of Neumann type hypernodes/faces.
     *
     * Read the indices ot the hypernodes/faces that are of Neumann type and therefore do not
     * contain degrees of freedom that are allowed to change within iterations of the iterative
     * solver and other processes. In contrast, these degrees of freedom are set by the user.
     *
     * The user creates a vector that contains the coefficients of the corresponding degrees of
     * freedom (read by this function) and defines the Dirichlet values by this choice. The
     * remaining elements of the global vector of unknowns (which is \b not the vector \c indices
     * are supposed to be zero).
     *
     * \param   indices       A \c std::vector containing the (global) indices of Dirichlet type
     *                        hypernodes/faces.
     **********************************************************************************************/
    void read_neumann_indices( const std::vector<unsigned int>& indices )
    {
      neumann_indices_.resize(indices.size());
      for (unsigned int i = 0; i < indices.size(); ++i)
      {
        hy_assert( (dof_index_t) indices[i] >= 0
                      && (dof_index_t) indices[i] < hyper_graph_.n_global_dofs(), 
                   "All indices of Dirichlet nodes need to be larger than or equal to zero and "
                   << "smaller than the total amount of degrees of freedom." << std::endl
                   << "In this case, the index is " << indices[i] << " and the total amount of " <<
                   "hypernodes is " << hyper_graph_.n_global_dofs() << "." );
        neumann_indices_[i] = (dof_index_t) indices[i];
      }
      std::sort(neumann_indices_.begin(),neumann_indices_.end());
    }
    /*!*********************************************************************************************
     * \brief   Returns vector of appropriate size for the predefined problem.
     *
     * Returns a vector containing only the value zero, but of the size \f$n\f$ which is also the
     * number which is returned if \c size_of_system() is evaluated.
     *
     * \retval  zero          A \c std::vector of the correct size for the unknowns of the given
     *                        problem.
     **********************************************************************************************/
    template < typename dof_value_t = double >  std::vector<dof_value_t> return_zero_vector () const
    { return std::vector<dof_value_t>(hyper_graph_.n_global_dofs(), 0.); }
    /*!*********************************************************************************************
     * \brief   Evaluate condensed matrix-vector product.
     *
     * Function that evaluates the condensed, matrix-free version of the matrix-vector product
     * \f$A x = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method, \f$x\f$ is the
     * vector of parameters to define the skeletal variable \f$\lambda\f$, and \f$y\f$ is the
     * resulting vector, which has the same size as the input vector \f$x\f$.
     *
     * \param   x_vec         A \c std::vector containing the input vector \f$x\f$.
     * \retval  y_vec         A \c std::vector containing the product \f$y = Ax\f$.
     **********************************************************************************************/
    template < typename hyNode_index_t = dof_index_t, typename dof_value_t >
    std::vector<dof_value_t> matrix_vector_multiply ( const std::vector<dof_value_t>& x_vec ) const
    {
      constexpr unsigned int hyEdge_dim       = TopologyT::hyEdge_dim();
      constexpr unsigned int n_dofs_per_node  = LocalSolverT::n_glob_dofs_per_node();
      
      std::vector<dof_value_t> vec_Ax( x_vec.size() , 0.);
      std::array<hyNode_index_t, 2*hyEdge_dim> hyEdge_hyNodes;
      std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;
      
      // Do matrix--vector multiplication by iterating over all hyperedges.
      std::for_each( hyper_graph_.begin() , hyper_graph_.end() , [&](auto hyEdge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyEdge.topology.get_hyNode_indices();
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyEdge_dofs[hyNode] = 
            hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec);
        
        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr
        ( 
          not_uses_geometry
          < LocalSolverT,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>
            ( std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim>& )
          >::value
        )
          hyEdge_dofs = local_solver_.numerical_flux_from_lambda(hyEdge_dofs);
        else  hyEdge_dofs = local_solver_.numerical_flux_from_lambda(hyEdge_dofs, hyEdge.geometry);
        
        // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyper_graph_.hyNode_factory().add_to_dof_values
            (hyEdge_hyNodes[hyNode], vec_Ax, hyEdge_dofs[hyNode]);
      });
      
      // Set all Dirichlet values to zero.
      for ( dof_index_t i = 0 ; i < dirichlet_indices_.size() ; ++i )
      {
        hy_assert( dirichlet_indices_[i] >= 0 
                     && dirichlet_indices_[i] < hyper_graph_.n_global_dofs() ,
                   "All indices of Dirichlet nodes need to be larger than or equal to zero and "
                   << "smaller than the total amount of degrees of freedom." << std::endl
                   << "In this case, the index is " << dirichlet_indices_[i] << " and the total " <<
                   "amount of hypernodes is " << hyper_graph_.n_global_dofs() << "." );
        vec_Ax[dirichlet_indices_[i]] = 0.;
      }
    
      return vec_Ax;
    }
    /*!*********************************************************************************************
     * \brief   Calculate errors.
     *
     * Function that evaluates the global right-hand side (implemented wthin the local solver) and
     * adds the result to the function argument.
     *
     * \param   x_vec         A \c std::vector containing the input vector \f$x\f$.
     * \retval  error         A \c std::vector containing the errors.
     **********************************************************************************************/
    template < typename hyNode_index_t = dof_index_t, typename dof_value_t >
    std::vector<dof_value_t> calculate_errors ( const std::vector<dof_value_t>& x_vec ) const
    {
      if constexpr ( !LocalSolverT::advanced_functions() )
      {
        hy_assert( false , "This functionality is not implemented for the local solver!" );
        return x_vec;
      }

      constexpr unsigned int hyEdge_dim       = TopologyT::hyEdge_dim();
      constexpr unsigned int n_dofs_per_node  = LocalSolverT::n_glob_dofs_per_node();
      constexpr unsigned int n_errors         = LocalSolverT::n_errors();
      
      std::vector<dof_value_t> result(n_errors, 0.);

      if constexpr ( n_errors == 0 )  return result;

      std::array<dof_value_t, n_errors> loc_error;
      std::array<hyNode_index_t, 2*hyEdge_dim> hyEdge_hyNodes;
      std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;
      
      // Calculate errors by iteration over all hyperedges.
      std::for_each( hyper_graph_.begin() , hyper_graph_.end() , [&](auto hyEdge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyEdge.topology.get_hyNode_indices();
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyEdge_dofs[hyNode] = 
            hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec);
        
        // Turn degrees of freedom of x_vec that have been stored locally into local errors.
        if constexpr ( LocalSolverT::use_geometry() )
          loc_error = local_solver_.calc_loc_error(hyEdge_dofs,hyEdge.geometry);
        else  loc_error = local_solver_.calc_loc_error(hyEdge_dofs);
        
        // Fill the vector of errors.
        for ( unsigned int err = 0 ; err < n_errors ; ++err )  result[err] += loc_error[err];
      });
    }
    /*!*********************************************************************************************
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
     **********************************************************************************************/
    dof_index_t size_of_system() const
    { return hyper_graph_.n_global_dofs(); }
    /*!*********************************************************************************************
     * \brief   Set plot option and return old plot option.
     *
     * Function to set and / or read the current plot option.
     *
     * \param   option        A \c std::string containing the plot option to be considered.
     * \param   value         A \c std::string containing the new value of the considered option.
     *                        If empty, the old value is kept.
     * \retval  opt_value     A \c std::string containing the value of the plot option.
     **********************************************************************************************/
    std::string& plot_option( const std::string& option, std::string& value = "" )
    {
      if (value == "")                            ;
      else if (option == "outputDir")             plot_options.outputDir = value;
      else if (option == "fileName")              plot_options.fileName = value;
      else if (option == "fileEnding")            plot_options.fileEnding = value;
      else if (option == "fileNumber")            plot_options.fileNumber = stoi(value);
      else if (option == "printFileNumber")       plot_options.printFileNumber =
                                                    (value == "true" || value == "1");
      else if (option == "incrementFileNumber")   plot_options.incrementFileNumber =
                                                    (value == "true" || value == "1");
      else hy_assert( 0 == 1 , "This plot option has not been defined (yet)." );
  
      if (option == "outputDir")                  value = plot_options.outputDir;
      else if (option == "fileName")              value = plot_options.fileName;
      else if (option == "fileEnding")            value = plot_options.fileEnding;
      else if (option == "fileNumber")            value = std::to_string(plot_options.fileNumber);
      else if (option == "printFileNumber")       value = std::to_string
                                                            (plot_options.printFileNumber);
      else if (option == "incrementFileNumber")   value = std::to_string
                                                            (plot_options.incrementFileNumber);
      else hy_assert( 0 == 1 , "This plot option has not been defined (yet)." );
  
      return value;
    }
    /*!*********************************************************************************************
     * \brief   Plot solution in vtu format.
     *
     * Function that plots the solution of the problem to a predefined file.
     *
     * \param   lambda        A vector of unknowns containing the data vector.
     * \retval  file          A file in the output directory.
     **********************************************************************************************/
    template < typename dof_value_t >
    void plot_solution( const std::vector<dof_value_t>& lambda ) const
    { plot(hyper_graph_, local_solver_, lambda , plot_options ); }
}; // end of class AbstractProblem
