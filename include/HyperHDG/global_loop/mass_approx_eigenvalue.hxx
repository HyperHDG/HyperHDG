#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/compile_time_tricks.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/plot.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <algorithm>
#include <array>
#include <cmath>


/*!*************************************************************************************************
 * \brief   This is an abstract example problem class.
 * 
 * This file provides an MassEigenvalue class defining an abstract problem. This abstract problem
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
class MassEigenvalue
{
  private:
    /*!*********************************************************************************************
     * \brief   Instantiation of a hypergraph.
     **********************************************************************************************/
    HDGHyperGraph
    < LocalSolverT::n_glob_dofs_per_node(), TopologyT, GeometryT, NodeDescriptorT, typename LocalSolverT::data_type >
    hyper_graph_;
    /*!*********************************************************************************************
     * \brief   Vector containing the indices of Dirichlet type nodes.
     **********************************************************************************************/
    std::vector<dof_index_t> dirichlet_indices_;
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
    MassEigenvalue
    ( 
      const typename TopologyT::constructor_value_type&    construct_topo,
      const typename GeometryT::constructor_value_type&    construct_geom,
      const typename LocalSolverT::constructor_value_type& construct_loc_sol
    )
    : hyper_graph_ ( construct_topo, construct_geom ), local_solver_ ( construct_loc_sol )
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
    MassEigenvalue
    ( 
      const typename TopologyT::constructor_value_type&    construct_topo,
      const typename LocalSolverT::constructor_value_type& construct_loc_sol
    )
    : hyper_graph_ ( construct_topo ), local_solver_ ( construct_loc_sol )
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
    MassEigenvalue ( const typename TopologyT::constructor_value_type& construct_topo )
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
     * \todo    ALL!
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
    template < typename hyNode_index_t = dof_index_t >
    std::vector<unsigned int> dirichlet_indices( )
    {
      constexpr unsigned int hyEdge_dim       = TopologyT::hyEdge_dim();
      std::array<dof_index_t, LocalSolverT::n_glob_dofs_per_node()> dof_indices;
      std::array<hyNode_index_t, 2*hyEdge_dim> hyEdge_hyNodes;
      std::array<unsigned int, 2 * hyEdge_dim> hyNode_types;
      
      std::for_each( hyper_graph_.begin() , hyper_graph_.end() , [&](auto hyper_edge)
      {
        hyNode_types = local_solver_.node_types(hyper_edge);
        hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();

        for (unsigned int i = 0; i < hyNode_types.size(); ++i)  if (hyNode_types[i] == 1)
        {
          dof_indices = hyper_graph_.hyNode_factory().get_dof_indices(hyEdge_hyNodes[i]);
          for (unsigned int j = 0; j < dof_indices.size(); ++j)
            dirichlet_indices_.push_back(dof_indices[j]);
        }
      });

      std::sort(dirichlet_indices_.begin(),dirichlet_indices_.end());
      auto last = std::unique(dirichlet_indices_.begin(),dirichlet_indices_.end());
      dirichlet_indices_.erase(last, dirichlet_indices_.end());

      return dirichlet_indices_;
    }
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
    std::vector<dof_value_t> matrix_vector_multiply
    ( const std::vector<dof_value_t>& x_vec, const dof_value_t time = 0. )
    {
      constexpr unsigned int hyEdge_dim       = TopologyT::hyEdge_dim();
      constexpr unsigned int n_dofs_per_node  = LocalSolverT::n_glob_dofs_per_node();
      
      std::vector<dof_value_t> vec_Ax( x_vec.size() , 0.);
      std::array<hyNode_index_t, 2*hyEdge_dim> hyEdge_hyNodes;
      std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;
      
      // Do matrix--vector multiplication by iterating over all hyperedges.
      std::for_each( hyper_graph_.begin() , hyper_graph_.end() , [&](auto hyper_edge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyEdge_dofs[hyNode] = 
            hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec);
        
        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr
        ( 
          not_uses_geometry
          < LocalSolverT,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>
            ( std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>& )
          >::value
        )
          hyEdge_dofs = local_solver_.numerical_flux_from_lambda(hyEdge_dofs, time);
        else  hyEdge_dofs = local_solver_.numerical_flux_from_lambda(hyEdge_dofs, hyper_edge, time);
        
        // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyper_graph_.hyNode_factory().add_to_dof_values
            (hyEdge_hyNodes[hyNode], vec_Ax, hyEdge_dofs[hyNode]);
      });
    
      return vec_Ax;
    }
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
    std::vector<dof_value_t> mass_matrix_multiply
    ( const std::vector<dof_value_t>& x_vec , const dof_value_t time = 0. )
    {
      constexpr unsigned int hyEdge_dim       = TopologyT::hyEdge_dim();
      constexpr unsigned int n_dofs_per_node  = LocalSolverT::n_glob_dofs_per_node();

      std::vector<dof_value_t> vec_Ax( x_vec.size() , 0.);
      std::array<hyNode_index_t, 2*hyEdge_dim> hyEdge_hyNodes;
      std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;

      // Do matrix--vector multiplication by iterating over all hyperedges.
      std::for_each( hyper_graph_.begin() , hyper_graph_.end() , [&](auto hyper_edge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyEdge_dofs[hyNode] =
            hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec);

        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr
        (
          not_uses_geometry
          < LocalSolverT,
            std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>
            ( std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * TopologyT::hyEdge_dim()>& )
          >::value
        )
        {
          hyEdge_dofs = local_solver_.numerical_flux_from_mass(hyEdge_dofs, time);
        }
        else
        {
          hyEdge_dofs = local_solver_.numerical_flux_from_mass(hyEdge_dofs, hyper_edge, time);
        }
        // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyper_graph_.hyNode_factory().add_to_dof_values
            (hyEdge_hyNodes[hyNode], vec_Ax, hyEdge_dofs[hyNode]);
      });
    
      return vec_Ax;
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
    std::string plot_option( const std::string& option, std::string value = "" )
    {
      return set_plot_option(plot_options, option, value);
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
    void plot_solution( const std::vector<dof_value_t>& lambda, const dof_value_t time = 0. )
    { plot(hyper_graph_, local_solver_, lambda , plot_options, time ); }
}; // end of class MassEigenvalue