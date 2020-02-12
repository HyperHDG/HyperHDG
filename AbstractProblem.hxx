/*!*************************************************************************************************
 * \file    AbstractProblem.hxx
 * \brief   This file provides the class AbstractProblem.
 *
 * This file provides a file \c AbstractProblem class defining an abstract problem.
 * 
 * This file is an .hxx file, since class and function are compiled "dynamically" depending on the
 * considered problem in Python or C++ code. Dynamically means, that either, when the C++ problem
 * or Python's ClassWrapper are compiled, the relevant template parameters for the respective class
 * and functions of this file are deduced and the needed versions are compiled.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#ifndef ABSTRACTPROBLEM_HXX
#define ABSTRACTPROBLEM_HXX

// Includes needed for external communication.
// These also would ben included when splitted in .C and .h files.
#include <HDGHyperGraph.hxx>
#include <Plot.hxx>
/// \todo Those do not have any business in "AbstractProblem.hxx"
#include <LSol_Diffusion.hxx>
#include <LSol_Elasticity.hxx>

// Includes solely needed for implementation of the different functions.
// These would not be included when splitted in .C and .h files.
#include <HyAssert.hxx>
#include <algorithm>
#include <array>

/*!*************************************************************************************************
 * \brief   This is an abstract example problem class.
 * 
 * \todo    Update these descriptions!
 *
 * This class contains functions to define and solve Poisson's equation, i.e.,
 * \f[
 *  - \Delta u = 0 \quad \text{ in } \Omega, \qquad u = u_\text D \quad \text{ on } \partial \Omega
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dim = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dim = 2) or volumes (\c hyEdge_dim = 3).
 *
 * \todo  The loop in matrix_vector_multiply() only combines properties of HyperGraph with local
 *        solvers, right? Dirichlet boundary conditions? Post filtering!
 *        -> I believe that we have to discuss, how to do this best. Note that the .C file now
 *        contains a for_each loop (cf. HDGHyperGraph.h)!
 * 
 * \todo  We should discuss, whether or not it makes sense to turn this class into an abstract class
 *        that receives a HyperGraph Topology, Geometry, and a LocalSolver as template parameters.
 * 
 * \todo  We should rewrite this explanation appropriately and think whether this is general enough.
 *        (With explanation, I mean this definition and the following function explanations, etc.)
 *
 * \tparam  hyEdge_dim    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dim     The dimension of the space, the object is located in. This number should
 *                        be larger than or equal to hyEdge_dim.
 * \tparam  poly_degree   The polynomial degree of test and trial functions.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template <class TopologyT, class GeometryT, class LocalSolverT>
class AbstractProblem
{
  private:
    HDGHyperGraph<LocalSolverT::n_glob_dofs_per_node(), TopologyT, GeometryT> hyper_graph_;
    std::vector<dof_index_type> dirichlet_indices_;
    const LocalSolverT  local_solver_;
    PlotOptions   plot_options;
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
    AbstractProblem( const typename TopologyT::constructor_value_type& construct_topo,
                     const typename GeometryT::constructor_value_type& construct_geom,
                     const typename LocalSolverT::constructor_value_type& construct_loc_sol )
    : hyper_graph_  ( construct_topo, construct_geom ),
      local_solver_ ( construct_loc_sol )
    {
      static_assert( TopologyT::hyEdge_dim() == GeometryT::hyEdge_dimension() ,
                     "Hyperedge dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::space_dim() == GeometryT::space_dimension() ,
                     "Space dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dimension() ,
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
    AbstractProblem( const typename TopologyT::constructor_value_type& construct_topo,
                     const typename LocalSolverT::constructor_value_type& construct_loc_sol )
    : hyper_graph_  ( construct_topo ),
      local_solver_ ( construct_loc_sol )
    {
      static_assert( TopologyT::hyEdge_dim() == GeometryT::hyEdge_dimension() ,
                     "Hyperedge dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::space_dim() == GeometryT::space_dimension() ,
                     "Space dimension of topology and geometry must be equal!" );
      static_assert( TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dimension() ,
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
        hy_assert( (dof_index_type) indices[i] >= 0
                      && (dof_index_type) indices[i] < hyper_graph_.n_global_dofs(), 
                   "All indices of Dirichlet nodes need to be larger than or equal to zero and "
                   << "smaller than the total amount of degrees of freedom." << std::endl
                   << "In this case, the index is " << indices[i] << " and the total amount of " <<
                   "hypernodes is " << hyper_graph_.n_global_dofs() << "." );
        dirichlet_indices_[i] = (dof_index_type) indices[i];
      }
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
    std::vector<dof_value_t> return_zero_vector( ) const
    {
      return std::vector<dof_value_t>(hyper_graph_.n_global_dofs(), 0.);
    }
    /*!*********************************************************************************************
     * \brief   Evaluate condensed matrix-vector product.
     * 
     * \todo    Decide whether \c auto \c hyperedge should be const!
     *
     * Function that evaluates the condensed, matrix-free version of the matrix-vector product
     * \f$A x = y\f$, where \f$A\f$ is the condensed matrix of the LDG-H method, \f$x\f$ is the
     * vector of parameters to define the skeletal variable \f$\lambda\f$, and \f$y\f$ is the
     * resulting vector, which has the same size as the input vector \f$x\f$.
     *
     * \param   x_vec         A \c std::vector containing the input vector \f$x\f$.
     * \retval  y_vec         A \c std::vector containing the product \f$y = Ax\f$.
     **********************************************************************************************/
    template < typename hyNode_index_t = unsigned int >
    std::vector<dof_value_t> matrix_vector_multiply( const std::vector<dof_value_t>& x_vec ) const
    {
      constexpr unsigned int hyEdge_dim       = TopologyT::hyEdge_dim();
      constexpr unsigned int n_dofs_per_node  = LocalSolverT::n_glob_dofs_per_node();
      
      std::vector<dof_value_t> vec_Ax( x_vec.size() , 0.);
      std::array<hyNode_index_t, 2*hyEdge_dim> hyEdge_hyNodes;
      std::array<std::array<dof_value_t, n_dofs_per_node>, 2 * hyEdge_dim> hyEdge_dofs;
      
      // Do matrix--vector multiplication by iterating over all hyperedges.
      std::for_each( hyper_graph_.begin() , hyper_graph_.end() , [&](const auto hyEdge)
      {
        // Fill x_vec's degrees of freedom of a hyperedge into hyEdge_dofs array.
        hyEdge_hyNodes = hyEdge.topology.get_hyNode_indices();
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyEdge_dofs[hyNode] = 
            hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[hyNode], x_vec);
        
        // Turn degrees of freedom of x_vec that have been stored locally into those of vec_Ax.
        if constexpr ( LocalSolverT::use_geometry() )
          hyEdge_dofs = local_solver_.numerical_flux_from_lambda(hyEdge_dofs, hyEdge.geometry);
        else  hyEdge_dofs = local_solver_.numerical_flux_from_lambda(hyEdge_dofs);
        
        // Fill hyEdge_dofs array degrees of freedom into vec_Ax.
        for ( unsigned int hyNode = 0 ; hyNode < hyEdge_hyNodes.size() ; ++hyNode )
          hyper_graph_.hyNode_factory().add_to_dof_values
            (hyEdge_hyNodes[hyNode], vec_Ax, hyEdge_dofs[hyNode]);
      });
      
      // Set all Dirichlet values to zero.
      for ( dof_index_type i = 0 ; i < dirichlet_indices_.size() ; ++i )
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
    dof_index_type size_of_system() const
    {
      return hyper_graph_.n_global_dofs();
    }
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
    void plot_solution( const std::vector<dof_value_t>& lambda ) const
    {
      plot(hyper_graph_, local_solver_, lambda , plot_options );
    }
}; // end of class AbstractProblem

/*!*************************************************************************************************
 * \brief   This is an example problem.
 *
 * This class contains functions to define and solve Poisson's equation, i.e.,
 * \f[
 *  - \Delta u = 0 \quad \text{ in } \Omega, \qquad u = u_\text D \quad \text{ on } \partial \Omega
 * \f]
 * in a spatial domain \f$\Omega \subset \mathbb R^d\f$. Here, \f$d\f$ is the spatial dimension
 * \c space_dim, \f$\Omega\f$ is a regular graph (\c hyEdge_dim = 1) or hypergraph whose
 * hyperedges are surfaces (\c hyEdge_dim = 2) or volumes (\c hyEdge_dim = 3).
 *
 * \tparam  hyEdge_dim    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dim     The dimension of the space, the object is located in. This number should
 *                        be larger than or equal to hyEdge_dim.
 * \tparam  poly_degree   The polynomial degree of test and trial functions.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dim, unsigned int space_dim, unsigned int poly_degree>
using DiffusionProblemRegularNaive = 
AbstractProblem < Topology::Cubic< hyEdge_dim, space_dim >,
                  Geometry::UnitCube< hyEdge_dim, space_dim >,
                  Diffusion_TensorialUniform < hyEdge_dim, poly_degree, 2 * poly_degree >
                >;

/*!*************************************************************************************************
 * \brief   This is an example problem.
 *
 * \todo    This has not yet been fully implemented!
 *
 * \tparam  hyEdge_dim    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dim     The dimension of the space, the object is located in. This number should
 *                        be larger than or equal to hyEdge_dim.
 * \tparam  poly_degree   The polynomial degree of test and trial functions.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
/*template <unsigned int hyEdge_dim, unsigned int space_dim, unsigned int poly_degree>
using ElasticityProblemRegular = 
AbstractProblem < Topology::Cubic< hyEdge_dim, space_dim >,
                  Geometry::UnitCube< hyEdge_dim, space_dim >,
                  ElasticitySolver_RegularQuad < hyEdge_dim, space_dim, poly_degree, 2*poly_degree >
                >;
*/
/*!*************************************************************************************************
 * \brief   This is an example problem.
 *
 * \todo    This has not yet been fully implemented!
 *
 * \tparam  hyEdge_dim    Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for
 *                        PDEs defined on surfaces, and 3 is for PDEs defined on volumes.
 * \tparam  space_dim     The dimension of the space, the object is located in. This number should
 *                        be larger than or equal to hyEdge_dim.
 * \tparam  poly_degree   The polynomial degree of test and trial functions.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2019--2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2019--2020.
 **************************************************************************************************/
template <unsigned int hyEdge_dim, unsigned int space_dim, unsigned int poly_degree>
using ElasticityProblemFile = 
AbstractProblem < Topology::File< hyEdge_dim, space_dim >,
                  Geometry::File< hyEdge_dim, space_dim >,
                  ElasticRods_TensorialUniform < hyEdge_dim, space_dim, poly_degree, 2*poly_degree >
                >;

#endif // end of ifndef ABSTRACTPROBLEM_HXX
