/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "AbstractProblem.h"
#include "HyAssert.h"
#include <algorithm>
#include <array>

using namespace std;
#include "AbstractProblem.inst"


template <class TopologyT, class GeometryT, class LocalSolverT>
AbstractProblem<TopologyT,GeometryT,LocalSolverT>::
AbstractProblem(const typename TopologyT::constructor_value_type& construct_topo,
                const typename GeometryT::constructor_value_type& construct_geom,
                const typename LocalSolverT::constructor_value_type& construct_loc_sol)
: hyper_graph_  ( construct_topo, construct_geom ),
  local_solver_ ( construct_loc_sol ),
  plot_options  ( hyper_graph_, local_solver_ )
{
  static_assert( TopologyT::hyperedge_dimension() == GeometryT::hyperedge_dimension() ,
                 "Hyperedge dimension of topology and geometry must be equal!"  );
  static_assert( TopologyT::space_dimension() == GeometryT::space_dimension() ,
                 "Space dimension of topology and geometry must be equal!" );
  static_assert( TopologyT::hyperedge_dimension() == LocalSolverT::hyperedge_dimension() ,
                 "Hyperedge dimension of hypergraph and local solver must be equal!" );
}


/*
 * Function read_dirichlet_indices(..) receives a vector of signed integers, since it is part of the
 * code that interacts directly with the Python interface which does not allow for unsigned integers
 * and therefore includes a comparison between signed and unsigned integers in the assertion.
 */
template <class TopologyT, class GeometryT, class LocalSolverT>
void AbstractProblem<TopologyT,GeometryT,LocalSolverT>::
read_dirichlet_indices(std::vector<int> indices)
{
  dirichlet_indices_.resize(indices.size());
  for (unsigned int i = 0; i < indices.size(); ++i)
  {
    hy_assert( indices[i] >= 0 && indices[i] < hyper_graph_.num_of_hypernodes(), 
               "All indices of Dirichlet nodes need to be larger than or equal to zero and smaller "
               << "than the total amount of hypernodes." << endl << "In this case, the index is " <<
               indices[i] << " and the total amount of hypernodes is " << 
               hyper_graph_.num_of_hypernodes() << "." );
    dirichlet_indices_[i] = indices[i];
  }
}


template <class TopologyT, class GeometryT, class LocalSolverT>
vector<double> AbstractProblem<TopologyT,GeometryT,LocalSolverT>::
return_zero_vector()
{
  return vector<double>(hyper_graph_.num_of_global_dofs(), 0.);
}


template <class TopologyT, class GeometryT, class LocalSolverT>
vector<double> AbstractProblem<TopologyT,GeometryT,LocalSolverT>::
matrix_vector_multiply(vector<double> x_vec)
{
  constexpr unsigned int hyperedge_dim = TopologyT::hyperedge_dimension();
  constexpr unsigned int poly_degree = LocalSolverT::polynomial_degree();
  
  vector<double> vec_Ax(x_vec.size(), 0.);
  array<unsigned int, 2*hyperedge_dim> hyperedge_hypernodes;
  array< array<double, compute_n_dofs_per_node(hyperedge_dim, poly_degree)> , 2*hyperedge_dim >
    local_result, hyperedge_dofs;
  
  for_each( hyper_graph_.begin(), hyper_graph_.end(), [&](const auto hyperedge)
  {
    hyperedge_hypernodes = hyperedge.topology.get_hypernode_indices();
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyperedge_dofs[hypernode] = 
        hyper_graph_.hypernode_factory().get_dof_values(hyperedge_hypernodes[hypernode], x_vec);
    local_result = local_solver_.numerical_flux_from_lambda(hyperedge_dofs);
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyper_graph_.hypernode_factory().add_to_dof_values
        (hyperedge_hypernodes[hypernode], vec_Ax, local_result[hypernode]);
  });
  
  for(unsigned int i = 0; i < dirichlet_indices_.size(); ++i) 
    hyper_graph_.hypernode_factory().set_dof_values(dirichlet_indices_[i], vec_Ax, 0.);
    
  return vec_Ax;
}


template <class TopologyT, class GeometryT, class LocalSolverT>
int AbstractProblem<TopologyT,GeometryT,LocalSolverT>::
size_of_system()
{
  return hyper_graph_.num_of_global_dofs();
}


template <class TopologyT, class GeometryT, class LocalSolverT>
std::string AbstractProblem<TopologyT,GeometryT,LocalSolverT>::
plot_option(std::string option, std::string value)
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
  else if (option == "fileNumber")            value = to_string(plot_options.fileNumber);
  else if (option == "printFileNumber")       value = to_string(plot_options.printFileNumber);
  else if (option == "incrementFileNumber")   value = to_string(plot_options.incrementFileNumber);
  else hy_assert( 0 == 1 , "This plot option has not been defined (yet)." );
  
  return value;
}


template <class TopologyT, class GeometryT, class LocalSolverT>
void AbstractProblem<TopologyT,GeometryT,LocalSolverT>::
plot_solution(std::vector<double> lambda)
{
  plot(lambda, plot_options);
}
