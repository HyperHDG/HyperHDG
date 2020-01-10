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
#include <algorithm>
#include <array>
#include <cassert>

using namespace std;
#include "AbstractProblem.inst"


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
DiffusionProblemRegular<hyperedge_dim,space_dim,poly_degree>::
DiffusionProblemRegular(vector<int> num_elements)
: hyper_graph_(num_elements,num_elements),
  local_solver_(1.), plot_options(hyper_graph_, local_solver_)
{ }


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
void DiffusionProblemRegular<hyperedge_dim,space_dim,poly_degree>::
read_dirichlet_indices(std::vector<int> indices)
{
  dirichlet_indices_.resize(indices.size());
  for (unsigned int i = 0; i < indices.size(); ++i)
  {
    assert( indices[i] >= 0 && indices[i] < hyper_graph_.num_of_hypernodes() );
    dirichlet_indices_[i] = indices[i];
  }
}


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
vector<double> DiffusionProblemRegular<hyperedge_dim,space_dim,poly_degree>::
return_zero_vector()
{
  return vector<double>(hyper_graph_.num_of_global_dofs(), 0.);
}


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
vector<double> DiffusionProblemRegular<hyperedge_dim,space_dim,poly_degree>::
matrix_vector_multiply(vector<double> x_vec)
{
  vector<double> vec_Ax(x_vec.size(), 0.);
  array< array<double, compute_n_dofs_per_node(hyperedge_dim, poly_degree)> , 2*hyperedge_dim > local_result, hyperedge_dofs;
  array<unsigned int, 2*hyperedge_dim> hyperedge_hypernodes;
  
  for_each( hyper_graph_.begin(), hyper_graph_.end(), [&](const auto hyperedge)
  {
    hyperedge_hypernodes = hyperedge.topology.get_hypernode_indices();
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyperedge_dofs[hypernode] = hyper_graph_.hypernode_factory().get_dof_values(hyperedge_hypernodes[hypernode], x_vec);
    local_result = local_solver_.numerical_flux_from_lambda(hyperedge_dofs);
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyper_graph_.hypernode_factory().add_to_dof_values(hyperedge_hypernodes[hypernode], vec_Ax, local_result[hypernode]);
  });
  
  for(unsigned int i = 0; i < dirichlet_indices_.size(); ++i) 
    hyper_graph_.hypernode_factory().set_dof_values(dirichlet_indices_[i], vec_Ax, 0.);
    
  return vec_Ax;
}


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
int DiffusionProblemRegular<hyperedge_dim,space_dim,poly_degree>::
size_of_system()
{
  return hyper_graph_.num_of_global_dofs();
}


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
std::string DiffusionProblemRegular<hyperedge_dim,space_dim,poly_degree>::
plot_option(std::string option, std::string value)
{
  if (value == "")                            ;
  else if (option == "outputDir")             plot_options.outputDir = value;
  else if (option == "fileName")              plot_options.fileName = value;
  else if (option == "fileEnding")            plot_options.fileEnding = value;
  else if (option == "fileNumber")            plot_options.fileNumber = stoi(value);
  else if (option == "printFileNumber")       plot_options.printFileNumber = (value == "true" || value == "1");
  else if (option == "incrementFileNumber")   plot_options.incrementFileNumber = (value == "true" || value == "1");
  else assert( 0 == 1 );
  
  if (option == "outputDir")                  value = plot_options.outputDir;
  else if (option == "fileName")              value = plot_options.fileName;
  else if (option == "fileEnding")            value = plot_options.fileEnding;
  else if (option == "fileNumber")            value = to_string(plot_options.fileNumber);
  else if (option == "printFileNumber")       value = to_string(plot_options.printFileNumber);
  else if (option == "incrementFileNumber")   value = to_string(plot_options.incrementFileNumber);
  else assert( 0 == 1 );
  
  return value;
}

template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int poly_degree>
void DiffusionProblemRegular<hyperedge_dim,space_dim,poly_degree>::
plot_solution(std::vector<double> lambda)
{
  plot(lambda, plot_options);
}
