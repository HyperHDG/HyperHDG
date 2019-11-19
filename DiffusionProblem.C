/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "DiffusionProblem.h"
//#include <cmath>
#include <iostream>
#include <array>
#include <cassert>

using namespace std;


template class DiffusionProblemRegular<1,1,1>;
template class DiffusionProblemRegular<1,2,1>;
template class DiffusionProblemRegular<1,3,1>;
//template class DiffusionProblemRegular<2,2,1>;
//template class DiffusionProblemRegular<2,3,1>;
//template class DiffusionProblemRegular<3,3,1>;


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
DiffusionProblemRegular<hyperedge_dim,space_dim,polynomial_degree>::
DiffusionProblemRegular(vector<int> num_elements)
: hyper_graph_topology(HyperGraph_Cubic< hyperedge_dim, space_dim >(num_elements)),
  local_solver(1.)
{
  cout << "Amount of HyperEdges = " << hyper_graph_topology.num_of_hyperedges() << endl;
  for(unsigned int i = 0; i < hyper_graph_topology.num_of_hyperedges(); ++i)
  {
    const HyperEdge_Cubic<hyperedge_dim,space_dim> hyperedge = hyper_graph_topology.get_hyperedge(i);
    const array<joint_index_type, 2*hyperedge_dim> indices = hyperedge.get_joint_indices();
    cout << i << "   ";
    for(unsigned int j = 0; j < indices.size(); ++j)  cout << indices[j] << "  ";
    cout << endl;
  }
}


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
void DiffusionProblemRegular<hyperedge_dim,space_dim,polynomial_degree>::
read_dirichlet_indices(std::vector<int> indices)
{
  dirichlet_indices.resize(indices.size());
  for (unsigned int i = 0; i < indices.size(); ++i)
  {
    assert( indices[i] >= 0 && indices[i] < hyper_graph_topology.num_of_vertices() );
    dirichlet_indices[i] = indices[i];
  }
}


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
vector<double> DiffusionProblemRegular<hyperedge_dim,space_dim,polynomial_degree>::
return_zero_vector()
{
  return vector<double>(hyper_graph_topology.num_of_global_dofs(), 0.);
}


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
vector<double> DiffusionProblemRegular<hyperedge_dim,space_dim,polynomial_degree>::
matrix_vector_multiply(vector<double> x_vec)
{
  vector<double> vec_Ax(x_vec.size(), 0.);
  array< array<double, local_dof_amount_node(hyperedge_dim, polynomial_degree)> , 2*hyperedge_dim > local_result, hyperedge_dofs;
  array<unsigned int, 2*hyperedge_dim> hyperedge_joints;
  
  for (unsigned int hyperedge = 0; hyperedge < hyper_graph_topology.num_of_hyperedges(); ++hyperedge)
  {
    hyperedge_joints = hyper_graph_topology.get_hyperedge(hyperedge).get_joint_indices();
    for (unsigned int joint = 0; joint < hyperedge_joints.size(); ++joint)
      hyperedge_dofs[joint] = hyper_graph_topology.hypernode_factory().get_dof_values(hyperedge_joints[joint], x_vec);
    local_result = local_solver.numerical_flux_from_lambda(hyperedge_dofs);
    for (unsigned int joint = 0; joint < hyperedge_joints.size(); ++joint)
      hyper_graph_topology.hypernode_factory().add_to_dof_values(hyperedge_joints[joint], vec_Ax, local_result[joint]);
  }
  
  for(unsigned int i = 0; i < dirichlet_indices.size(); ++i) 
    hyper_graph_topology.hypernode_factory().set_dof_values(dirichlet_indices[i], vec_Ax, 0.);
    
  return vec_Ax;
}


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
int DiffusionProblemRegular<hyperedge_dim,space_dim,polynomial_degree>::
size_of_system()
{
  return hyper_graph_topology.num_of_global_dofs();
}
