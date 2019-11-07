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
#include <cmath>
#include <iostream>

using namespace std;


template class DiffusionProblemRegular<1,1>;
template class DiffusionProblemRegular<1,2>;
template class DiffusionProblemRegular<1,3>;
template class DiffusionProblemRegular<2,2>;
template class DiffusionProblemRegular<2,3>;
template class DiffusionProblemRegular<3,3>;


template <unsigned int hyperedge_dim, unsigned int space_dim>
DiffusionProblemRegular<hyperedge_dim,space_dim>::
DiffusionProblemRegular(vector<int> num_elements, int polynomial_degree)
: hyper_graph_topology
  (
    JointGetter_RegularQuad<hyperedge_dim,space_dim>
      (pow((polynomial_degree + 1), hyperedge_dim-1), num_elements[0], num_elements[1], num_elements[2]),
    HyperGraph_Cubic< hyperedge_dim, space_dim >
      (num_elements[0], num_elements[1], num_elements[2])
  ) ,
  local_solver(polynomial_degree,2,1.)
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
  dirichlet_indices.push_back(0);
  dirichlet_indices.push_back(8);
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
vector<double> DiffusionProblemRegular<hyperedge_dim,space_dim>::
return_zero_vector()
{
  return vector<double>(hyper_graph_topology.num_of_global_dofs(), 0.);
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
vector<double> DiffusionProblemRegular<hyperedge_dim,space_dim>::
matrix_vector_multiply(vector<double> x_vec)
{
  vector<double> vec_Ax(x_vec.size(), 0.);
  vector< vector<double> > local_result, hyperedge_dofs;
  array<unsigned int, 2*hyperedge_dim> hyperedge_joints;
  
  for (unsigned int hyperedge = 0; hyperedge < hyper_graph_topology.num_of_hyperedges(); ++hyperedge)
  {
    hyperedge_joints = hyper_graph_topology.get_hyperedge(hyperedge).get_joint_indices();
    hyperedge_dofs.resize(hyperedge_joints.size());
    for (unsigned int joint = 0; joint < hyperedge_joints.size(); ++joint)
      hyperedge_dofs[joint] = hyper_graph_topology.get_joint(hyperedge_joints[joint]).get_dof_values(x_vec);
    local_result = local_solver.numerical_flux_from_lambda(hyperedge_dofs);
    for (unsigned int joint = 0; joint < hyperedge_joints.size(); ++joint)
      hyper_graph_topology.get_joint(hyperedge_joints[joint]).add_to_dof_values(vec_Ax, local_result[joint]);
  }
  
  for(unsigned int i = 0; i < dirichlet_indices.size(); ++i) 
    hyper_graph_topology.get_joint(dirichlet_indices[i]).set_dof_values(vec_Ax, 0.);
    
  return vec_Ax;
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
int DiffusionProblemRegular<hyperedge_dim,space_dim>::
size_of_system()
{
  return hyper_graph_topology.num_of_global_dofs();
}
