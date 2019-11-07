#include "HyperGraphTopology.h"
#include "DiffusionSolver.h"
#include <iostream>
#include <cmath>

#include "DiffusionProblem.h"

using namespace std;

/*
template <unsigned int hyperedge_dim, unsigned int space_dim>
class DiffusionProblemRegular
{
  private:
    HyperGraphTopology < JointGetter_RegularQuad<hyperedge_dim,space_dim>,
                         hyperedgeGetter_RegularQuad<hyperedge_dim,space_dim>,
                         Joint_RegularQuad,
                         hyperedge_RegularQuad<hyperedge_dim,space_dim> >
                       hyper_graph_topology;
    vector<double> matrix_vector_argument, matrix_vector_result;
    DiffusionSolver_RegularQuad<hyperedge_dim> local_solver;
    vector<unsigned int> Dirichlet_indices;
  public:
    DiffusionProblemRegular(const vector<unsigned int>& num_elements, const unsigned int polynomials);
    vector<double> return_zeros(vector<double> input);
/*    
    typedef MeshWorker::IntegrationInfo<dim> CellInfo;

    InteriorPenaltyProblem(const double _lambda, const double _mu, const unsigned int _num_row,
                           const unsigned int _num_column, const bool _homogenization,  const unsigned int _extent,
                           const unsigned int _communication_level, const bool _analytic_test);
    void run(unsigned int n_steps); */
//};

/*
template <unsigned int hyperedge_dim, unsigned int space_dim>
DiffusionProblemRegular<hyperedge_dim,space_dim>::
DiffusionProblemRegular(const vector<unsigned int>& num_elements, const unsigned int polynomial_degree)
: hyper_graph_topology
  (
    JointGetter_RegularQuad<hyperedge_dim,space_dim>
      (pow((hyperedge_dim + 1) * (polynomial_degree + 1), hyperedge_dim), num_elements[0], num_elements[1], num_elements[2]),
    hyperedgeGetter_RegularQuad<hyperedge_dim,space_dim>
      (num_elements[0], num_elements[1], num_elements[2])
  ) ,
  matrix_vector_argument(hyper_graph_topology.num_of_global_dofs(), 0.),
  matrix_vector_result(hyper_graph_topology.num_of_global_dofs(), 0.),
  local_solver(polynomial_degree,1,1.)
{
  cout << "Amount of hyperedges = " << hyper_graph_topology.num_of_hyperedges() << endl;
  for(unsigned int i = 0; i < hyper_graph_topology.num_of_hyperedges(); ++i)
  {
    const hyperedge_RegularQuad<hyperedge_dim,space_dim> hyperedge = hyper_graph_topology.get_hyperedge(i);
    const vector<joint_index_type> indices = hyperedge.get_joint_indices();
    cout << i << "   ";
    for(unsigned int j = 0; j < indices.size(); ++j)  cout << indices[j] << "  ";
    cout << endl;
  }
}


template <unsigned int hyperedge_dim, unsigned int space_dim>
vector<double> DiffusionProblemRegular<hyperedge_dim,space_dim>::
return_zeros(vector<double> input)
{
  return vector<double>(input.size(), 0.);
}

*/


int main(int argc, char *argv[])
{
  const unsigned int hyperedge_dim = 1, space_dim = 2;
  const unsigned int polynomial_degree = 1; //TODO Change to polynomial degree and deduce amount of ansatz functions
  vector<int> num_elements(3, 2);
   
//  DiffusionProblemRegular<hyperedge_dim,space_dim> diffusion_problem(num_elements, polynomial_degree);
//  diffusion_problem.run();
  
  DiffusionProblemRegular<hyperedge_dim,space_dim> test_problem(num_elements, polynomial_degree);
  vector<double> vec = test_problem.return_zero_vector();
  vec[0] = 1.;
  vector<double> result = test_problem.matrix_vector_multiply(vec);
  
  for (unsigned int i = 0; i < result.size(); ++i)
    cout << vec[i] << "  " << result[i] << endl;
  
  return 0;
}

