#include <iostream>
#include <cmath>

#include "AbstractProblem.h"

using namespace std;


int main(int argc, char *argv[])
{
  const unsigned int hyperedge_dim = 1, space_dim = 2;
  const unsigned int polynomial_degree = 1; //TODO Change to polynomial degree and deduce amount of ansatz functions
  vector<int> num_elements(2, 2);
   
  DiffusionProblemRegular<hyperedge_dim, space_dim, polynomial_degree> diffusion_problem(num_elements);
  vector<double> vectorDirichlet = diffusion_problem.return_zero_vector();
  vector<int> dirichlet_indices = { 0 , 8 };
  diffusion_problem.read_dirichlet_indices(dirichlet_indices);
  vector<double> result = diffusion_problem.matrix_vector_multiply(vectorDirichlet);
  
//  DiffusionProblemRegular<hyperedge_dim,space_dim> diffusion_problem(num_elements, polynomial_degree);
//  diffusion_problem.run();
  
//  DiffusionProblemRegular<hyperedge_dim,space_dim> test_problem(num_elements, polynomial_degree);
//  vector<double> vec = test_problem.return_zero_vector();
//  vec[0] = 1.;
//  vector<double> result = test_problem.matrix_vector_multiply(vec);
  
  for (unsigned int i = 0; i < result.size(); ++i)
    cout << "  " << result[i] << endl;
  
  return 0;
}

