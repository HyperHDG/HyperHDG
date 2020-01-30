#include "../AbstractProblem.hxx"
#include "../SparseLinearAlgebra.hxx"
#include "../ReadDomain.hxx"
#include <iostream>

#include <string>

using namespace std;
using namespace SparseLA;


int main(int argc, char *argv[])
{
  vector<int> num_elements = { 4 , 2 , 2 };
  
  std::string filename = "domains/SimpleTriangle.geo";
  DiffusionProblemFileNoGeo<1,2,1> diffusion_problem(filename, 1.);
  
  vector<double> vectorDirichlet = diffusion_problem.return_zero_vector();
  vectorDirichlet[0] = 1.;
  vectorDirichlet[vectorDirichlet.size()-1] = 0.;
  
  vector<int> index_vector = { 0 , ((int) vectorDirichlet.size())-1 };
  diffusion_problem.read_dirichlet_indices(index_vector);
  
  vector<double> vectorRHS = diffusion_problem.matrix_vector_multiply(vectorDirichlet);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)  vectorRHS[i] *= -1.;
  
  int num_of_iterations = 0;
  vector<double> solution = conjugate_gradient( vectorRHS, diffusion_problem, num_of_iterations );
  solution = linear_combination(1., solution, 1., vectorDirichlet);
  
  for (unsigned int i = 0; i < solution.size(); ++i)  cout << solution[i] << endl;
  
  cout << "SUCCESS" << endl;
   
  
  
  return 0;
}
