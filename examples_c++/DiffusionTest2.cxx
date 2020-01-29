#include "../AbstractProblem.hxx"
#include "../SparseLinearAlgebra.hxx"
#include "../ReadDomain.hxx"
#include <iostream>

#include <string>

using namespace std;


int main(int argc, char *argv[])
{
  vector<int> num_elements = { 4 , 2 , 2 };
  
  DiffusionProblemRegularNaive<1,3,1> diffusion_problem(num_elements, num_elements, 1.);
  
  vector<double> vectorDirichlet = diffusion_problem.return_zero_vector();
  vectorDirichlet[0] = 1.;
  vectorDirichlet[vectorDirichlet.size()-1] = 0.;
  
  vector<int> index_vector = { 0 , ((int) vectorDirichlet.size())-1 };
  diffusion_problem.read_dirichlet_indices(index_vector);
  
  vector<double> vectorRHS = diffusion_problem.matrix_vector_multiply(vectorDirichlet);
  for (unsigned int i = 0; i < vectorRHS.size(); ++i)  vectorRHS[i] *= -1.;
  
  
  std::string filename = "domains/SimpleTriangle.geo";
  DomainInfo<1,2> di = read_domain<1,2>(filename);
  
  cout << "SUCCESS" << endl;
   
  
  
  return 0;
}
