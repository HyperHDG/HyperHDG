#include <tensor_rank1.hxx>
#include <iostream>

int main()
{
  Tensor::Rank1<2,3> t1;

  for (unsigned int i : { 0, 1, 2})
    {
      t1.factor(0)[i] = i+1;
      t1.factor(1)[i] = i+2;
    }
  for (unsigned int i : { 0, 1, 2})
    {
      for (unsigned int j : { 0, 1, 2})
	{
	  std::cout << t1({i,j}) << '\t';
	}
      std::cout << '\n';
    }
}
