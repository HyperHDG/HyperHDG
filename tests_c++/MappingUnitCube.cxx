
#include <Geom_UnitCube.hxx>

template <unsigned int edim, unsigned int sdim, unsigned int npts>
void test()
{
  std::cout << "test< " << edim << " , " << sdim << " , " << npts << ">\n";
  std::vector<unsigned int> sizes (sdim);
  for (unsigned int d=0;d<sdim;++d)
    sizes[d] = 1<<d;

  Geometry::UnitCube<edim,sdim> geometry(sizes);

  unsigned int n=1;
  if (edim == sdim)
    for (unsigned int d=0;d<sdim;++d)
      n *= sizes[d];

  std::array<float, npts> p1d;
  for (unsigned int i=0;i<npts;++i)
    p1d[i] = 1.*i/(npts-1);
  
  for (unsigned int i=0;i<n;++i)
    {
      auto edge = geometry[i];
      auto mapping = edge.mapping_tensor(p1d);

      std::cout << "Cell " << i << "\n";
      for (unsigned int j=0;j<mapping.size();++j)
	std::cout << mapping.lexicographic(j) << "\n";
    }
  std::cout << "\n";
}

int main()
{
  test<2,2,2> ();
  test<2,2,3> ();
  test<3,3,2> ();
}
