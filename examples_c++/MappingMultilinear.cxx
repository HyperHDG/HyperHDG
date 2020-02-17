
#include <tensor_mapping.hxx>
#include <Point.hxx>
#include <Hypercube.hxx>

template <unsigned int dim, unsigned int np>
void brick()
{
  std::cout << "brick< " << dim << " , " << np << " >\n";
  const unsigned int sdim = dim;
  std::array<Point<dim>, Hypercube<dim>::n_vertices()> vertices;

  for (unsigned int i=0;i<vertices.size();++i)
    {
      vertices[i][0] = -1. + 2*(i%2);
      if constexpr (dim >= 2) vertices[i][1] = -2. + 4*((i/2)%2);
      if constexpr (dim >= 3) vertices[i][2] = -3. + 6*((i/4)%2);
    }
  
  std::array<double, np> points;

  for (unsigned int i=0;i<points.size();++i)
    points[i] = 1.*i/(np-1);

  Tensor::MappingMultilinear<sdim,dim,np> mapping(vertices, points);

  for (unsigned int i=0;i<mapping.size();++i)
    {
      std::array<unsigned int, dim> indices;
      indices[0] = i%np;
      if constexpr (dim>=2) indices[1] = (i/np)%np;
      if constexpr (dim>=3) indices[2] = (i/(np*np))%np;
      
      if constexpr (dim>=2) if (i!= 0 && i%np == 0) std::cout << "\n";
      std::cout << "Index ";
      for (auto j: indices) std::cout << j;
      std::cout << "\tCoord " << mapping(indices) << "\n";
    }
}

int main()
{
  brick<1,5>();
  brick<2,5>();
  brick<3,3>();
}
