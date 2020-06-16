#include <HyperHDG/DenseLA.hxx>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

template < unsigned int dim, typename float_t = float, typename index_t = unsigned int >
void make_epsilon_neighborhood_graph ( std::string& filename )
{
  hy_assert( filename.substr(filename.size()-4, filename.size()) == ".pts" ,
             "The given file needs to be a .pts file for this function to be applicable!" );
  
  struct Pair
  {
    const index_t left, right;
    Pair(const index_t left_, const index_t right_) : left(left_), right(right_) { }
  };

  std::ifstream infile(filename);
  std::istringstream linestream;
  std::string line, keyword, equal_sign;

  unsigned int space_dim;
  float_t epsilon;
  index_t index;

  std::vector< Point<dim,float_t> > points;
  std::vector< index_t > search;
  std::vector< Pair > connections;
  Point<dim> pt;

  while ( keyword != "Space_Dim" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> space_dim;

  hy_assert( keyword == "Space_Dim" ,
             "The keyword Space_Dim has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( space_dim == dim ,
             "Space_Dim in " << filename << " is " << space_dim << ", but should be "
             << dim << "!" );

  while ( keyword != "Epsilon" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> epsilon;

  hy_assert( keyword == "Epsilon" ,
             "The keyword Space_Dim has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( epsilon > 0 ,
             "Epsilon needs to be larger than zeros, but is " << epsilon << "!" );

  while ( keyword != "POINTS:" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  
  hy_assert( keyword == "POINTS:" ,
             "The keyword 'POINTS:' has not been found in the file " << filename << "!" );

  while ( std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    for (unsigned int dimension = 0; dimension < dim; ++dimension)  linestream >> pt[dimension];
    points.push_back( pt );
  }

  std::sort( points.begin(), points.end() );
  hy_assert( std::adjacent_find( points.begin(), points.end() ) == points.end() , 
             "Points must be unique in given file!" );

  std::vector< bool > bool_vec(points.size(), false);
  search.push_back(0);
  bool_vec[0] = true;

  while( !search.empty() )
  {
    index = search.back();
    search.pop_back();
    for 
    ( 
      index_t ind = index;
      ind >= 0 && ind < points.size() && std::abs(points[index][0] - points[ind][0]) < epsilon; 
      --ind
    )
      if
      ( 
        norm_2(points[index] - points[ind]) < epsilon && ind != index
        && ( !bool_vec[ind] || std::find(search.begin(), search.end(), ind) != search.end() )
      )
      {
        search.push_back(ind);
        bool_vec[ind] = true;
        connections.push_back(Pair(index,ind));
      }
    for 
    ( 
      index_t ind = index;
      ind >= 0 && ind < points.size() && std::abs(points[index][0] - points[ind][0]) < epsilon; 
      ++ind
    )
      if
      ( 
        norm_2(points[index] - points[ind]) < epsilon && ind != index
        && ( !bool_vec[ind] || std::find(search.begin(), search.end(), ind) != search.end() )
      )
      {
        search.push_back(ind);
        bool_vec[ind] = true;
        connections.push_back(Pair(index,ind));
      }
  }

  infile.close();
  filename = filename + ".geo";
  std::ofstream outfile(filename); 

  outfile << "Space_Dim     = " << dim << ";" << std::endl;
  outfile << "HyperEdge_Dim = 1;" << std::endl << std::endl;
  outfile << "N_Points      = " << points.size() << ";" << std::endl;
  outfile << "N_HyperNodes  = " << points.size() << ";" << std::endl;
  outfile << "N_HyperEdges  = " << connections.size() << ";" << std::endl << std::endl;
  outfile << "POINTS:" << std::endl;
  for (unsigned int i = 0; i < points.size(); ++i)
    outfile << points[i];
  outfile << std::endl << "HYPERNODES_OF_HYPEREDGES:" << std::endl;
  for (unsigned int i = 0; i < connections.size(); ++i)
    outfile << connections[i].left << " " << connections[i].right << std::endl;
  outfile << "TYPES_OF_HYPERFACES:" << std::endl;
  for (unsigned int i = 0; i < connections.size(); ++i)
    outfile << connections[i].left << " " << connections[i].right << std::endl;
  outfile << "POINTS_OF_HYPEREDGES:" << std::endl;
  for (unsigned int i = 0; i < connections.size(); ++i)
    outfile << connections[i].left << " " << connections[i].right << std::endl;

  outfile.close();

  for (unsigned int i = 0; i < bool_vec.size(); ++i)
    hy_assert( bool_vec[i] , "All points need to belong to one graph!" );
}