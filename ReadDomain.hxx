#ifndef READDOMAIN_HXX
#define READDOMAIN_HXX


#include "TypeDefs.hxx"
#include "Point.hxx"
#include "HyAssert.hxx"
#include "HyperCube.hxx"

#include <array>
#include <vector>

#include <fstream>
#include <string>
#include <sstream>



#include <iostream>

template < unsigned int hyEdge_dim, unsigned int space_dim >
struct DomainInfo
{
  std::vector< Point<space_dim> > points;
  std::vector< std::array< pt_index_t, 1 << (hyEdge_dim-1) > > points_hyNode;
  std::vector< std::array< hyNode_index_t, 2 * hyEdge_dim > > hyNodes_hyEdge;
  
  DomainInfo ( pt_index_t n_points, hyNode_index_t n_hyNodes, hyEdge_index_t n_hyEdges)
  : points(n_points), points_hyNode(n_hyNodes), hyNodes_hyEdge(n_hyEdges) { }
};






template < unsigned int hyEdge_dim, unsigned int space_dim >
DomainInfo<hyEdge_dim,space_dim> read_domain( const std::string& filename )
{
  hy_assert( filename.substr(filename.size()-4, filename.size()) == ".geo" ,
             "The given file needs to be a .geo file, since no other input file types are currently"
             << " implemented." );
  
  std::ifstream infile(filename);
  
  std::istringstream linestream;
  std::string line, keyword, equal_sign;
  unsigned int Space_Dim, HyperEdge_Dim;
  pt_index_t N_Points = 0, pt_iter;
  hyNode_index_t N_HyperNodes = 0, hyNode_iter;
  hyEdge_index_t N_HyperEdges = 0, hyEdge_iter;
  Point<space_dim> pt;
  
  while ( keyword != "Space_Dim" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword >> equal_sign >> Space_Dim;
  }
  
  hy_assert( keyword == "Space_Dim" ,
             "The keyword Space_Dim has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( Space_Dim == space_dim ,
             "Space_Dim in " << filename << " is " << Space_Dim << ", but should be "
             << space_dim << "!" );
  
  while ( keyword != "HyperEdge_Dim" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword >> equal_sign >> HyperEdge_Dim;
  }
  
  hy_assert( keyword == "HyperEdge_Dim" ,
             "The keyword HyperEdge_Dim has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( HyperEdge_Dim == hyEdge_dim ,
             "HyperEdge_Dim in " << filename << " is " << Space_Dim << ", but should be "
             << hyEdge_dim << "!" );
  
  while ( keyword != "N_Points" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword >> equal_sign >> N_Points;
  }
  
  hy_assert( keyword == "N_Points" ,
             "The keyword N_Points has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( N_Points != 0 ,
             "The value of N_Points has not been set correctly!" );
  
  while ( keyword != "N_HyperNodes" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword >> equal_sign >> N_HyperNodes;
  }
  
  hy_assert( keyword == "N_HyperNodes" ,
             "The keyword N_Points has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( N_HyperNodes != 0 ,
             "The value of N_HyperNodes has not been set correctly!" );
  
  while ( keyword != "N_HyperEdges" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword >> equal_sign >> N_HyperEdges;
  }
  
  hy_assert( keyword == "N_HyperEdges" ,
             "The keyword N_Points has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( N_HyperEdges != 0 ,
             "The value of N_HyperEdges has not been set correctly!" );
  
  DomainInfo<hyEdge_dim,space_dim> domain_info(N_Points, N_HyperNodes, N_HyperEdges);
  
  while ( keyword != "POINTS:" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  
  hy_assert( keyword == "POINTS:" ,
             "The keyword 'POINTS:' has not been found in the file " << filename << "!" );
  
  for ( pt_iter = 0; pt_iter < N_Points && std::getline(infile, line); ++pt_iter )
  {
    linestream = std::istringstream(line);
    for (unsigned int dim = 0; dim < space_dim; ++dim)  linestream >> pt[dim];
    domain_info.points[pt_iter] = pt;
  }
  
  hy_assert( pt_iter == N_Points ,
             "Not all points have been added to the list!" );
  
  while ( keyword != "HYPERNODES:" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  
  hy_assert( keyword == "HYPERNODES:" ,
             "The keyword 'HYPERNODES:' has not been found in the file " << filename << "!" );
  
  for ( hyNode_iter = 0; hyNode_iter < N_HyperNodes && std::getline(infile, line); ++hyNode_iter )
  {
    linestream = std::istringstream(line);
    for (unsigned int i = 0; i < domain_info.points_hyNode[hyNode_iter].size(); ++i)
      linestream >> domain_info.points_hyNode[hyNode_iter][i];
  }
  
  hy_assert( hyNode_iter == N_HyperNodes ,
             "Not all hypernodes have been added to the list!" );
  
  while ( keyword != "HYPEREDGES:" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  
  hy_assert( keyword == "HYPEREDGES:" ,
             "The keyword 'HYPEREDGES:' has not been found in the file " << filename << "!" );
  
  for ( hyEdge_iter = 0; hyEdge_iter < N_HyperEdges && std::getline(infile, line); ++hyEdge_iter )
  {
    linestream = std::istringstream(line);
    for (unsigned int i = 0; i < domain_info.hyNodes_hyEdge[hyEdge_iter].size(); ++i)
      linestream >> domain_info.hyNodes_hyEdge[hyEdge_iter][i];
  }
  
  hy_assert( hyEdge_iter == N_HyperEdges ,
             "Not all hyperedges have been added to the list!" );

  return domain_info;

}

#endif // end of ifndef READDOMAIN_HXX
