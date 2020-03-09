/*!*************************************************************************************************
 * \file    ReadDomain.hxx
 * \brief   This file reads domains from different file types and returns the uniform DomainInfo.
 *
 * \todo    DomainInfo assumes that hyperedges have 2 * hyEdge_dim hyperedges and that every
 *          hyperedge has 2 ^ (hyEdge_dim) points/vertices. This does neither allow for triangular
 *          hypernodes (needed for simplicial hyperedges in 3D) nor for tringular hyperedges (2D)!
 *          Is this ok or should this be generalized?
 * 
 * \todo    Do the remaining doxygen for this file (when corresponding topology and geometry are
 *          available) and ensure compatability with the rest of the code.
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#pragma once // Ensure that file is included only once in a single compilation.

#include <HyperHDG/Point.hxx>
#include <HyperHDG/HyAssert.hxx>

#include <algorithm>
#include <array>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>


template < class T >
bool is_unique(const std::vector<T>& vec)
{
  std::vector<T> x(vec);                                          // Copy vector enabling reorder.
  std::sort( x.begin(), x.end() );                                // Sort vector in O(N log N).
  return ( std::adjacent_find( x.begin(), x.end() ) == x.end() ); // Check for duplicates in O(n).
}


template 
< unsigned int hyEdge_dim, unsigned int space_dim, typename hyEdge_index_t = unsigned int,
  typename hyNode_index_t = unsigned int, typename pt_index_t = unsigned int >
struct DomainInfo
{
  hyEdge_index_t                                              n_hyEdges;
  hyNode_index_t                                              n_hyNodes;
  pt_index_t                                                  n_points;
  std::vector< Point<space_dim> >                             points;
  std::vector< std::array< hyNode_index_t, 2 * hyEdge_dim > > hyNodes_hyEdge; // 2 * hyEdge_dim
  std::vector< std::array< pt_index_t, 1 << hyEdge_dim > >    points_hyEdge;  // 2 ^ hyEdge_dim
  
  DomainInfo (const pt_index_t n_points, const hyEdge_index_t n_hyEdge,
              const hyNode_index_t n_hyNode, const pt_index_t n_point)
  : n_hyEdges(n_hyEdge), n_hyNodes(n_hyNode), n_points(n_point),
    points(n_points), hyNodes_hyEdge(n_hyEdges), points_hyEdge(n_hyEdges) { }
  
  bool check_consistency()
  {
    bool consistent = ( hyNodes_hyEdge.size() == hyNodes_hyEdge.size() );
    hy_assert( consistent ,
               "The sizes of hyNodes_hyEdge and hyNodes_hyEdge need to be equal!" );
    
    std::for_each( hyNodes_hyEdge.begin(), hyNodes_hyEdge.end(),
                   [&]( std::array< hyNode_index_t, 2 * hyEdge_dim > hyEdge )
    {
      for (unsigned int i = 0; i < hyEdge.size(); ++i)
      {
        consistent = ( hyEdge[i] < n_hyNodes && hyEdge[i] >= 0 );
        hy_assert( consistent ,
                   "At least one hypernode index is invalid!" );
       }
    });
    
    std::for_each( points_hyEdge.begin(), points_hyEdge.end(),
                   [&]( std::array< pt_index_t, 1 << hyEdge_dim > hyEdge )
    {
      for (unsigned int i = 0; i < hyEdge.size(); ++i)
      {
        consistent = ( hyEdge[i] < n_points && hyEdge[i] >= 0 );
        hy_assert( consistent ,
                   "At least one point index is invalid!" );
      }
    });
    
    consistent = is_unique(points);
    hy_assert( consistent ,
               "DomainInfo.points contains duplicate points!" );
    if ( !consistent ) return false;
    
    consistent = is_unique(hyNodes_hyEdge);
    hy_assert( consistent ,
               "DomainInfo.hyNodes_hyEdge contains duplicate hypernode!" );
    if ( !consistent ) return false;
    
    consistent = is_unique(points_hyEdge);
    hy_assert( consistent ,
               "DomainInfo.points_hyEdge contains duplicate points!" );
    if ( !consistent ) return false;
    
    return true;
  } // end of check_consistency
}; // end of struct DomainInfo


template
< unsigned int hyEdge_dim, unsigned int space_dim, typename hyEdge_index_t = unsigned int,
  typename hyNode_index_t = unsigned int, typename pt_index_t = unsigned int >
DomainInfo<hyEdge_dim,space_dim> read_domain_geo( const std::string& filename )
{
  hy_assert( filename.substr(filename.size()-4, filename.size()) == ".geo" ,
             "The given file needs to be a .geo file for this function to be applicable!" );
  
  std::ifstream infile(filename);
  std::istringstream linestream;
  std::string line, keyword, equal_sign;
  
  unsigned int    Space_Dim, HyperEdge_Dim;
  pt_index_t      N_Points = 0, pt_iter;
  hyEdge_index_t  N_HyperEdges = 0, hyEdge_iter;
  hyNode_index_t  N_HyperNodes = 0;
  Point<space_dim> pt;
  
  while ( keyword != "Space_Dim" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> Space_Dim;
  
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
    linestream >> keyword;
  }
  linestream >> equal_sign >> HyperEdge_Dim;
  
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
    linestream >> keyword;
  }
  linestream >> equal_sign >> N_Points;
  
  hy_assert( keyword == "N_Points" ,
             "The keyword N_Points has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( N_Points != 0 ,
             "The value of N_Points has not been set correctly!" );
  
  while ( keyword != "N_HyperNodes" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> N_HyperNodes;
  
  hy_assert( keyword == "N_HyperNodes" ,
             "The keyword N_Points has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( N_HyperNodes != 0 ,
             "The value of N_HyperNodes has not been set correctly!" );
  
  while ( keyword != "N_HyperEdges" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  linestream >> equal_sign >> N_HyperEdges;
  
  hy_assert( keyword == "N_HyperEdges" ,
             "The keyword N_Points has not been found in the file " << filename << "!" );
  hy_assert( equal_sign == "=" ,
             "The keyword " << keyword << " has not been followd by = symbol!" );
  hy_assert( N_HyperEdges != 0 ,
             "The value of N_HyperEdges has not been set correctly!" );
  
  DomainInfo<hyEdge_dim,space_dim> domain_info(N_Points, N_HyperEdges, N_HyperNodes, N_Points);
  
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
  
  while ( keyword != "HYPERNODES_OF_HYPEREDGES:" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  
  hy_assert( keyword == "HYPERNODES_OF_HYPEREDGES:" ,
             "The keyword 'HYPERNODES_OF_HYPEREDGES:' has not been found in the file "
             << filename << "!" );
  
  for ( hyEdge_iter = 0; hyEdge_iter < N_HyperEdges && std::getline(infile, line); ++hyEdge_iter )
  {
    linestream = std::istringstream(line);
    for (unsigned int i = 0; i < domain_info.hyNodes_hyEdge[hyEdge_iter].size(); ++i)
      linestream >> domain_info.hyNodes_hyEdge[hyEdge_iter][i];
  }
  
  hy_assert( hyEdge_iter == N_HyperEdges ,
             "Not all hyperedges have been added to the list!" );
  
  while ( keyword != "POINTS_OF_HYPEREDGES:" && std::getline(infile, line) )
  {
    linestream = std::istringstream(line);
    linestream >> keyword;
  }
  
  hy_assert( keyword == "POINTS_OF_HYPEREDGES:" ,
             "The keyword 'POINTS_OF_HYPEREDGES:' has not been found in the file "
             << filename << "!" );
  
  for ( hyEdge_iter = 0; hyEdge_iter < N_HyperEdges && std::getline(infile, line); ++hyEdge_iter )
  {
    linestream = std::istringstream(line);
    for (unsigned int i = 0; i < domain_info.points_hyEdge[hyEdge_iter].size(); ++i)
      linestream >> domain_info.points_hyEdge[hyEdge_iter][i];
  }
  
  hy_assert( hyEdge_iter == N_HyperEdges ,
             "Not all hyperedges have been added to the list!" );

  return domain_info;
} // end of read_domain_geo


template < unsigned int hyEdge_dim, unsigned int space_dim >
DomainInfo<hyEdge_dim,space_dim> read_domain( const std::string& filename )
{
  hy_assert( filename.substr(filename.size()-4, filename.size()) == ".geo" ,
             "The given file needs to be a .geo file, since no other input file types are currently"
             << " implemented." );
  
  DomainInfo<hyEdge_dim,space_dim> domain_info = read_domain_geo<hyEdge_dim,space_dim>(filename);
  
  hy_assert( domain_info.check_consistency() ,
             "Domain info appears to be inconsistent!" << std::endl
             << "This assertion is never to be thrown since it can only be caused by internal "
             << "assertions of DomainInfo.check_consistency()!" );
  
  return domain_info;
} // end of read_domain
