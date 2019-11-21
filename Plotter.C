/* ------------------------------------------------------------------------------------------------------
 *
 * This file is part of EP2 of the STRUCTURES initiative of the University of Heidelberg.
 * It solves a PDE that is solely defined on a graph using the HDG method.
 * 
 * Definition of a class for point objects providing getter functions only.
 *
 * ------------------------------------------------------------------------------------------------------
 *
 * Author: Andreas Rupp, University of Heidelberg, 2019
 */


#include "Plotter.h"
#include "TypeDefs.h"
#include <fstream>
#include <iomanip>
#include <cmath>


using namespace std;


template class Plotter <1, 1, 1>;
template class Plotter <1, 2, 1>;
template class Plotter <1, 3, 1>;
template class Plotter <2, 2, 1>;
template class Plotter <2, 3, 1>;
template class Plotter <3, 3, 1>;


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
Plotter<hyperedge_dim,space_dim,polynomial_degree>::
Plotter(HyperGraph< local_dof_amount_node(hyperedge_dim, polynomial_degree), HyperGraph_Cubic< hyperedge_dim, space_dim >, HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >& hyper_graph,
        DiffusionSolver_RegularQuad<hyperedge_dim, polynomial_degree, 2 * polynomial_degree>& local_solver)
: hyper_graph_(hyper_graph), local_solver_(local_solver) { }


template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
void Plotter<hyperedge_dim,space_dim,polynomial_degree>::
plot(vector<double> lambda)
{
  hyperedge_index_type num_of_hyperedges = hyper_graph_.num_of_hyperedges();
  unsigned int points_per_hyperedge = pow(2, hyperedge_dim);
  point_index_type num_of_points = points_per_hyperedge * num_of_hyperedges;
//  HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim> hyperedge_geometry;
  
  static_assert (hyperedge_dim <= 3);
  unsigned int element_id;
  if constexpr (hyperedge_dim == 1)       element_id = 3;
  else if constexpr (hyperedge_dim == 2)  element_id = 8;
  else if constexpr (hyperedge_dim == 3)  element_id = 11;
  
  Point<space_dim> point;
  ofstream myfile;
  
  
	myfile.open("output/example.vtu");
	myfile << "<?xml version=\"1.0\"?>"  << endl;
	myfile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << endl;
	myfile << "  <UnstructuredGrid>" << endl;
	myfile << "    <Piece NumberOfPoints=\"" << num_of_points << "\" NumberOfCells= \"" << num_of_hyperedges << "\">" << endl;
	myfile << "      <Points>" << endl;
	myfile << "        <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;
  
  for (hyperedge_index_type he_number = 0; he_number < num_of_hyperedges; ++he_number)
  {
    HyperEdge_Cubic_UnitCube<hyperedge_dim, space_dim> hyperedge_geometry = hyper_graph_.get_hyperedge_geometry(he_number);
    for (unsigned int pt_number = 0; pt_number < points_per_hyperedge; ++pt_number)
    {
      myfile << "        ";
      point = hyperedge_geometry.point(pt_number);
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        myfile << "  " << fixed << scientific << setprecision(3) << point[dim];
      for (unsigned int dim = space_dim; dim < 3; ++dim)
        myfile << "  0.0";
      myfile << endl;
    }
  }
  
  myfile << "        </DataArray>" << endl;
  myfile << "      </Points>" << endl;
	myfile << "      <Cells>" << endl;
	myfile << "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">" << endl;
  myfile << "        ";
  
	for (point_index_type pt_number = 0; pt_number < num_of_points; ++pt_number)
    myfile << "  " << pt_number;
  myfile << endl;
  
  myfile << "        </DataArray>" << endl;
  myfile << "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">" << endl;
  myfile << "        ";
  
  for (point_index_type pt_number = points_per_hyperedge; pt_number <= num_of_points; pt_number += points_per_hyperedge)
    myfile << "  " << pt_number;
  myfile << endl;
  
  myfile << "        </DataArray>" << endl;
	myfile << "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
  myfile << "        ";
  
  for (hyperedge_index_type he_number = 0; he_number < num_of_hyperedges; ++he_number)
    myfile << "  " << element_id;
  myfile << endl;
  
  myfile << "        </DataArray>" << endl;
	myfile << "      </Cells>" << endl;
  
  
  myfile << "      <PointData Scalars=\"" << "example_scalar" << "\" Vectors=\"" << "example_vector" << "\">" << endl;
  myfile << "        <DataArray type=\"Float32\" Name=\"" << "dual" << "\" NumberOfComponents=\"" << hyperedge_dim << "\" format=\"ascii\">" << endl;
    
  array< array<double, local_dof_amount_node(hyperedge_dim, polynomial_degree)> , 2*hyperedge_dim > hyperedge_dofs;
  array<unsigned int, 2*hyperedge_dim> hyperedge_hypernodes;
  array<double, corners_amount(hyperedge_dim)> local_primal;
  array< array<double, hyperedge_dim> , corners_amount(hyperedge_dim) > local_dual;
  
  for (hyperedge_index_type he_number = 0; he_number < num_of_hyperedges; ++he_number)
  {
    hyperedge_hypernodes = hyper_graph_.get_hyperedge(he_number).get_hypernode_indices();
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyperedge_dofs[hypernode] = hyper_graph_.hypernode_factory().get_dof_values(hyperedge_hypernodes[hypernode], lambda);
    local_dual = local_solver_.dual_in_corners_from_lambda(hyperedge_dofs);
    myfile << "      ";
    for (unsigned int corner = 0; corner < corners_amount(hyperedge_dim); ++corner)
    {
      myfile << "  ";
      for (unsigned int dim = 0; dim < hyperedge_dim; ++dim)
        myfile << "  " << local_dual[corner][dim];
    }
    myfile << endl;
  }

  myfile << "        </DataArray>" << endl;
  myfile << "        <DataArray type=\"Float32\" Name=\"" << "primal" << "\" NumberOfComponents=\"1\" format=\"ascii\">" << endl;
  
  
  for (hyperedge_index_type he_number = 0; he_number < num_of_hyperedges; ++he_number)
  {
    hyperedge_hypernodes = hyper_graph_.get_hyperedge(he_number).get_hypernode_indices();
    for (unsigned int hypernode = 0; hypernode < hyperedge_hypernodes.size(); ++hypernode)
      hyperedge_dofs[hypernode] = hyper_graph_.hypernode_factory().get_dof_values(hyperedge_hypernodes[hypernode], lambda);
    local_primal = local_solver_.primal_in_corners_from_lambda(hyperedge_dofs);
    myfile << "        ";
    for (unsigned int corner = 0; corner < corners_amount(hyperedge_dim); ++corner)
      myfile << "  " << local_primal[corner];
    myfile << endl;
  }

  myfile << "        </DataArray>" << endl;
  myfile << "      </PointData>" << endl;
  

	myfile << "    </Piece>" << endl;
	myfile << "  </UnstructuredGrid>" << endl;
	myfile << "</VTKFile>" << endl;
  myfile.close();
}
