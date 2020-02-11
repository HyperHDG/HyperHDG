/*!*************************************************************************************************
 * \file    examples_c++/DiffusionTest.C
 * \brief   File that tests several aspects of the C++ implementation against a given reference
 *          solution obtained with the Python interface.
 *
 * This file implements an alternative to Executable.py (which usses the Cython interface).
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/

#include <Topo_Cubic.hxx>
#include <Geom_UnitCube.hxx>
#include <HDGHyperGraph.hxx>
#include <LSol_Diffusion.hxx>
#include <Plot.hxx>
#include <iostream>

#include <memory>

using namespace std;

/*!*************************************************************************************************
 * \brief   Output cubic geometries
 *
 * \authors   Guido Kanschat, University of Heidelberg, 2020.
 * \authors   Andreas Rupp, University of Heidelberg, 2020.
 **************************************************************************************************/
template<int edge_dim, int space_dim>
void test()
{
  vector<unsigned int> num_elements;

  if (space_dim >= 1) num_elements.push_back(4);
  if (space_dim >= 2) num_elements.push_back(3);
  if (space_dim >= 3) num_elements.push_back(2);

  typedef Topology::Cubic<edge_dim,space_dim> Topo;
  typedef Geometry::UnitCube<edge_dim,space_dim> Geo;
  auto topo = std::make_shared<const Topo> (num_elements);
  auto geo = std::make_shared<const Geo> (*topo);

  Diffusion_TensorialUniform<edge_dim,1,2> lsolver(1.); // Quadrature must be sufficient!
  HDGHyperGraph<Diffusion_TensorialUniform<edge_dim,1,1>::n_glob_dofs_per_node(),Topo,Geo>
    hdg_graph(topo, geo); // Must be according to the local solver!
  vector<double> vectorDirichlet(hdg_graph.n_global_dofs());
  
  PlotOptions pltop;
  std::string name = "plot1-";
    name += std::to_string(edge_dim) + std::string("-") + std::to_string(space_dim);
  pltop.fileName = name;
  pltop.printFileNumber = false;
  plot(hdg_graph, lsolver, vectorDirichlet, pltop);
}

int main()
{
  test<1,1>();
  test<1,2>();
  test<2,2>();
  test<1,3>();
  test<2,3>();
  test<3,3>();
  return 0;
}
