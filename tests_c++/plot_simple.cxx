#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/local_solver/template.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/plot.hxx>

#include <iostream>
#include <memory>

using namespace std;

/*!*************************************************************************************************
 * \brief   Output cubic geometries
 *
 * \authors   Guido Kanschat, Heidelberg University, 2020.
 * \authors   Andreas Rupp, Heidelberg University, 2020.
 **************************************************************************************************/
template <int edge_dim, int space_dim>
void test()
{
  vector<unsigned int> num_elements;

  if (space_dim >= 1)
    num_elements.push_back(2);
  if (space_dim >= 2)
    num_elements.push_back(3);
  if (space_dim >= 3)
    num_elements.push_back(4);

  typedef Topology::Cubic<edge_dim, space_dim> Topo;
  typedef Geometry::UnitCube<edge_dim, space_dim> Geo;
  typedef NodeDescriptor::Cubic<edge_dim, space_dim> Node;

  auto topo = std::make_shared<Topo>(num_elements);
  auto geo = std::make_shared<Geo>(*topo);
  auto node = std::make_shared<Node>(num_elements);

  typedef LocalSolver::Template<edge_dim, double> SolverType;
  SolverType lsolver;
  HDGHyperGraph<SolverType::n_glob_dofs_per_node(), Topo, Geo, Node, typename SolverType::data_type>
    hdg_graph(topo, geo, node);  // Must be according to the local solver!
  vector<double> vectorDirichlet(hdg_graph.n_global_dofs());

  PlotOptions pltop;
  pltop.scale = 0.8;
  pltop.plot_edge_boundaries = true;
  pltop.plot_edges = true;
  pltop.boundary_scale = 0.9;
  std::string name = "plot1-";
  name += std::to_string(edge_dim) + std::string("-") + char('a' + space_dim - edge_dim);
  pltop.fileName = name;
  pltop.printFileNumber = false;
  plot(hdg_graph, lsolver, vectorDirichlet, pltop);
}

int main()
{
  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<1, 3>();
  test<2, 3>();
  test<3, 3>();
  return 0;
}
