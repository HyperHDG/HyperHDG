#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/local_solver/template.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/plot.hxx>

#include <iostream>
#include <memory>

using namespace std;

template <int edge_dim, int space_dim, int nsub>
void output_cubic()
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

  typedef LocalSolver::Template<edge_dim, double> Solver;
  Solver lsolver;
  typedef HDGHyperGraph<Solver::n_glob_dofs_per_node(), Topo, Geo, Node, typename Solver::data_type>
    HDG;
  HDG hdg_graph(topo, geo, node);
  vector<double> vectorDirichlet(hdg_graph.n_global_dofs());

  PlotOptions pltop;
  pltop.scale = .8;
  pltop.plot_edge_boundaries = true;
  pltop.numbers = true;
  std::string name = "plot2-";
  name += std::to_string(edge_dim) + std::string("-") + char('a' + space_dim - edge_dim);
  pltop.fileName = name;
  pltop.printFileNumber = false;
  plot_vtu<HDG, Solver, vector<double>, double, nsub>(hdg_graph, lsolver, vectorDirichlet, pltop);
}

int main()
{
  output_cubic<1, 1, 2>();
  output_cubic<1, 2, 2>();
  output_cubic<1, 3, 2>();
  output_cubic<2, 2, 2>();
  output_cubic<2, 3, 3>();
  output_cubic<3, 3, 2>();
  return 0;
}
