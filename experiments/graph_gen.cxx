#include <CLI/CLI.hpp>
#include <print>
#include <string>
#include <fstream>
#include "geobin.hxx"
#include <fmtlog/fmtlog.h>

namespace {

geobin::Graph generate_grid_graph(const char* path, geobin::u64 n) {
  geobin::Graph graph;
  int N = n*n;
  int M = 2 * (n-1) * n;
  double h = 1. / (n-1);
  geobin::Point p = {0, 0, 0};
  geobin::Prop cprop;
  std::fill(cprop.begin(), cprop.end(), 1);

  for (geobin::u64 i = 0; i < n; i++) {
    for (geobin::u64 j = 0; j < n; j++) {
      graph.vertices.push_back(p);
      if (i+1 < n)
        graph.edges.push_back({i*n+j, (i+1)*n+j});
      if (j+1 < n)
        graph.edges.push_back({i*n+j, i*n+j+1});

      p[1] += h;
    }
    p[1] = 0.;
    p[0] += h;
  }

  for (geobin::u64 e = 0; e < graph.edges.size(); e++)
    graph.edge_props.push_back(cprop);

  geobin::compute_types(graph);

  return graph;
}

}

int main(int argc, char** argv) {
  CLI::App app("synthetic graph generation");
  fmtlog::startPollingThread(1);

  std::string output_path = "graph";
  app.add_option("output", output_path, "the output path of the graph")->required();

  geobin::u64 n = 10;
  app.add_option("n", n, "related to number of vertices")->required();

  CLI11_PARSE(app, argc, argv);

  logi("args");
  logi("  output_path={}", output_path);
  logi("  n={}", n);

  geobin::Graph graph = generate_grid_graph(output_path.c_str(), n);
  geobin::serialize_bin(output_path.c_str(), graph);
}
