#include <cstdio>
#include <string>
#include "geobin.hxx"

geobin::Graph generate_grid_graph(const char* path, geobin::u64 n) {
  geobin::Graph graph;
  double h = 1. / (n-1);
  geobin::Point p = {0, 0, 0};
  geobin::Prop cprop;
  std::fill(cprop.begin(), cprop.end(), 1);

  for (geobin::u64 i = 0; i < n; i++) {
    for (geobin::u64 j = 0; j < n; j++) {
      p[0] = j*h;
      p[1] = i*h;
      graph.vertices.push_back(p);
      if (i+1 < n)
        graph.edges.push_back({i*n+j, (i+1)*n+j});
      if (j+1 < n)
        graph.edges.push_back({i*n+j, i*n+j+1});
    }
  }

  for (geobin::u64 e = 0; e < graph.edges.size(); e++)
    graph.edge_props.push_back(cprop);

  geobin::compute_types(graph);

  return graph;
}

int usage(int argc, char** argv) {
  fprintf(stderr, "ERROR: usage: %s <n> <path>\n", argv[0]);
  return 1;
}

int main(int argc, char** argv) {
  const char* out;
  geobin::u64 n;

  if (argc < 3) return usage(argc, argv);
  if (1 != sscanf(argv[1], "%lu", &n)) return usage(argc, argv);
  out = argv[2];

  printf("args\n");
  printf("  out=%s\n", out);
  printf("  n=%lu\n", n);

  geobin::Graph graph = generate_grid_graph(out, n);
  geobin::serialize_bin(out, graph);
}
