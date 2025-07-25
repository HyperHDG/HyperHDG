#include <print>

#include "geobin.hxx"
#include <KaHIP/lib/data_structure/graph_access.h>
#include <fmtlog/fmtlog.h>

int main(int argc, char** argv) {
  if (argc < 2) {
    std::println(stderr, "ERROR: usage: {} <input_network>", argv[0]);
    return 1;
  }

  geobin::GraphEdgeList graph_edge_list = geobin::deserialize_bin(argv[1]);
  logi("nodes = {}", graph_edge_list.vertices.size());
  logi("edges = {}", graph_edge_list.vertices.size());
  logi("done");
}
