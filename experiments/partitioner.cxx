#include <print>
#include <fstream>

#include "geobin.hxx"
#include <KaHIP/app/configuration.h>
#include <KaHIP/app/balance_configuration.h>
#include <KaHIP/lib/partition/graph_partitioner.h>
#include <KaHIP/lib/data_structure/graph_access.h>
#include <KaHIP/lib/tools/random_functions.h>
#include <KaHIP/lib/tools/quality_metrics.h>
#include <fmtlog/fmtlog.h>

namespace {

void serialize_graph_partition_vtu(
  const geobin::GraphEdgeList& graph,
  const std::vector<PartitionID>& partition,
  const char* file_path
) {
  std::ofstream file(file_path);
  if (!file.is_open()) {
    loge("could not opt file {}", file_path);
    return;
  }

  std::print(file, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  std::print(file, "  <UnstructuredGrid>\n");
  std::print(file, "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n", graph.vertices.size(), graph.edges.size());
  std::print(file, "      <Points>\n");
  std::print(file, "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (const auto& point : graph.vertices)
      std::print(file, "          {} {} {}\n", point[0], point[1], point[2]);
  std::print(file, "        </DataArray>\n");
  std::print(file, "      </Points>\n");
  std::print(file, "      <Cells>\n");
  std::print(file, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  std::print(file, "          ");
  for (const auto& edge : graph.edges)
      std::print(file, "{} {} ", edge.first, edge.second);
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  // offsets: the cumulative sum of the number of points in each cell.
  // for VTK_LINE (2 points per cell), this will be 2, 4, 6, ...
  std::print(file, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  std::print(file, "          ");
  for (size_t i = 0; i < graph.edges.size(); ++i)
      std::print(file, "{} ", (i + 1) * 2);
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  std::print(file, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  std::print(file, "          ");
  for (size_t i = 0; i < graph.edges.size(); ++i)
      std::print(file, "3 "); // VTK_LINE type
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  std::print(file, "      </Cells>\n");
  std::print(file, "      <PointData>\n");
  std::print(file, "        <DataArray type=\"Int32\" Name=\"PartitionID\" format=\"ascii\">\n");
  std::print(file, "          ");
  for (PartitionID id : partition)
      std::print(file, "{} ", id);
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  std::print(file, "      </PointData>\n");
  std::print(file, "    </Piece>\n");
  std::print(file, "  </UnstructuredGrid>\n");
  std::print(file, "</VTKFile>\n");
}

}

int main(int argc, char** argv) {
  fmtlog::startPollingThread(1);

  if (argc < 4) {
    std::println(stderr, "ERROR: usage: {} <input> <num_partitions> <output>", argv[0]);
    return 1;
  }

  logi("reading graph");

  geobin::GraphEdgeList graph_edge_list = geobin::deserialize_bin(argv[1]);
  PartitionID num_partitions = std::stol(argv[2]);
  const char* output_path = argv[3];

  // TODO: test if output_path is writable

  logi("constructing KaHIP graph_acc");

  // is directed -> 2*edges
  graph_access graph_acc;
  graph_acc.start_construction(graph_edge_list.vertices.size(), 2*graph_edge_list.edges.size());

  std::vector<std::vector<geobin::ID>> graph_adjacency(graph_edge_list.vertices.size());
  for (const geobin::Edge& edge : graph_edge_list.edges) {
    graph_adjacency[edge.first].push_back(edge.second);
    graph_adjacency[edge.second].push_back(edge.first);
  }

  for (NodeID n = 0; n < graph_edge_list.vertices.size(); n++) {
    NodeID nn = graph_acc.new_node();
    graph_acc.setNodeWeight(nn, 1);
    // graph_acc.setPartitionIndex(nn, 0);
    for (const geobin::ID neighbor : graph_adjacency[n]) {
      EdgeID e = graph_acc.new_edge(nn, neighbor);
      graph_acc.setEdgeWeight(e, 1);
    }
  }

  graph_acc.finish_construction();

  logi("KaHIP partitioner");

  graph_partitioner partitioner;
  PartitionConfig partition_config;
  configuration cfg;
  cfg.strong(partition_config);
  partition_config.k = num_partitions;
  partition_config.seed = 0;
  srand(partition_config.seed);
  random_functions::setSeed(partition_config.seed);
  graph_acc.set_partition_count(partition_config.k);
  balance_configuration bc;
  bc.configurate_balance(partition_config, graph_acc);

  partitioner.perform_partitioning(partition_config, graph_acc);

  quality_metrics qm;
  logi("KaHIP metrics");
  logi("  cut = {}", qm.edge_cut(graph_acc));
  logi("  bnd = {}", qm.boundary_nodes(graph_acc));
  logi("  bal = {}", qm.balance(graph_acc));

  std::vector<PartitionID> partition(graph_edge_list.vertices.size());
  forall_nodes(graph_acc, n) {
    partition[n] = graph_acc.getPartitionIndex(n);
  } endfor

  logi("serializing graph partitio to vtu");

  serialize_graph_partition_vtu(graph_edge_list, partition, output_path);
}
