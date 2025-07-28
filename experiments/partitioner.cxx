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
#include <CLI/CLI.hpp>
#include <bxzstr.hpp>

namespace {

void serialize_graph_partition_vtu(
  const geobin::GraphEdgeList& graph,
  const std::vector<PartitionID>& partition,
  const char* file_path
) {
  std::ofstream file(std::format("{}.vtu", file_path));
  if (!file.is_open()) {
    loge("could not opt file {}.vtu", file_path);
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

void serialize_domains(const char* path, const std::vector<std::vector<NodeID>>& domains) {
  bxz::ofstream file(std::format("{}.dom.zstd", path), bxz::zstd);

  geobin::DomainsHeader hdr = {
    .magic = "DOMAIN1",
    .idsize = sizeof(NodeID),
    .n_domains = domains.size(),
  };
  file.write((char*)&hdr, sizeof(hdr));

  std::vector<geobin::DataTable> domain_tables(domains.size());
  geobin::u64 offset = sizeof(hdr) + domain_tables.size() * sizeof(geobin::DataTable);
  for (PartitionID p = 0; p < domains.size(); p++) {
    strncpy(domain_tables[p].name, "DOMAIN", sizeof(domain_tables[p].name));
    domain_tables[p].offset = offset;
    domain_tables[p].size = domains[p].size() * sizeof(NodeID);
    offset += domain_tables[p].size;
    file.write((char*)&domain_tables[p], sizeof(geobin::DataTable));
  }

  for (geobin::u64 p = 0; p < domains.size(); p++)
    file.write((char*)domains[p].data(), domain_tables[p].size);
}

}

int main(int argc, char** argv) {
  fmtlog::startPollingThread(1);
  CLI::App app("partitioning tool");
  argv = app.ensure_utf8(argv);
  namespace fs = std::filesystem;

  // ARGUMENTS

  std::string input_path;
  app.add_option("input", input_path, "input path to the network to be partitioned");

  PartitionID partitions;
  app.add_option("partitions", partitions, "number of partitions to create");

  std::string output_path = std::format("{}.dom.zstd", input_path);
  app.add_option("output", output_path, "output path to the domain to be partitioned");

  // OPTIONS

  std::string vtu_output_path;
  app.add_option("--vtu", vtu_output_path, "serialize the graph to a vtu file including the partition");

  geobin::u64 hops = 3;
  app.add_option("--hops", hops, "number of hops to enlarge partitions by");

  std::string backend_str = "kahip";
  app.add_option("-b,--backend", backend_str, "the partitioner backend to use, must be one of (kahip|metis)");

  std::string test_overlap;
  app.add_option("--test-overlap", test_overlap, "test the overlap algorithm");

  logi("args");
  logi("  input_path={}");
  logi("  output_path={}");
  logi("  vtu_output_path={}");
  logi("  hops={}");
  logi("  backend_str={}");
  logi("  test_overlap={}");

  CLI11_PARSE(app,argc,argv);

  logi("reading graph");

  graph_access graph_acc;

  static_assert(std::is_same<NodeID, geobin::ID>());
  geobin::GraphEdgeList graph_edge_list = geobin::deserialize_bin(input_path.c_str());
  graph_edge_list.to_access(graph_acc);

  logi("partitioner backend");

  // TODO: support metis

  graph_partitioner partitioner;
  PartitionConfig partition_config;
  configuration cfg;
  cfg.strong(partition_config);
  partition_config.k = partitions;
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

  if (!vtu_output_path.empty()) {
    logi("serializing graph partition to vtu");
    serialize_graph_partition_vtu(graph_edge_list, partition, vtu_output_path.c_str());
  }

  logi("make partition overlap");

  std::vector<std::vector<NodeID>> domains(partitions);
  forall_nodes(graph_acc, n) {
    domains[partition[n]].push_back(n);
  } endfor

  std::vector<geobin::u8> visited(graph_edge_list.vertices.size());
  for (PartitionID p = 0; p < partitions; p++) {
    std::fill(visited.begin(), visited.end(), 0); // slowest?
    for (const NodeID& n : domains[p])
      visited[n] = 1;

    // frontier marker
    domains[p].push_back((NodeID)-1);

    geobin::u64 hop = 0;
    for (geobin::u64 bfs_front = 0; bfs_front < domains[p].size() && hop < hops; bfs_front++) {
      const NodeID n = domains[p][bfs_front];

      // if we see a frontier marker, then hop is complete
      if (n == (NodeID)-1) {
        hop++;
        domains[p].push_back((NodeID)-1);
        continue;
      }

      // all non visited (hence other partition) neighbors are added to the overlapping domain
      forall_out_edges(graph_acc, e, n) {
        NodeID nn = graph_acc.getEdgeTarget(e);
        if (!visited[nn]) {
          domains[p].push_back(nn);
          visited[nn] = 1;
        }
      } endfor
    }

    if (hop < hops)
      logi("domain {}: bfs terminated early after {} hops", p, hop, hops);
  }

  // remove frontier markers
  for (PartitionID p = 0; p < partitions; p++) {
    geobin::u64 offset = 0;
    for (geobin::u64 i = 0; i+offset < domains[p].size(); i++) {
      if (domains[p][i+offset] == (NodeID)-1)
        offset++;
      if (i+offset < domains[p].size())
        domains[p][i] = domains[p][i+offset];
    }
    domains[p].resize(domains[p].size()-offset);
  }

  logi("writing overlapping partition");

  serialize_domains(output_path.c_str(), domains);

  if (!test_overlap.empty()) {
    std::vector<PartitionID> fake_partition(graph_edge_list.vertices.size(), 0);
    for (NodeID n : domains[0])
      fake_partition[n] = 1;
    serialize_graph_partition_vtu(graph_edge_list, fake_partition, test_overlap.c_str());
  }
}
