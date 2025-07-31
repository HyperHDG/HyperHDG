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

  geobin::DataTable ioffset_table = {
    .name = "IOFFSET",
    .offset = sizeof(geobin::DomainsHeader),
    .size = (1+domains.size())*sizeof(NodeID),
  };

  geobin::DataTable domains_table = {
    .name = "DOMAINS",
    .offset = ioffset_table.offset + ioffset_table.size,
    .size = 0,
  };
  for (geobin::u64 p = 0; p < domains.size(); p++)
    domains_table.size += domains[p].size() * sizeof(NodeID);

  geobin::DomainsHeader hdr = {
    .magic = "DOMAIN1",
    .idsize = sizeof(NodeID),
    .n_domains = domains.size(),
    .tables = {
      ioffset_table,
      domains_table,
    }
  };
  file.write((char*)&hdr, sizeof(hdr));

  NodeID ioff = 0;
  std::vector<NodeID> ioffsets(1+domains.size());
  for (geobin::u64 p = 0; p < domains.size(); p++) {
    ioffsets[p] = ioff;
    ioff += domains[p].size();
  }
  ioffsets[domains.size()] = ioff;
  file.write((char*)&ioffsets[0], ioffsets.size()*sizeof(NodeID));

  for (PartitionID p = 0; p < domains.size(); p++)
    file.write((char*)&domains[p][0], domains[p].size()*sizeof(NodeID));
}

struct SimpleStats {
  double min, max, sum, avg, stddev;
};

void compute_stats(double* values, size_t n, SimpleStats* stats) {
  if (n == 0)
    return;
  stats->sum = stats->max = stats->avg = stats->stddev = 0;
  stats->min = values[0];
  for (geobin::u64 p = 0; p < n; p++) {
    stats->sum += values[p];
    stats->max = std::max(stats->max, values[p]);
    stats->min = std::min(stats->min, values[p]);
  }
  stats->avg = stats->sum / n;
  for (geobin::u64 p = 0; p < n; p++) {
    double d = values[p] - stats->avg;
    stats->stddev += d*d;
  }
  stats->stddev = std::sqrt(1./(n-1) * stats->stddev);
}

void print_stats(const char* msg, SimpleStats* stats) {
  logi("{}", msg);
  logi("  min={}", stats->min);
  logi("  max={}", stats->max);
  logi("  sum={}", stats->sum);
  logi("  avg={}", stats->avg);
  logi("  std={}", stats->stddev);
}

}

int main(int argc, char** argv) {
  fmtlog::startPollingThread(1);
  CLI::App app("partitioning tool");
  argv = app.ensure_utf8(argv);
  namespace fs = std::filesystem;

  // ARGUMENTS

  std::string input_path;
  app.add_option("input", input_path, "input path to the network to be partitioned")
    ->required();

  PartitionID partitions;
  app.add_option("partitions", partitions, "number of partitions to create")
    ->required();

  std::string output_path = std::format("{}.dom.zstd", input_path);
  app.add_option("output", output_path, "output path to the domain to be partitioned")
    ->required();

  // OPTIONS

  std::string vtu_output_path;
  app.add_option("--vtu", vtu_output_path, "serialize the graph to a vtu file including the partition");

  geobin::u64 hops = 3;
  app.add_option("--hops", hops, "number of hops to enlarge partitions by");

  std::string backend_str = "kahip";
  app.add_option("-b,--backend", backend_str, "the partitioner backend to use, must be one of (kahip|metis)");

  std::string test_overlap;
  app.add_option("--test-overlap", test_overlap, "test the overlap algorithm");

  CLI11_PARSE(app,argc,argv);

  logi("args");
  logi("  input_path={}", input_path);
  logi("  partitions={}", partitions);
  logi("  output_path={}", output_path);
  logi("  vtu_output_path={}", vtu_output_path);
  logi("  hops={}", hops);
  logi("  backend_str={}", backend_str);
  logi("  test_overlap={}", test_overlap);

  logi("reading graph");

  graph_access graph_acc;

  static_assert(std::is_same<NodeID, geobin::ID>());
  geobin::GraphEdgeList graph_edge_list = geobin::deserialize_bin(input_path.c_str());
  // map types in edges to types per node
  // NOTE: we assume that the type for each node is independent of the edge the node is in
  std::vector<NodeID> node_types(graph_edge_list.vertices.size(), (NodeID)-1);
  for (geobin::u64 e = 0; e < graph_edge_list.edges.size(); e++) {
    const geobin::Edge edge = graph_edge_list.edges[e];
    const geobin::Edge edge_types = graph_edge_list.types[e];
    node_types[edge.first] = edge_types.first;
    node_types[edge.second] = edge_types.second;
  }

  graph_edge_list.to_access(graph_acc);

  logi("graph stats");
  logi("  nodes={}", graph_acc.number_of_nodes());
  logi("  edges={}", graph_acc.number_of_edges()/2);
  //
  // TODO: set default hops based on some graph stats like diameter, girth, etc

  logi("partitioner backend");

  std::vector<PartitionID> partition(graph_edge_list.vertices.size());


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

  forall_nodes(graph_acc, n) {
    partition[n] = graph_acc.getPartitionIndex(n);
  } endfor

  if (!vtu_output_path.empty()) {
    logi("serializing graph partition to vtu");
    serialize_graph_partition_vtu(graph_edge_list, partition, vtu_output_path.c_str());
  }

  logi("make partition overlap and filter boundary nodes");

  std::vector<std::vector<NodeID>> domains(partitions);
  forall_nodes(graph_acc, n) {
    domains[partition[n]].push_back(n);
  } endfor

  std::vector<double> sizes_before(partitions);
  for (PartitionID p = 0; p < partitions; p++)
    sizes_before[p] = domains[p].size();
  SimpleStats stats_before;
  compute_stats(sizes_before.data(), partitions, &stats_before);
  print_stats("sizes before", &stats_before);

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
        if (domains[p].back() == (NodeID)-1) {
          logi("domain {}: bfs terminated early after {} hops", p, hop, hops);
          logi("  no new nodes added in last hop");
          break;
        }
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
  }

  // remove frontier markers AND dirichlet nodes (type 1)
  for (PartitionID p = 0; p < partitions; p++) {
    geobin::u64 offset = 0;
    for (geobin::u64 i = 0; i+offset < domains[p].size(); ) {
      NodeID node = domains[p][i+offset];
      if (node == (NodeID)-1 || node_types[node] == 1) {
        offset++;
      } else {
        domains[p][i] = node;
        i++;
      }
    }
    domains[p].resize(domains[p].size()-offset);
  }

  // TODO: integrate Q1 partitioning here to better compare
  // TODO: check if cython builds optimized...

  std::vector<double> sizes_after(partitions);
  for (PartitionID p = 0; p < partitions; p++)
    sizes_after[p] = domains[p].size();
  SimpleStats stats_after;
  compute_stats(sizes_after.data(), partitions, &stats_after);
  print_stats("sizes after", &stats_after);

  std::vector<double> sizes_fractions(partitions);
  for (PartitionID p = 0; p < partitions; p++)
    sizes_fractions[p] = sizes_after[p] / sizes_before[p];
  SimpleStats stats_frac;
  compute_stats(sizes_fractions.data(), partitions, &stats_frac);
  print_stats("after/before = ", &stats_frac);

  logi("writing overlapping partition");

  serialize_domains(output_path.c_str(), domains);

  if (!test_overlap.empty()) {
    std::vector<PartitionID> fake_partition(graph_edge_list.vertices.size(), 0);
    for (NodeID n : domains[0])
      fake_partition[n] = 1;
    serialize_graph_partition_vtu(graph_edge_list, fake_partition, test_overlap.c_str());
  }
}
