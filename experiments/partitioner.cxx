#include <print>
#include <fstream>
#include <limits>

#include "libpartition.hxx"
#include "geobin.hxx"
#include "stats.hxx"
#include <fmtlog/fmtlog.h>
#include <CLI/CLI.hpp>

namespace {

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

  geobin::ID partitions;
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
  app.add_option("-b,--backend", backend_str, "the partitioner backend to use, must be one of (kahip|metis|naive)");

  int kahip_mode = 2;
  app.add_option("--kahip-mode", kahip_mode, "the mode to run the kahip backend in, one of (0|1|2) representing FAST,ECO,STRONG, default=STRONG");
  bool kahip_no_suppress_output = false;
  app.add_option("--kahip-no-suppress-output", kahip_no_suppress_output, "do not suppress kahip backend output");
  int kahip_seed = 0;
  app.add_option("--kahip-seed", kahip_seed, "kahip seed, default=0");

  std::string test_overlap;
  app.add_option("--test-overlap", test_overlap, "test the overlap algorithm");

  geobin::ID partitions_z = 1;
  app.add_option("--partitions-z", partitions_z, "set the number of paritions in z direction when using backend 'naive'");

  double imbalance = 0.03;
  app.add_option("--imbalance", imbalance, "set the desired maximum imbalance in the algebraic partitioners");

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

  geobin::Graph graph = geobin::deserialize_bin(input_path.c_str());
  geobin::ID nverts = graph.vertices.size();
  geobin::ID nedges = graph.edges.size();

  logi("graph stats");
  logi("  nodes={}", nverts);
  logi("  edges={}", nedges);
  //
  // TODO: set default hops based on some graph stats like diameter, girth, etc
  //   maybe using 2-approximation or 3/2-approximation of diameter

  logi("partitioner backend");

  geobin::ID edgecut;
  std::vector<geobin::ID> partition(nverts);
  libpartition::PartConfig config = {libpartition::str_to_backend.at(frozen::string(backend_str)), !kahip_no_suppress_output, kahip_seed, kahip_mode, partitions_z};
  libpartition::do_partition(&graph, &partitions, &imbalance, partition.data(), &edgecut, &config);

  // TODO: compute balance
  logi("partition metrics");
  logi("  cut = {}", edgecut);
  logi("  bal = {}", -1);

  if (!vtu_output_path.empty()) {
    logi("serializing graph partition to vtu");
    serialize_graph_partition_vtu(graph, partition, vtu_output_path.c_str());
  }

  logi("make partition overlap and filter boundary nodes");

  std::vector<std::vector<geobin::ID>> domains(partitions);
  for (geobin::ID n = 0; n < nverts; n++)
    domains[partition[n]].push_back(n);

  std::vector<double> sizes_before(partitions);
  for (geobin::ID p = 0; p < partitions; p++)
    sizes_before[p] = domains[p].size();
  SimpleStats stats_before;
  compute_stats(sizes_before.data(), partitions, &stats_before);
  print_stats("sizes before", &stats_before);

  std::vector<geobin::u8> visited(graph.vertices.size());
  for (geobin::ID p = 0; p < partitions; p++) {
    std::fill(visited.begin(), visited.end(), 0); // slowest?
    for (const geobin::ID& n : domains[p])
      visited[n] = 1;

    // frontier marker
    domains[p].push_back((geobin::ID)-1);

    geobin::ID hop = 0;
    for (geobin::ID bfs_front = 0; bfs_front < domains[p].size() && hop < hops; bfs_front++) {
      const geobin::ID n = domains[p][bfs_front];

      // if we see a frontier marker, then hop is complete
      if (n == (geobin::ID)-1) {
        hop++;
        if (domains[p].back() == (geobin::ID)-1) {
          logi("domain {}: bfs terminated early after {} hops", p, hop, hops);
          logi("  no new nodes added in last hop");
          break;
        }
        domains[p].push_back((geobin::ID)-1);
        continue;
      }

      // all non visited (hence other partition) neighbors are added to the overlapping domain
      for (geobin::ID i = graph.xadj[n]; i < graph.xadj[n+1]; i++) {
        geobin::ID nn = graph.adjncy[i]; // neighbor
         if (!visited[nn]) {
          domains[p].push_back(nn);
          visited[nn] = 1;
        }
      }
    }
  }

  // remove frontier markers AND dirichlet nodes (type 1)
  for (geobin::ID p = 0; p < partitions; p++) {
    geobin::u64 offset = 0;
    for (geobin::u64 i = 0; i+offset < domains[p].size(); ) {
      geobin::ID node = domains[p][i+offset];
      if (node == (geobin::ID)-1 || graph.node_types[node] == 1) {
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
  for (geobin::ID p = 0; p < partitions; p++)
    sizes_after[p] = domains[p].size();
  SimpleStats stats_after;
  compute_stats(sizes_after.data(), partitions, &stats_after);
  print_stats("sizes after", &stats_after);

  std::vector<double> sizes_fractions(partitions);
  for (geobin::ID p = 0; p < partitions; p++)
    sizes_fractions[p] = sizes_after[p] / sizes_before[p];
  SimpleStats stats_frac;
  compute_stats(sizes_fractions.data(), partitions, &stats_frac);
  print_stats("after/before = ", &stats_frac);

  logi("writing overlapping partition");

  geobin::serialize_domains(output_path.c_str(), domains);

  if (!test_overlap.empty()) {
    std::vector<geobin::ID> fake_partition(graph.vertices.size(), 0);
    for (geobin::ID n : domains[0])
      fake_partition[n] = 1;
    serialize_graph_partition_vtu(graph, fake_partition, test_overlap.c_str());
  }
}
