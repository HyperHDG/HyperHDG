#include <print>
#include <fstream>
#include <limits>

#include "geobin.hxx"
#include "stats.hxx"
#include <kaHIP_interface.h>
#include <fmtlog/fmtlog.h>
#include <CLI/CLI.hpp>
#include <bxzstr.hpp>
#include <metis.h>

namespace {

void serialize_graph_partition_vtu(
  const geobin::Graph& graph,
  const std::vector<geobin::ID>& partition,
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
  for (geobin::ID id : partition)
      std::print(file, "{} ", id);
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  std::print(file, "      </PointData>\n");
  std::print(file, "    </Piece>\n");
  std::print(file, "  </UnstructuredGrid>\n");
  std::print(file, "</VTKFile>\n");
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

  if (backend_str == "kahip") {
    static_assert(sizeof(int) == sizeof(geobin::ID));
    kaffpa((int*)&nverts, NULL, (int*)graph.xadj.data(), NULL, (int*)graph.adjncy.data(), (int*)&partitions, &imbalance, !kahip_no_suppress_output, kahip_seed, kahip_mode, (int*)&edgecut, (int*)partition.data());
  } else if (backend_str == "metis") {
    static_assert(sizeof(idx_t) == sizeof(geobin::ID));
    METIS_PartGraphKway((idx_t*)&nverts, (idx_t*)&nedges, (idx_t*)graph.xadj.data(), (idx_t*)graph.adjncy.data(), NULL, NULL, NULL, (idx_t*)&partitions, NULL, NULL, NULL, (idx_t*)&edgecut, (idx_t*)partition.data());
  } else if (backend_str == "parhip") {
    loge("backend=parhip unsupported as of yet");
  }else if (backend_str == "naive") {
    // TODO: move this into header
    geobin::Point min_p = {std::numeric_limits<geobin::Real>::max()}, max_p = {std::numeric_limits<geobin::Real>::min()};
    for (geobin::u64 n = 0; n < graph.vertices.size(); n++) {
      const geobin::Point& p = graph.vertices[n];
      for (geobin::u64 i = 0; i < 3; i++) {
        min_p[i] = std::min(min_p[i], p[i]);
        max_p[i] = std::max(max_p[i], p[i]);
      }
    }

    geobin::Real eps = 1e-10;
    std::array<geobin::ID, 3> partitions3d = {(geobin::ID)std::sqrt(partitions/partitions_z), (geobin::ID)std::sqrt(partitions/partitions_z), partitions_z};
    for (geobin::u64 n = 0; n < graph.vertices.size(); n++) {
      const geobin::Point& p = graph.vertices[n];
      geobin::ID pid = 0;
      for (geobin::u64 i = 0; i < 3; i++) {
        pid *= partitions3d[i];
        pid += p[i]/((1+eps)*(max_p[i]-min_p[i])) * partitions3d[i]; // truncate
      }
      assert(pid < partitions);
      partition[n] = pid;
    }

    edgecut = 0;
    for (const geobin::Edge& edge : graph.edges) {
      if (partition[edge.first] != partition[edge.second])
        edgecut++;
    }
  } else {
    loge("ERROR: unsupported backend '{}'", backend_str);
    return 1;
  }

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
