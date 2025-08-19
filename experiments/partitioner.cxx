#include <print>
#include <fstream>
#include <limits>
#include <deque>
#include <random>

#include "libpartition.hxx"
#include "geobin.hxx"
#include "stats.hxx"
#include <CLI/CLI.hpp>
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/stopwatch.h>
#include <spdlog/cfg/helpers.h>

namespace {

void print_stats(const char* msg, SimpleStats* stats, const nlohmann::json& data) {
  auto lg = spdlog::get("logger");
  lg->info(msg)({
      {"min", stats->min}, {"max", stats->max}, {"sum", stats->sum}, {"avg", stats->avg}, {"std", stats->stddev}, data
  });
}

using UUID = std::array<geobin::u8,16>;
UUID generate_uuid() {
  UUID uuid = {0};
  uint32_t* u = (uint32_t*) uuid.data();
  std::random_device rd;
  std::mt19937 gen(rd());

  // A 128-bit UUID is 16 bytes. We'll generate it as four 32-bit integers.
  std::uniform_int_distribution<uint32_t> dis(0, 0xFFFFFFFF);
  for (size_t i = 0; i < 4; i++)
    u[i] = dis(gen);
  return uuid;
}

std::string uuid_to_str(const UUID& uuid) {
  char buf[32+1];
  uint32_t* u = (uint32_t*) uuid.data();
  for (size_t i = 0; i < 32; i += 8) {
    snprintf(buf+i, 32+1, "%08X", *u);
    u++;
  }
  return std::string(buf);
}

geobin::ID diameter_2approx(geobin::Graph& graph) {
  geobin::ID nverts = graph.vertices.size();
  std::mt19937 gen(std::random_device{}());
  std::uniform_int_distribution<> dist(0, nverts-1);
  geobin::ID arbitrary_start = dist(gen);

  std::vector<geobin::u8> visited(nverts, 0);
  std::deque<geobin::ID> deque {arbitrary_start};
  geobin::ID farthest = arbitrary_start;
  while (!deque.empty()) {
    geobin::ID n = deque.front();
    deque.pop_front();
    for (geobin::ID i = graph.xadj[n]; i < graph.xadj[n+1]; i++) {
      farthest = graph.adjncy[i];
      if (visited[farthest] != 1) {
        deque.push_back(farthest);
        visited[farthest] = 1;
      }
    }
  }

  geobin::ID distance = 0;
  deque.push_back(farthest);
  deque.push_back((geobin::ID)-1);
  visited[farthest] = 2;
  while (!deque.empty()) {
    geobin::ID n = deque.front();
    deque.pop_front();
    if (n == (geobin::ID)-1) {
      if (deque.empty())
        break;
      distance++;
      deque.push_back((geobin::ID)-1);
      continue;
    }
    for (geobin::ID i = graph.xadj[n]; i < graph.xadj[n+1]; i++) {
      geobin::ID nn = graph.adjncy[i];
      if (visited[nn] != 2) {
        deque.push_back(nn);
        visited[nn] = 2;
      }
    }
  }

  return distance;
}

}

int main(int argc, char** argv) {
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

  geobin::u64 delta = 2;
  app.add_option("--delta", delta, "overlap parameter delta = number of hops to enlarge partitions by");

  std::string backend = "kahip";
  app.add_option("-b,--backend", backend, "the partitioner backend to use, must be one of (kahip|metis|naive)");

  int kahip_mode = 2;
  app.add_option("--kahip-mode", kahip_mode, "the mode to run the kahip backend in, one of (0|1|2) representing FAST,ECO,STRONG, default=STRONG");
  bool kahip_suppress_output = true;
  app.add_option("--kahip-suppress-output", kahip_suppress_output, "suppress kahip backend output, default=true");
  int kahip_seed = 0;
  app.add_option("--kahip-seed", kahip_seed, "kahip seed, default=0");

  std::string test_overlap;
  app.add_option("--test-overlap", test_overlap, "test the overlap algorithm");

  geobin::ID partitions_z = 1;
  app.add_option("--partitions-z", partitions_z, "set the number of paritions in z direction when using backend 'naive'");

  double imbalance = 0.03;
  app.add_option("--imbalance", imbalance, "set the desired maximum imbalance in the algebraic partitioners");

  std::string log_file;
  app.add_option("--log-file", log_file, "log to file path instead of to stdout");

  std::string log_level = "info";
  app.add_option("--log-level", log_level, "set log level, one of  (trace|debug|info|warn|err)");

  bool extended_stats = false;
  app.add_option("--extended-stats", extended_stats, "compute extended (expensive) stats, default=false");

  bool square = true;
  app.add_option("--square", square, "square the number of requested partitions, default=true");

  std::string srunid;
  app.add_option("--runid", srunid, "set runid, default=random");

  CLI11_PARSE(app,argc,argv);

  auto lg = log_file.empty() ? spdlog::stdout_logger_st("logger") : spdlog::basic_logger_st("logger", log_file);
  spdlog::stopwatch sw;
  spdlog::cfg::helpers::load_levels(log_level);

  if (square)
    partitions = partitions*partitions;

  if (srunid.empty()) {
    UUID uuid = generate_uuid();
    srunid = uuid_to_str(uuid);
  }
  const nlohmann::json runid = {"runid", srunid};

  lg->info("args")({
      {"input_path", input_path},
      {"partitions", partitions},
      {"output_path", output_path},
      {"vtu_output_path", vtu_output_path},
      {"delta", delta},
      {"backend", backend},
      {"test_overlap", test_overlap},
      {"kahip_suppress_output", kahip_suppress_output},
      {"kahip_mode", kahip_mode},
      {"kahip_seed", kahip_seed},
      {"log_level", log_level},
      {"log_file", log_file},
      {"imbalance", imbalance},
      {"partitions_z", partitions_z},
      {"square", square},
      runid,
  });

  lg->info("reading graph...");

  sw.reset();
  geobin::Graph graph = geobin::deserialize_bin(input_path.c_str());
  geobin::ID nverts = graph.vertices.size();
  geobin::ID nedges = graph.edges.size();

  lg->debug("reading graph")({{"time", sw.elapsed().count()}, runid});
  lg->debug("graph stats")({{"nodes", nverts}, {"edges", nedges}, runid});

  if (extended_stats) {
    std::vector<double> degrees(nverts, 0);
    for (geobin::ID n = 0; n < nverts; n++)
      degrees[n] = (double)(graph.xadj[n+1] - graph.xadj[n]);
    SimpleStats stats_degrees;
    compute_stats(degrees.data(), nverts, &stats_degrees);
    print_stats("degree stats", &stats_degrees, runid);

    lg->info("diameter_2approx")({{"diameter_2approx", diameter_2approx(graph)}, runid});
  }

  // TODO: set default hops based on some graph stats like diameter, girth, etc
  //   maybe using 2-approximation or 3/2-approximation of diameter
  //   NOTE2: or maybe not? maybe 2 is fine?

  lg->info("partitioner backend...");

  sw.reset();
  geobin::ID edgecut;
  std::vector<geobin::ID> partition(nverts);
  libpartition::PartConfig config = {backend.c_str(), kahip_suppress_output, kahip_seed, kahip_mode, partitions_z};
  libpartition::do_partition(&graph, &partitions, &imbalance, partition.data(), &edgecut, &config);

  lg->debug("partitioner_backend")({{"time", sw.elapsed().count()}, runid});

  std::vector<std::vector<geobin::ID>> domains(partitions);
  for (geobin::ID n = 0; n < nverts; n++)
    domains[partition[n]].push_back(n);

  std::vector<double> sizes_before(partitions);
  for (geobin::ID p = 0; p < partitions; p++) {
    sizes_before[p] = domains[p].size();
    lg->trace("partition size")({{"id",p},{"size", sizes_before[p]}});
  }
  SimpleStats stats_before;
  compute_stats(sizes_before.data(), partitions, &stats_before);

  lg->debug("partition metrics")({runid, {"cut", edgecut}, {"bal", stats_before.max / ((double)nverts/partitions)}});

  if (!vtu_output_path.empty()) {
    lg->info("serializing graph partition to vtu");
    serialize_graph_partition_vtu(graph, partition, vtu_output_path.c_str());
  }

  lg->info("make partition overlap and filter boundary nodes...");

  print_stats("sizes before", &stats_before, runid);
  lg->debug("bal before")({runid, {"bal", stats_before.max / ((double)nverts/partitions)}});

  sw.reset();
  libpartition::make_domains_overlap(graph, domains, delta);
  lg->debug("make_domains_overlap")({runid, {"time", sw.elapsed().count()}});

  std::vector<double> sizes_after(partitions);
  for (geobin::ID p = 0; p < partitions; p++)
    sizes_after[p] = domains[p].size();
  SimpleStats stats_after;
  compute_stats(sizes_after.data(), partitions, &stats_after);
  print_stats("sizes after", &stats_after, runid);
  lg->debug("bal after")({runid, {"bal", stats_after.max / ((double)nverts/partitions)}});

  std::vector<double> sizes_fractions(partitions);
  for (geobin::ID p = 0; p < partitions; p++)
    sizes_fractions[p] = sizes_after[p] / sizes_before[p];
  SimpleStats stats_frac;
  compute_stats(sizes_fractions.data(), partitions, &stats_frac);
  print_stats("sizes after/before", &stats_frac, runid);

  std::vector<double> part_overlap(nverts, 0);
  for (geobin::ID p = 0; p < partitions; p++)
    for (geobin::ID n : domains[p])
      part_overlap[n] += 1;
  SimpleStats stats_overlap;
  compute_stats(part_overlap.data(), nverts, &stats_overlap);
  print_stats("overlap stats pointwise", &stats_overlap, runid);
  lg->debug("overlap total")({runid, {"overlap_total", stats_after.sum / nverts}});

  lg->info("writing overlapping partition");

  geobin::serialize_domains(output_path.c_str(), domains);

  if (!test_overlap.empty()) {
    std::vector<geobin::ID> fake_partition(graph.vertices.size(), 0);
    for (geobin::ID n : domains[0])
      fake_partition[n] = 1;
    serialize_graph_partition_vtu(graph, fake_partition, test_overlap.c_str());
  }
}
