#include <print>
#include <fstream>
#include <limits>

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

void print_stats(const char* msg, SimpleStats* stats) {
  auto lg = spdlog::get("logger");
  lg->info(msg)({
      {"min", stats->min}, {"max", stats->max}, {"sum", stats->sum}, {"avg", stats->avg}, {"std", stats->stddev}
  });
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

  std::string backend = "KaFFPa";
  app.add_option("-b,--backend", backend, "the partitioner backend to use, must be one of (kahip|metis|naive)");

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

  std::string log_file;
  app.add_option("--log-file", log_file, "log to file path instead of to stdout");

  std::string log_level = "info";
  app.add_option("--log-level", log_level, "set log level, one of  (trace|debug|info|warn|err)");

  CLI11_PARSE(app,argc,argv);

  auto lg = log_file.empty() ? spdlog::stdout_logger_st("logger") : spdlog::basic_logger_st("logger", log_file);
  spdlog::stopwatch sw;
  spdlog::cfg::helpers::load_levels(log_level);

  lg->info("args")({
      {"input_path", input_path},
      {"partitions", partitions},
      {"output_path", output_path},
      {"vtu_output_path", vtu_output_path},
      {"delta", delta},
      {"backend", backend},
      {"test_overlap", test_overlap}
  });

  lg->info("reading graph...");

  sw.reset();
  geobin::Graph graph = geobin::deserialize_bin(input_path.c_str());
  geobin::ID nverts = graph.vertices.size();
  geobin::ID nedges = graph.edges.size();

  lg->debug("reading graph")({{"time", sw.elapsed().count()}});
  lg->debug("graph stats")({{"nodes", nverts}, {"edges", nedges}});

  // TODO: set default hops based on some graph stats like diameter, girth, etc
  //   maybe using 2-approximation or 3/2-approximation of diameter
  //   NOTE2: or maybe not? maybe 2 is fine?

  lg->info("partitioner backend...");

  sw.reset();
  geobin::ID edgecut;
  std::vector<geobin::ID> partition(nverts);
  libpartition::PartConfig config = {libpartition::str_to_backend.at(frozen::string(backend)), !kahip_no_suppress_output, kahip_seed, kahip_mode, partitions_z};
  libpartition::do_partition(&graph, &partitions, &imbalance, partition.data(), &edgecut, &config);

  lg->debug("partitioner_backend")({{"time", sw.elapsed().count()}});

  std::vector<std::vector<geobin::ID>> domains(partitions);
  for (geobin::ID n = 0; n < nverts; n++)
    domains[partition[n]].push_back(n);

  std::vector<double> sizes_before(partitions);
  for (geobin::ID p = 0; p < partitions; p++)
    sizes_before[p] = domains[p].size();
  SimpleStats stats_before;
  compute_stats(sizes_before.data(), partitions, &stats_before);

  lg->debug("partition metrics")({{"cut", edgecut}, {"bal", stats_before.max / ((double)nverts/partitions)}});

  if (!vtu_output_path.empty()) {
    lg->info("serializing graph partition to vtu");
    serialize_graph_partition_vtu(graph, partition, vtu_output_path.c_str());
  }

  lg->info("make partition overlap and filter boundary nodes...");

  print_stats("sizes before", &stats_before);

  sw.reset();
  libpartition::make_domains_overlap(graph, domains, delta);
  lg->debug("make_domains_overlap")({{"time", sw.elapsed().count()}});

  std::vector<double> sizes_after(partitions);
  for (geobin::ID p = 0; p < partitions; p++)
    sizes_after[p] = domains[p].size();
  SimpleStats stats_after;
  compute_stats(sizes_after.data(), partitions, &stats_after);
  print_stats("sizes after", &stats_after);
  lg->debug("bal after")({{"bal", stats_after.max / ((double)nverts/partitions)}});

  std::vector<double> sizes_fractions(partitions);
  for (geobin::ID p = 0; p < partitions; p++)
    sizes_fractions[p] = sizes_after[p] / sizes_before[p];
  SimpleStats stats_frac;
  compute_stats(sizes_fractions.data(), partitions, &stats_frac);
  print_stats("sizes after/before", &stats_frac);

  std::vector<double> part_overlap(nverts, 0);
  for (geobin::ID p = 0; p < partitions; p++)
    for (geobin::ID n : domains[p])
      part_overlap[n] += 1;
  SimpleStats stats_overlap;
  compute_stats(part_overlap.data(), nverts, &stats_overlap);
  print_stats("overlap stats pointwise", &stats_overlap);
  lg->debug("overlap total")({{"overlap_total", stats_after.sum / nverts}});

  lg->info("writing overlapping partition");

  geobin::serialize_domains(output_path.c_str(), domains);

  if (!test_overlap.empty()) {
    std::vector<geobin::ID> fake_partition(graph.vertices.size(), 0);
    for (geobin::ID n : domains[0])
      fake_partition[n] = 1;
    serialize_graph_partition_vtu(graph, fake_partition, test_overlap.c_str());
  }
}
