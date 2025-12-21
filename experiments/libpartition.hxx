#ifndef LIBPARTITION_HXX
#define LIBPARTITION_HXX

#include <type_traits>
#include "geobin.hxx"
#include <cassert>
#include <limits>
#include <cmath>

namespace libpartition {

struct PartConfig {
  const char* backend;
  bool kahip_suppress_output;
  int  kahip_seed;
  int  kahip_mode;
  geobin::ID naive_partitions_z;
};

geobin::ID compute_edgecut(geobin::Graph* graph, geobin::ID* partition);

void naive_geometric_partition(geobin::Graph* graph, geobin::ID* npartition, geobin::ID* partition, geobin::ID* edgecut, geobin::ID* npartitions_z);

int do_partition(geobin::Graph* graph, geobin::ID* npartition, double* imbalance, geobin::ID* partition, geobin::ID* edgecut, PartConfig* config);

void make_domains_overlap(geobin::Graph& graph, std::vector<std::vector<geobin::ID>>& domains, geobin::ID delta);

}


#endif // LIBPARTITION_HXX
