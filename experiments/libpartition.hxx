#ifndef LIBPARTITION_HXX
#define LIBPARTITION_HXX

#include <type_traits>
#include "geobin.hxx"
#include <cassert>
#include <limits>
#include <cmath>
#include <frozen/unordered_map.h>
#include <frozen/string.h>

namespace libpartition {

#define PART_BACKEND_LIST(X) \
    X(KaFFPa) \
    X(METIS) \
    X(NAIVE)

enum class Backend {
#define X(elem) elem,
    PART_BACKEND_LIST(X)
#undef X
    COUNT
};
constexpr auto backend_to_str = frozen::make_unordered_map<Backend, const char*>({
#define X(elem) {Backend::elem, #elem},
  PART_BACKEND_LIST(X)
#undef X
});
constexpr auto str_to_backend = frozen::make_unordered_map<frozen::string, Backend> ({
#define X(elem) {#elem, Backend::elem},
  PART_BACKEND_LIST(X)
#undef X
});

struct PartConfig {
  Backend backend;
  bool kahip_suppress_output;
  int  kahip_seed;
  int  kahip_mode;
  geobin::ID naive_partitions_z;
};

geobin::ID compute_edgecut(geobin::Graph* graph, geobin::ID* partition);

void naive_geometric_partition(geobin::Graph* graph, geobin::ID* npartition, geobin::ID* partition, geobin::ID* edgecut, geobin::ID* npartitions_z);

void do_partition(geobin::Graph* graph, geobin::ID* npartition, double* imbalance, geobin::ID* partition, geobin::ID* edgecut, PartConfig* config);

}


#endif // LIBPARTITION_HXX
