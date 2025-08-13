#include "libpartition.hxx"

#include <kaHIP_interface.h>
#include <metis.h>
#include <print>

namespace libpartition {

geobin::ID compute_edgecut(geobin::Graph* graph, geobin::ID* partition) {
  geobin::ID edgecut = 0;
  for (const geobin::Edge& edge : graph->edges) {
    if (partition[edge.first] != partition[edge.second])
      edgecut++;
  }
  return edgecut;
}

void naive_geometric_partition(geobin::Graph* graph, geobin::ID* npartition, geobin::ID* partition, geobin::ID* edgecut, geobin::ID* npartitions_z) {
  geobin::Point min_p = {std::numeric_limits<geobin::Real>::max()}, max_p = {std::numeric_limits<geobin::Real>::min()};
  for (geobin::u64 n = 0; n < graph->vertices.size(); n++) {
    const geobin::Point& p = graph->vertices[n];
    for (geobin::u64 i = 0; i < 3; i++) {
      min_p[i] = std::min(min_p[i], p[i]);
      max_p[i] = std::max(max_p[i], p[i]);
    }
  }

  geobin::Real eps = 1e-10;
  geobin::ID npart = *npartition;
  geobin::ID npart_z = *npartitions_z;
  geobin::ID npart_xy = (geobin::ID)std::sqrt(npart/npart_z);
  std::array<geobin::ID, 3> partitions3d = {npart_xy, npart_xy, npart_z};

  for (geobin::u64 n = 0; n < graph->vertices.size(); n++) {
    const geobin::Point& p = graph->vertices[n];
    geobin::ID pid = 0;
    for (geobin::u64 i = 0; i < 3; i++) {
      pid *= partitions3d[i];
      pid += p[i]/((1+eps)*(max_p[i]-min_p[i])) * partitions3d[i]; // truncate
    }
    assert(pid < npart);
    partition[n] = pid;
  }

  *edgecut = compute_edgecut(graph, partition);
}

void do_partition(geobin::Graph* graph, geobin::ID* npartition, double* imbalance, geobin::ID* partition, geobin::ID* edgecut, PartConfig* config) {
  geobin::ID nverts = graph->vertices.size();
  geobin::ID nedges = graph->vertices.size();
  using kidx_t = int;
  static_assert(sizeof(kidx_t) == sizeof(geobin::ID));

  switch (config->backend) {
    case Backend::KaFFPa:
      kaffpa((kidx_t*)&nverts, NULL, (kidx_t*)graph->xadj.data(), NULL, (kidx_t*)graph->adjncy.data(), (kidx_t*)npartition, imbalance, config->kahip_suppress_output, config->kahip_seed, config->kahip_mode, (kidx_t*)edgecut, (idx_t*)partition);
      break;
    case Backend::METIS:
      METIS_PartGraphKway((idx_t*)&nverts, (idx_t*)&nedges, (idx_t*)graph->xadj.data(), (idx_t*)graph->adjncy.data(), NULL, NULL, NULL, (idx_t*)npartition, NULL, NULL, NULL, (idx_t*)edgecut, (idx_t*)partition);
      break;
    case Backend::NAIVE:
      naive_geometric_partition(graph, npartition, partition, edgecut, &config->naive_partitions_z);
      break;
    default:
      std::println("ERROR: unsupported backend '{}'", backend_to_str.at(config->backend));
      break;
  }
}



}
