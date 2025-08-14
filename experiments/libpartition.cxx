#include "libpartition.hxx"

#ifdef HYPERDG_USE_KaHIP
#include <kaHIP_interface.h>
#endif
#ifdef HYPERHDG_USE_METIS
#include <metis.h>
#endif
#ifdef HYPERHDG_USE_PARHIP
#include <parhip_interface.h>
#endif



#include <print>
#include <spdlog/spdlog.h>

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
  const char* backend = config->backend;

  if (strcmp(backend, "naive")) {
    naive_geometric_partition(graph, npartition, partition, edgecut, &config->naive_partitions_z);
  }
  #ifdef HYPERHDG_USE_KaHIP
  else if (strcmp(backend, "kahip") == 0) {
    kaffpa((kidx_t*)&nverts, NULL, (kidx_t*)graph->xadj.data(), NULL, (kidx_t*)graph->adjncy.data(),
           (kidx_t*)npartition, imbalance,
           config->kahip_suppress_output, config->kahip_seed, config->kahip_mode,
           (kidx_t*)edgecut, (idx_t*)partition);
  }
  #endif
  #ifdef HYPERHDG_USE_METIS
  else if (strcmp(backend, "metis")) {
    METIS_PartGraphKway((idx_t*)&nverts, (idx_t*)&nedges, (idx_t*)graph->xadj.data(), (idx_t*)graph->adjncy.data(), NULL, NULL, NULL, (idx_t*)npartition, NULL, NULL, NULL, (idx_t*)edgecut, (idx_t*)partition);
  }
  #endif
  #ifdef HYPERHDG_USE_PARHIP
  else if (strcmp(backend, "parhip")) {
    /*
     * Use a parallel naive partition to distribute the vertices between the processes,
     * then reorder graph such that partitions are consecutive (this must be sequential?)
     */
    MPI_Init(NULL, NULL);
    MPI_Comm* comm = &MPI_COMM_WORLD;
    ParHIPPartitionKWay(idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, idxtype *adjwgt,
      int *nparts, double* imbalance,
      bool suppress_output, int seed, int mode,
      int *edgecut, idxtype *part, comm);
  }
  #endif
  else {
    spdlog::get("logger")->info("unsupported backend");
  }
}

void make_domains_overlap(geobin::Graph& graph, std::vector<std::vector<geobin::ID>>& domains, geobin::ID delta) {
  std::vector<geobin::u8> visited(graph.vertices.size());
  for (geobin::ID p = 0; p < domains.size(); p++) {
    if (domains[p].size() == 0)
      continue;
    std::fill(visited.begin(), visited.end(), 0); // slowest? -> use partition id to track
    for (const geobin::ID& n : domains[p])
      visited[n] = 1;

    // frontier marker
    domains[p].push_back((geobin::ID)-1);

    geobin::ID hop = 0;
    for (geobin::ID bfs_front = 0; bfs_front < domains[p].size() && hop < delta; bfs_front++) {
      const geobin::ID n = domains[p][bfs_front];

      // if we see a frontier marker, then hop is complete
      if (n == (geobin::ID)-1) {
        hop++;
        if (domains[p].back() == (geobin::ID)-1) {
          spdlog::get("logger")->warn("bfs terminated early")({{"rounds_completed", hop-1}, {"domain_idx", p}, {"rounds_requested", delta}});
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
  for (geobin::ID p = 0; p < domains.size(); p++) {
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
}

}
