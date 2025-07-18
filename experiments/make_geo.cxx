#include <filesystem>
#include <print>
#include <vector>
#include <cstdint>
#include <utility> // for pair
#include <string>
#include <fstream>
#include <format>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include <nanoflann.hpp>
#include <fmtlog/fmtlog.h>

namespace {

using u64 = uint64_t;
using Real = double;
using ID = u64;
using ID = u64;

struct Connection {
  ID f1;
  ID f2;
  Real a1;
  Real a2;
};

using Point = std::array<Real, 3>;
using Edge = std::pair<ID, ID>;
using Prop = std::array<Real, 12>;

std::vector<Point> read_nodes(const char* path) {
  std::vector<Point> nodes;
  std::fstream nodes_file(path);
  if (!nodes_file) {
    std::println("error: couldn't open {}", path);
    return {};
  }
  for (std::string line; std::getline(nodes_file, line); ) {
    ID id; Real x, y, z;
    if (4 == sscanf(line.c_str(), "%zu,%lf,%lf,%lf", &id, &x, &y, &z)) {
      nodes.push_back({x,y,z});
    } else {
      if (strcmp("Id,x,y,z", line.c_str()) == 0)
        continue;
      else {
        std::println(stderr, "error: {}: couldn't parse line '{}'", path, line);
        return {};
      }
      continue;
    }
  }
  return nodes;
}

std::vector<Edge> read_fibers(const char* path) {
  std::vector<Edge> fibers;
  std::fstream fibers_file(path);
  if (!fibers_file) {
    std::println(stderr, "error: couldn't open {}", path);
    return {};
  }
  for (std::string line; std::getline(fibers_file, line); ) {
    ID f; ID u, v;
    if (3 == sscanf(line.c_str(), "%zu,%zu,%zu", &f, &u, &v)) {
      fibers.push_back({u,v});
    } else {
      if (strcmp("Id,node1,node2", line.c_str()) == 0)
        continue;
      else {
        std::println(stderr, "error: {}: couldn't parse line '{}'", path, line);
        return {};
      }
      continue;
    }
  }

  return fibers;
}

std::vector<Connection> read_connections(const char* path) {
  std::vector<Connection> connections;
  std::fstream connections_file(path);
  if (!connections_file) {
    std::println(stderr, "error: couldn't open {}", path);
    return {};
  }
  for (std::string line; std::getline(connections_file, line); ) {
    ID c,f1,f2; Real a1,a2;
    if (5 == sscanf(line.c_str(), "%zu,%zu,%zu,%lf,%lf", &c, &f1, &f2, &a1, &a2)) {
      connections.push_back({f1,f2,a1,a2});
    } else {
      if (strcmp("Id,fiber1,fiber2,a1,a2", line.c_str()) == 0)
        continue;
      else {
        std::println(stderr, "error: {}: couldn't parse line '{}'", path, line);
        return {};
      }
      continue;
    }
  }

  return connections;
}

std::vector<Prop> read_props(const char* path) {
  std::vector<Prop> props;
  std::fstream fiber_props_file(path);
  if (!fiber_props_file) {
    std::println(stderr, "error: couldn't open {}", path);
    return {};
  }
  for (std::string line; std::getline(fiber_props_file, line); ) {
    ID f;
    Prop prop; // 6 structural + 2*3 normal
    int total_chars_read = 0, chars_read = 0;
    const char* cline = line.c_str();

    if (1 != sscanf(cline, "%zu%n,", &f, &chars_read)) {
      if (0 == strcmp(cline, "Id,EA,kG_1A,kG_2A,G_xI_x,E_1I_1,E_2I_2,n_11,n_12,n_13,n_21,n_22,n_23"))
        continue;
      std::println(stderr, "error: {}: sscanf format '%zu' invalid for '{}'", path, cline);
      return {};
    }

    total_chars_read += chars_read+1;

    for (u64 i = 0; i < 12; i++) {
      if (1 != sscanf(cline+total_chars_read, "%lf%n", &prop[i], &chars_read)) {
        std::println(stderr, "error: {}: sscanf format '%lf' invalid for '{}'", path, cline+total_chars_read);
        return {};
      }
      total_chars_read += chars_read+1; // skip ','
    }

    props.push_back(prop);
  }

  return props;
}

template <typename T>
struct PointCloud
{
  using Point = std::array<T,3>;
  using coord_t = T;

  std::vector<Point> pts;

  inline size_t kdtree_get_point_count() const { return pts.size(); }
  inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
    return pts[idx][dim];
  }

  template <class BBOX>
  bool kdtree_get_bbox(BBOX& /* bb */) const {
      return false;
  }
};

}

int main(int argc, char** argv) {
  fmtlog::startPollingThread(1);

  // verify args
  // argv[0] is executable path
  if (argc < 3) {
    std::println(stderr, "usage: {} <input_folder> <output_folder>", argv[0]);
    return 1;
  }
  std::string input_folder = argv[1];
  std::string output_path = argv[2];
  namespace fs = std::filesystem;

  if (!fs::is_directory(input_folder)) {
    loge("invalid argument <input_folder>, got '{}'", input_folder);
    return 1;
  }

  if (!fs::is_directory(output_path)) {
    std::ofstream file(std::format("{}.geo", output_path));
    if (!file) {
      loge("invalid argument <output_path>, got '{}'", output_path);
      loge("  neither directory, nor writable path");
      return 1;
    }
  }

  logi("reading data");

  std::vector<Point> nodes = read_nodes(std::format("{}/nodes.csv", input_folder).c_str());
  logi("nodes: {}", nodes.size());

  std::vector<Edge> fibers = read_fibers(std::format("{}/fibers.csv", input_folder).c_str());
  logi("fibers: {}", fibers.size());

  std::vector<Connection> connections = read_connections(std::format("{}/connections.csv", input_folder).c_str());
  logi("connections: {}", connections.size());

  std::vector<Prop> fiber_props = read_props(std::format("{}/fibersProps.csv", input_folder).c_str());
  logi("fibersProps: {}", fiber_props.size());

  std::vector<Prop> connection_props = read_props(std::format("{}/connectionsProp.csv", input_folder).c_str());
  logi("connectionProps: {}", connection_props.size());


  // collect all points to build fast KNN lookup datastructure
  PointCloud<Real> pcloud;
  pcloud.pts.reserve(connections.size()*2 + nodes.size());

  // mapping of fiber to nodeids of connection points
  std::vector<std::vector<u64>> fiber_segments(fibers.size());

  // first all connection points
  for (const Connection& con : connections) {
    Edge  e1 = fibers[con.f1], e2 = fibers[con.f2];
    Point p1, p2;
    Point e11 = nodes[e1.first];
    Point e12 = nodes[e1.second];
    Point e21 = nodes[e2.first];
    Point e22 = nodes[e2.second];

    for (u64 i = 0; i < 3; i++) {
      p1[i] = (1-con.a1) * e11[i] + con.a1*e12[i];
      p2[i] = (1-con.a1) * e21[i] + con.a1*e22[i];
    }

    u64 p1id = pcloud.pts.size();
    pcloud.pts.push_back(p1);
    u64 p2id = pcloud.pts.size();
    pcloud.pts.push_back(p2);

    fiber_segments[con.f1].push_back(p1id);
    fiber_segments[con.f2].push_back(p2id);
  }

  // next all fiber endpoints
  for (const Point& p : nodes)
    pcloud.pts.push_back(p);

  const u64 dim = 3, maxleaf = 10;
  using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<Real, PointCloud<Real>>,
    PointCloud<Real>,
    dim
  >;
  KDTree kdtree(dim, pcloud, {maxleaf});
  logi("building kdtree");
  kdtree.buildIndex();

  logi("merging close points");

  // create mapping from node id to node id with which to merge
  // init to 0 to not skip first (and thus all) node(s)
  std::vector<u64> merge_map(pcloud.pts.size(), (u64)-1);

  u64 nn_avg = 0, nn_max = 0;
  Real r = 1e-10;
  using Neighbor = nanoflann::ResultItem<uint32_t, Real>; // Neighbor = (id,distance)
  std::vector<Neighbor> neighbors;
  neighbors.reserve(1000); // reserve more than enough
  for (u64 nodeid = 0; nodeid < pcloud.pts.size(); nodeid++) {
    // else find close neighbors
    const double* p = pcloud.pts[nodeid].data();
    const u64 num_neighbors = kdtree.radiusSearch(p, r, neighbors);

    // if there were none, continue
    if (num_neighbors == 0)
      continue;

    // else map myself to neighbor with lowest idx
    u64 newid = nodeid;
    for (const Neighbor& neighbor : neighbors)
      newid = std::min(newid, (u64)neighbor.first);
    merge_map[nodeid] = newid;

    // count avg, max number of neighbors
    nn_avg += num_neighbors;
    nn_max = std::max(nn_max, num_neighbors);
  }

  logi("num neighbors: avg = {:.3f}, max = {}",
       (double)nn_avg / pcloud.pts.size(), nn_max);

  // create mapping from old node ids to consecutive new ids after merging
  u64 new_count = 0;
  std::vector<u64> new_ids(pcloud.pts.size(), 0);
  for (u64 old_id = 0; old_id < pcloud.pts.size(); old_id++) {
    // node is kept
    if (merge_map[old_id] == old_id) {
      new_ids[old_id] = new_count;
      new_count++;
    } else {
      u64 merge_partner_old = merge_map[old_id];
      new_ids[old_id] = new_ids[merge_partner_old];
    }
  }

  // create new node list
  std::vector<Point> vertices(new_count);
  for (u64 old_id = 0; old_id < pcloud.pts.size(); old_id++)
    vertices[new_ids[old_id]] = pcloud.pts[old_id];

  logi("create edge list");

  std::vector<Edge> edges;
  std::vector<Prop> edge_props;

  edges.reserve(connections.size()*2+fibers.size());

  // for every connection with id two nodes are created with 2*id, 2*id+1 which form the connection
  for (u64 cid = 0; cid < connections.size(); cid++) {
    u64 nid1 = new_ids[2*cid+0];
    u64 nid2 = new_ids[2*cid+1];
    // ignore self loops (TODO: multiedges?)
    if (nid1 != nid2) {
      edges.push_back({nid1, nid2});
      edge_props.push_back(connection_props[cid]);
    }
  }

  // add segmentation of fiber as edges
  for (u64 fid = 0; fid < fibers.size(); fid++) {
    Edge fiber = fibers[fid];
    const std::vector<u64>& cpoints = fiber_segments[fid];

    // TODO: test for multi-edges?
    auto add_edge_if_simple = [&](u64 u, u64 v){
      u64 nu = new_ids[u];
      u64 nv = new_ids[v];
      if (nu != nv) {
        edges.push_back({nu,nv});
        edge_props.push_back(fiber_props[fid]);
      }
    };

    if (cpoints.empty()) {
      add_edge_if_simple(fiber.first,fiber.second);
    } else {
      add_edge_if_simple(fiber.first, cpoints[0]);
      for (u64 i = 1; i < cpoints.size(); i++)
        add_edge_if_simple(cpoints[i-1], cpoints[i]);
      add_edge_if_simple(cpoints.back(), fiber.second);
    }
  }

  logi("generate output");

  // Calculate the bounding box (min/max x, y, z) for all vertices
  Real min_x = 1e10, min_y = 1e10, min_z = 1e10, max_x = 1e-10, max_y = 1e-10, max_z = 1e-10;
  for (const Point& vertex : vertices) {
      min_x = std::min(min_x, vertex[0]);
      min_y = std::min(min_y, vertex[1]);
      min_z = std::min(min_z, vertex[2]);
      max_x = std::max(max_x, vertex[0]);
      max_y = std::max(max_y, vertex[1]);
      max_z = std::max(max_z, vertex[2]);
  }

  if (fs::is_directory(output_path))
    output_path = std::format("{}/fiber_network_{}", output_path, edges.size());

  std::ofstream gfile(std::format("{}.geo", output_path));

  std::print(gfile, "# This file was auto-generated!\n\n");
  std::print(gfile, "Space_Dim     = 3;  # Dimension of space.\n");
  std::print(gfile, "HyperEdge_Dim = 3;  # Dimension of hyperedge (must be uniform).\n");
  std::print(gfile, "N_Points      = {}; # Number of vertices.\n", vertices.size());
  std::print(gfile, "N_HyperNodes  = {}; # Number of hypernodes.\n", vertices.size());
  std::print(gfile, "N_HyperEdges  = {}; # Number of hyperedges.\n", edges.size());

  std::print(gfile, "\nPOINTS:\n");
  for (const Point& vertex : vertices)
    std::print(gfile, "{} {} {}\n", vertex[0], vertex[1], vertex[2]);

  std::print(gfile, "\nHYPERNODES_OF_HYPEREDES:\n");
  for (const Edge& edge : edges)
    std::print(gfile, "{} {}\n", edge.first, edge.second);

  std::print(gfile, "\nTYPES_OF_HYPERFACES:\n");
  for (const Edge& edge : edges) {
    u64 left = 0, right = 0;
    Point& vertex = vertices[edge.first];
    if (vertex[0] - min_x < 1e-6 * (max_x - min_x) || max_x - vertex[0] < 1e-6 * (max_x - min_x) ||
        vertex[1] - min_y < 1e-6 * (max_y - min_y) || max_y - vertex[0] < 1e-6 * (max_y - min_y))
      left = 1;
    vertex = vertices[edge.second];
    if (vertex[0] - min_x < 1e-6 * (max_x - min_x) || max_x - vertex[0] < 1e-6 * (max_x - min_x) ||
        vertex[1] - min_y < 1e-6 * (max_y - min_y) || max_y - vertex[0] < 1e-6 * (max_y - min_y))
      right = 1;
    std::print(gfile, "{} {}", left, right);
  }

  std::println(gfile, "\nHYPEREDGE_PROPERTIES: 12\n");
  for (const Prop& prop : edge_props) {
    std::print(gfile, "{}", prop[0]);
    for (u64 i = 1; i < 12; i++)
      std::print(gfile, " {}", prop[i]);
    std::print(gfile, "\n");
  }

  std::ofstream pfile(std::format("{}_points.txt", output_path));
  for (const Point& vertex : vertices)
    std::print(pfile, "{} {} {}\n", vertex[0], vertex[1], vertex[2]);
}
