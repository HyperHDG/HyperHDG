#include <print>
#include <vector>
#include <cstdint>
#include <utility> // for pair
#include <string>
#include <fstream>
#include <format>
#include <cstdio>
#include <cstring>
#include <nanoflann.hpp>
#include <algorithm>
#include <fmtlog/fmtlog.h>

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


int main(int argc, char** argv) {
  if (argc < 3) {
    std::println(stderr, "usage: {} <input_folder> <output_folder>", argv[0]);
    return 1;
  }

  logi("hello world from fmtlog {}", 69420);

  // argv[0] is executable path
  const char* input_folder = argv[1];
  const char* output_folder = argv[2];

  // read data
  PointCloud<Real> nodes;

  nodes.pts = read_nodes(std::format("{}/nodes.csv", input_folder).c_str());
  logi("nodes: {}", nodes.pts.size());

  std::vector<Edge> fibers = read_fibers(std::format("{}/fibers.csv", input_folder).c_str());
  logi("fibers: {}", fibers.size());

  std::vector<Connection> connections = read_connections(std::format("{}/connections.csv", input_folder).c_str());
  logi("connections: {}", connections.size());

  std::vector<Prop> fiber_props = read_props(std::format("{}/fibersProps.csv", input_folder).c_str());
  logi("fibersProps: {}", fiber_props.size());

  std::vector<Prop> connection_props = read_props(std::format("{}/connectionsProp.csv", input_folder).c_str());
  logi("connectionProps: {}", connection_props.size());

  // first enter all point to be considered

  for (const Connection& con : connections) {
    Edge  e1 = fibers[con.f1], e2 = fibers[con.f2];
    Point p1, p2, e11 = nodes.pts[e1.first], e12 = nodes.pts[e1.second], e21 = nodes.pts[e2.first], e22 = nodes.pts[e2.second];

    for (u64 i = 0; i < 3; i++) {
      p1[i] = (1-con.a1) * e11[i] + con.a1*e12[i];
      p2[i] = (1-con.a1) * e21[i] + con.a1*e22[i];
    }

    nodes.pts.push_back(p1);
    nodes.pts.push_back(p2);
  }

  const u64 dim = 3, maxleaf = 10;
  using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
    nanoflann::L2_Simple_Adaptor<Real, PointCloud<Real>>,
    PointCloud<Real>,
    dim
  >;
  KDTree kdtree(dim, nodes, {maxleaf});
  logi("building kdtree");
  kdtree.buildIndex();

  logi("querying kdtree");
  u64 avg = 0, max = 0;
  Real r = 1e-10;
  using Neighbor = nanoflann::ResultItem<uint32_t, Real>;
  std::vector<Neighbor> neighbors;
  neighbors.reserve(1000); // reserve more than enough
  for (const Point& node : nodes.pts) {
    const u64 n = kdtree.radiusSearch(&node[0], r, neighbors);
    avg += n;
    max = std::max(max, n);
  }

  logi("num neighbors: avg = {:.3f}, max = {}",
       (double)avg / nodes.pts.size(), max);

}
