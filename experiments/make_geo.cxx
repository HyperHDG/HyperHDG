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
using u8 = uint8_t;
using Real = double;
using ID = u64;
using ID = u64;

struct Connection {
  ID f1;
  ID f2;
  Real a1;
  Real a2;
};

struct ConnectionPoint {
  ID nodeid;
  Real a;
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


Point interpolate(const Point& u, const Point& v, Real a) {
  Point res;
  for (u64 i = 0; i < 3; i++)
    res[i] = (1-a)*u[i] + a*v[i];
  return res;
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

struct GraphEdgeList {
  std::vector<Edge> edges;
  std::vector<Point> vertices;
  std::vector<Prop> edge_props;
  std::vector<std::pair<u8, u8>> types;
};

GraphEdgeList& compute_types(GraphEdgeList& graph) {
  graph.types.resize(0);

  // Calculate the bounding box (min/max x, y, z) for all vertices
  Real min_x = 1e10, min_y = 1e10, min_z = 1e10, max_x = 1e-10, max_y = 1e-10, max_z = 1e-10;
  for (const Point& vertex : graph.vertices) {
      min_x = std::min(min_x, vertex[0]);
      min_y = std::min(min_y, vertex[1]);
      min_z = std::min(min_z, vertex[2]);
      max_x = std::max(max_x, vertex[0]);
      max_y = std::max(max_y, vertex[1]);
      max_z = std::max(max_z, vertex[2]);
  }

  for (const Edge& edge : graph.edges) {
    u8 left = 0, right = 0;
    Point vertex = graph.vertices[edge.first];
    if (vertex[0] - min_x < 1e-6 * (max_x - min_x) || max_x - vertex[0] < 1e-6 * (max_x - min_x) ||
        vertex[1] - min_y < 1e-6 * (max_y - min_y) || max_y - vertex[1] < 1e-6 * (max_y - min_y))
      left = 1;
    vertex = graph.vertices[edge.second];
    if (vertex[0] - min_x < 1e-6 * (max_x - min_x) || max_x - vertex[0] < 1e-6 * (max_x - min_x) ||
        vertex[1] - min_y < 1e-6 * (max_y - min_y) || max_y - vertex[1] < 1e-6 * (max_y - min_y))
      right = 1;
    graph.types.push_back({left, right});
  }

  return graph;
}

void serialize_txt(const char* output_path, const GraphEdgeList& graph) {
  const std::vector<Point>& vertices = graph.vertices;
  const std::vector<Edge>& edges = graph.edges;
  const std::vector<Prop>& edge_props = graph.edge_props;

  std::ofstream gfile(std::format("{}.geo", output_path));

  std::print(gfile, "# This file was auto-generated!\n\n");
  std::print(gfile, "Space_Dim     = 3;  # Dimension of space.\n");
  std::print(gfile, "HyperEdge_Dim = 1;  # Dimension of hyperedge (must be uniform).\n");
  std::print(gfile, "N_Points      = {};  # Number of vertices.\n", vertices.size());
  std::print(gfile, "N_HyperNodes  = {};  # Number of hypernodes.\n", vertices.size());
  std::print(gfile, "N_HyperEdges  = {};  # Number of hyperedges.\n", edges.size());

  std::print(gfile, "\nPOINTS:\n");
  for (const Point& vertex : vertices)
    std::print(gfile, "{:.18e}  {:.18e}  {:.18e}\n", vertex[0], vertex[1], vertex[2]);

  std::print(gfile, "\nHYPERNODES_OF_HYPEREDGES:\n");
  for (const Edge& edge : edges)
    std::print(gfile, "{} {}\n", edge.first, edge.second);

  std::print(gfile, "\nTYPES_OF_HYPERFACES:\n");
  for (const auto& type : graph.types)
    std::print(gfile, "{} {}\n", type.first, type.second);

  std::print(gfile, "\nPOINTS_OF_HYPEREDGES:\n");
  for (const Edge& edge : edges)
    std::print(gfile, "{} {}\n", edge.first, edge.second);

  std::print(gfile, "\nHYPEREDGE_PROPERTIES: 12\n");
  for (const Prop& prop : edge_props) {
    std::print(gfile, "{:.18e}", prop[0]);
    for (u64 i = 1; i < 12; i++)
      std::print(gfile, "  {:.18e}", prop[i]);
    std::print(gfile, "\n");
  }

  std::ofstream pfile(std::format("{}_points.txt", output_path));
  for (const Point& vertex : vertices)
    std::print(pfile, "{:.18e} {:.18e} {:.18e}\n", vertex[0], vertex[1], vertex[2]);
}

struct DataTable {
  char name[8];
  u64 offset;     // from file start
  u64 size;       // in bytes
};

struct GeoBinHeader {
  char magic[8]; // should contain GEOBINxx
  u64 space_dim;
  u64 hyperedge_dim;
  u64 n_points;
  u64 n_hypernodes;
  u64 n_hyperedges;
  DataTable tables[5];
};

void serialize_bin(const char* output_path, const GraphEdgeList& graph) {
  u64 n = graph.vertices.size();
  u64 m = graph.edges.size();
  if (graph.edge_props.size() != m) {
    std::println(stderr, "WARNING: unequal number of graph edges and edge props provided {}!={}", m, graph.edge_props.size());
  }

  DataTable points = {
    .name = "POINTS\0", // extra \0
    .offset = sizeof(GeoBinHeader),
    .size = n * 3 * sizeof(Real),
  };

  DataTable hypernodes_of_hyperedges = {
    .name = "HYPNODE", // extra \0
    .offset = points.offset + points.size,
    .size = m * 2 * sizeof(ID),
  };

  DataTable types_of_hyperfaces = {
    .name = "TYPES\0\0", // extra \0
    .offset = hypernodes_of_hyperedges.offset + hypernodes_of_hyperedges.size,
    .size = m * 2 * sizeof(u8),
  };

  DataTable points_of_hyperedges = {
    .name = "POIHYPE", // extra \0
    .offset = types_of_hyperfaces.offset + types_of_hyperfaces.size,
    .size = m * 2 * sizeof(ID),
  };

  DataTable hyperedge_properties = {
    .name = "HYPPROP", // extra \0
    .offset = points_of_hyperedges.offset + points_of_hyperedges.size,
    .size = m * 12 * sizeof(Real),
  };

   GeoBinHeader header = {
    .magic = "GEOBIN1", // extra \0
    .space_dim = 3,
    .hyperedge_dim = 1,
    .n_points = n,
    .n_hypernodes = n,
    .n_hyperedges = m,
    .tables = {
      points,
      hypernodes_of_hyperedges,
      types_of_hyperfaces,
      points_of_hyperedges,
      hyperedge_properties
    },
  };

  std::ofstream file(output_path, std::ios::binary);
  file.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  file.write((char*)&header, sizeof(header));
  assert((u64)file.tellp() == points.offset);
  file.write((char*)&graph.vertices[0], points.size);
  assert((u64)file.tellp() == hypernodes_of_hyperedges.offset);
  file.write((char*)&graph.edges[0], hypernodes_of_hyperedges.size);
  assert((u64)file.tellp() == types_of_hyperfaces.offset);
  file.write((char*)&graph.types[0], types_of_hyperfaces.size);
  assert((u64)file.tellp() == points_of_hyperedges.offset);
  file.write((char*)&graph.edges[0], points_of_hyperedges.size);
  assert((u64)file.tellp() == hyperedge_properties.offset);
  file.write((char*)&graph.edge_props[0], hyperedge_properties.size);
}

GraphEdgeList deserialize_bin(const char* input_path) {
  GraphEdgeList graph;

  std::ifstream file(input_path, std::ios::binary);
  file.exceptions(std::ifstream::badbit | std::ifstream::failbit);
  std::array<char, sizeof(GeoBinHeader)> header_buf;
  file.read(header_buf.data(), sizeof(GeoBinHeader));
  GeoBinHeader* header = (GeoBinHeader*)header_buf.data();
  assert(0 == std::strcmp(header->magic, "GEOBIN1")); // extra trailing \0
  assert(3 == header->space_dim);
  assert(1 == header->hyperedge_dim);

  logi("deserialize_bin");
  logi("  n_points = {}", header->n_points);
  logi("  n_hypernodes = {}", header->n_hypernodes);
  logi("  n_hyperedges = {}", header->n_hyperedges);

  graph.vertices.resize(header->n_points);
  graph.edges.resize(header->n_hyperedges);
  graph.edge_props.resize(header->n_hyperedges);
  graph.types.resize(header->n_hyperedges);

  for (u64 i = 0; i < sizeof(header->tables)/sizeof(DataTable); i++)
    logi("  table[{}].name = {}", i, header->tables[i].name);

  DataTable* tables = header->tables;
  assert((u64)file.tellg() == tables[0].offset);
  file.read((char*)graph.vertices.data(), tables[0].size);
  assert((u64)file.tellg() == tables[1].offset);
  file.read((char*)graph.edges.data(), tables[1].size);
  assert((u64)file.tellg() == tables[2].offset);
  file.read((char*)graph.types.data(), tables[2].size);
  assert((u64)file.tellg() == tables[3].offset);
  file.seekg(tables[3].size, std::ios::cur); // skip next
  assert((u64)file.tellg() == tables[4].offset);
  file.read((char*)graph.edge_props.data(), tables[4].size);

  return graph;
}

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
  std::vector<std::vector<ConnectionPoint>> fiber_connections(fibers.size());

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
      p2[i] = (1-con.a2) * e21[i] + con.a2*e22[i];
    }

    u64 p1id = pcloud.pts.size();
    pcloud.pts.push_back(p1);
    u64 p2id = pcloud.pts.size();
    pcloud.pts.push_back(p2);

    fiber_connections[con.f1].push_back({p1id, con.a1});
    fiber_connections[con.f2].push_back({p2id, con.a2});
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

  logi("generating vertex/edge list");

  std::vector<Point> vertices;
  std::vector<u64> is_in_vertices(pcloud.pts.size(), 0);
  std::vector<u64> vertices_idx(pcloud.pts.size(), (u64)-1);

  std::vector<Edge> edges;

  using Neighbor = nanoflann::ResultItem<uint32_t, Real>; // Neighbor = (id,distance)
  std::vector<Neighbor> neighbors;
  neighbors.reserve(1000); // reserve more than enough

  auto find_merged_id = [&kdtree, &is_in_vertices, &neighbors, &vertices_idx](const Point& p, Real r) {
    // NOTE: nanoflann works with squared L2 distances
    Real d = 2*r;
    const u64 num_neighbors = kdtree.radiusSearch(p.data(), r*r, neighbors); // only consider neighbors closer than r
    u64 index = (u64)-1;
    assert(num_neighbors != 0); // point_a should be found
    for (u64 i = 0; i < num_neighbors; i++) {
      if (0 == is_in_vertices[neighbors[i].first]) // only consider neighbors already in vertices
        continue;
      if (neighbors[i].second < d) {               // of those, pick the closest
        index = vertices_idx[neighbors[i].first];
        d = neighbors[i].second;
      }
    }
    return index;
  };

  auto find_merged_id_and_insert = [&kdtree, &neighbors, &vertices, &is_in_vertices, &vertices_idx, &find_merged_id](u64 id, const Point& p, Real r) {
    u64 index = find_merged_id(p, r);
    if (index == (u64)-1) {
      index = vertices.size();
      vertices.push_back(p);
      is_in_vertices[id] = 1;
    } else {
      is_in_vertices[id] = 0;
      // logi("{} -> {}", id, index);
    }
    vertices_idx[id] = index;

    return index;
  };

  logi("  from connections");

  // if len(vertices) == 0: vertices = np.vstack((point_a, point_b))
  vertices.push_back(pcloud.pts[0]);
  vertices.push_back(pcloud.pts[1]);
  vertices_idx[0] = 0;
  vertices_idx[1] = 1;
  is_in_vertices[0] = 1;
  is_in_vertices[1] = 1;
  edges.push_back({0,1});

  for (u64 cid = 1; cid < connections.size(); cid++) {
    const Point& point_a = pcloud.pts[2*cid];
    const Point& point_b = pcloud.pts[2*cid+1];

    Real r = 1e-10;

    // index_a = np.argmin(np.linalg.norm(point_a - vertices, axis=1))
    // if np.linalg.norm(point_a - vertices[index_a]) > 1e-10:
    //   index_a = len(vertices)
    //   vertices = np.vstack((vertices, point_a))

    const u64 index_a = find_merged_id_and_insert(2*cid, point_a, r);

    // index_b = np.argmin(np.linalg.norm(point_b - vertices, axis=1))
    // if np.linalg.norm(point_b - vertices[index_b]) > 1e-10:
    //   index_b = len(vertices)
    //   vertices = np.vstack((vertices, point_b))

    const u64 index_b = find_merged_id_and_insert(2*cid+1, point_b, r);

    if (index_a != index_b)
      edges.push_back({index_a, index_b});
  }

  logi("  from fibers");

  // NOTE: the python code below does not respect the self loops filtered out above
  //   edges_prop = np.vstack((connectionsProp, np.array(edges_prop)))
  // but we replicate this behaviour.
  std::vector<Prop>& edge_props = connection_props;

  for (u64 fid = 0; fid < fibers.size(); fid++) {
    // helper = act_fibers[act_fibers[:,0] == index, 1]
    std::vector<ConnectionPoint>& helper = fiber_connections[fid];
    if (helper.empty())
      continue;

    const Point& point_a = nodes[fibers[fid].first];
    const Point& point_b = nodes[fibers[fid].second];
    const u64 off = 2*connections.size();

    const u64 endpoint_idx_a = find_merged_id_and_insert(off+fibers[fid].first, point_a, 1e-15);
    const u64 endpoint_idx_b = find_merged_id_and_insert(off+fibers[fid].second, point_b, 1e-15);

    std::sort(helper.begin(), helper.end(), [](ConnectionPoint p1, ConnectionPoint p2){ return p1.a < p2.a; });

    // helper = [0.] + helper + [1.]
    helper.insert(helper.begin(), {.nodeid = endpoint_idx_a, .a = 0});
    helper.insert(helper.end(),   {.nodeid = endpoint_idx_b, .a = 1});

    for (u64 k = 0; k < helper.size()-1; k++) {
      Point point_ab = interpolate(point_a, point_b, helper[k].a);
      Point point_ba = interpolate(point_a, point_b, helper[k+1].a);

      u64 index_a = find_merged_id(point_ab, 1e-10);
      if (index_a == (u64)-1) {
        std::println(stderr, "ERROR: expected connection point to be found");
        return 1;
      }
      u64 index_b = find_merged_id(point_ba, 1e-10);
      if (index_b == (u64)-1) {
        std::println(stderr, "ERROR: expected connection point to be found");
        return 1;
      }

      if (index_a != index_b) {
        edges.push_back({index_a, index_b});
        edge_props.push_back(fiber_props[fid]);
      }
    }
  }

  logi("output size");
  logi("  vertices.size = {}", vertices.size());
  logi("  edges.size    = {}", edges.size());

  if (fs::is_directory(output_path))
    output_path = std::format("{}/fiber_network_{}", output_path, edges.size());

  GraphEdgeList graph = { .edges = edges, .vertices = vertices,  .edge_props = edge_props, .types = {}};
  compute_types(graph);

  std::string binpath = std::format("{}.geo.bin", output_path);
  serialize_bin(binpath.c_str(), graph);
  GraphEdgeList graph2 = deserialize_bin(binpath.c_str());

  logi("generating txt");

  serialize_txt(output_path.c_str(), graph);
  serialize_txt(std::format("{}.geo2", output_path).c_str(), graph2);
}
