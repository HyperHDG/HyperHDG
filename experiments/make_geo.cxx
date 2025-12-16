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

#include "geobin.hxx"

using namespace geobin;

#define logi(...) std::println(stdout, __VA_ARGS__)
#define loge(...) std::println(stderr, __VA_ARGS__)

Point interpolate(const Point& u, const Point& v, Real a) {
  Point res;
  for (ID i = 0; i < 3; i++)
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

int usage(int argc, char** argv) {
  fprintf(stderr, "ERROR: usage: %s <input_folder> <output_folder> [txt]", argv[0]);
  return 1;
}

int main(int argc, char** argv) {
  if (argc < 3) return usage(argc, argv);

  std::string input_folder = argv[1];
  std::string output_path = argv[2];
  namespace fs = std::filesystem;

  if (!fs::is_directory(input_folder)) {
    loge("invalid argument <input_folder>, got '{}'", input_folder);
    return usage(argc, argv);
  }

  if (!fs::is_directory(output_path)) {
    std::string test_path = std::format("{}.test", output_path);
    std::ofstream file(test_path);
    if (!file) {
      loge("invalid argument <output_path>, got '{}'", output_path);
      loge("  neither directory, nor writable path");
      return usage(argc, argv);
    }
    fs::remove(test_path);
  }

  bool txt = false;
  if (argc == 4 && 0 == strcmp(argv[3], "txt")) {
    txt = true;
  }

  bool testbin = false;
  if (argc == 4 && 0 == strcmp(argv[3], "testbin")) {
    txt = true;
  }

  logi("reading data");

  std::vector<Point> nodes = read_nodes(std::format("{}/nodes.csv", input_folder).c_str());
  logi("nodes: {}", nodes.size());

  std::vector<geobin::Edge> fibers = read_fibers(std::format("{}/fibers.csv", input_folder).c_str());
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
    geobin::Edge  e1 = fibers[con.f1], e2 = fibers[con.f2];
    Point p1, p2;
    Point e11 = nodes[e1.first];
    Point e12 = nodes[e1.second];
    Point e21 = nodes[e2.first];
    Point e22 = nodes[e2.second];

    for (ID i = 0; i < 3; i++) {
      p1[i] = (1-con.a1) * e11[i] + con.a1*e12[i];
      p2[i] = (1-con.a2) * e21[i] + con.a2*e22[i];
    }

    ID p1id = pcloud.pts.size();
    pcloud.pts.push_back(p1);
    ID p2id = pcloud.pts.size();
    pcloud.pts.push_back(p2);

    fiber_connections[con.f1].push_back({p1id, con.a1});
    fiber_connections[con.f2].push_back({p2id, con.a2});
  }

  // next all fiber endpoints
  for (const Point& p : nodes)
    pcloud.pts.push_back(p);

  const ID dim = 3, maxleaf = 10;
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
  std::vector<ID> is_in_vertices(pcloud.pts.size(), 0);
  std::vector<ID> vertices_idx(pcloud.pts.size(), (ID)-1);

  std::vector<geobin::Edge> edges;

  using Neighbor = nanoflann::ResultItem<uint32_t, Real>; // Neighbor = (id,distance)
  std::vector<Neighbor> neighbors;
  neighbors.reserve(1000); // reserve more than enough

  auto find_merged_id = [&kdtree, &is_in_vertices, &neighbors, &vertices_idx](const Point& p, Real r) {
    // NOTE: nanoflann works with squared L2 distances
    Real d = 2*r;
    const ID num_neighbors = kdtree.radiusSearch(p.data(), r*r, neighbors); // only consider neighbors closer than r
    ID index = (ID)-1;
    assert(num_neighbors != 0); // point_a should be found
    for (ID i = 0; i < num_neighbors; i++) {
      if (0 == is_in_vertices[neighbors[i].first]) // only consider neighbors already in vertices
        continue;
      if (neighbors[i].second < d) {               // of those, pick the closest
        index = vertices_idx[neighbors[i].first];
        d = neighbors[i].second;
      }
    }
    return index;
  };

  auto find_merged_id_and_insert = [&kdtree, &neighbors, &vertices, &is_in_vertices, &vertices_idx, &find_merged_id](ID id, const Point& p, Real r) {
    ID index = find_merged_id(p, r);
    if (index == (ID)-1) {
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

  for (ID cid = 1; cid < connections.size(); cid++) {
    const Point& point_a = pcloud.pts[2*cid];
    const Point& point_b = pcloud.pts[2*cid+1];

    Real r = 1e-10;

    // index_a = np.argmin(np.linalg.norm(point_a - vertices, axis=1))
    // if np.linalg.norm(point_a - vertices[index_a]) > 1e-10:
    //   index_a = len(vertices)
    //   vertices = np.vstack((vertices, point_a))

    const ID index_a = find_merged_id_and_insert(2*cid, point_a, r);

    // index_b = np.argmin(np.linalg.norm(point_b - vertices, axis=1))
    // if np.linalg.norm(point_b - vertices[index_b]) > 1e-10:
    //   index_b = len(vertices)
    //   vertices = np.vstack((vertices, point_b))

    const ID index_b = find_merged_id_and_insert(2*cid+1, point_b, r);

    if (index_a != index_b)
      edges.push_back({index_a, index_b});
  }

  logi("  from fibers");

  // NOTE: the python code below does not respect the self loops filtered out above
  //   edges_prop = np.vstack((connectionsProp, np.array(edges_prop)))
  // but we replicate this behaviour.
  std::vector<Prop>& edge_props = connection_props;

  for (ID fid = 0; fid < fibers.size(); fid++) {
    // helper = act_fibers[act_fibers[:,0] == index, 1]
    std::vector<ConnectionPoint>& helper = fiber_connections[fid];
    if (helper.empty())
      continue;

    const Point& point_a = nodes[fibers[fid].first];
    const Point& point_b = nodes[fibers[fid].second];
    const ID off = 2*connections.size();

    const ID endpoint_idx_a = find_merged_id_and_insert(off+fibers[fid].first, point_a, 1e-15);
    const ID endpoint_idx_b = find_merged_id_and_insert(off+fibers[fid].second, point_b, 1e-15);

    std::sort(helper.begin(), helper.end(), [](ConnectionPoint p1, ConnectionPoint p2){ return p1.a < p2.a; });

    // helper = [0.] + helper + [1.]
    helper.insert(helper.begin(), {.nodeid = endpoint_idx_a, .a = 0});
    helper.insert(helper.end(),   {.nodeid = endpoint_idx_b, .a = 1});

    for (ID k = 0; k < helper.size()-1; k++) {
      Point point_ab = interpolate(point_a, point_b, helper[k].a);
      Point point_ba = interpolate(point_a, point_b, helper[k+1].a);

      ID index_a = find_merged_id(point_ab, 1e-10);
      if (index_a == (ID)-1) {
        std::println(stderr, "ERROR: expected connection point to be found");
        return 1;
      }
      ID index_b = find_merged_id(point_ba, 1e-10);
      if (index_b == (ID)-1) {
        std::println(stderr, "ERROR: expected connection point to be found");
        return 1;
      }

      if (index_a != index_b) {
        edges.push_back({index_a, index_b});
        edge_props.push_back(fiber_props[fid]);
      }
    }
  }

  logi("output");
  logi("  vertices.size = {}", vertices.size());
  logi("  edges.size    = {}", edges.size());
  logi("  txt           = {}", txt);

  if (fs::is_directory(output_path))
    output_path = std::format("{}/fiber_network_{}", output_path, edges.size());

  Graph graph = { .edges = edges, .vertices = vertices, .edge_props = edge_props, .types = {}, .node_types = {}, .xadj = {}, .adjncy = {}};
  compute_types(graph);

  if (txt) {
    serialize_txt(output_path.c_str(), graph);
  } else {
    serialize_bin(output_path.c_str(), graph);
    if (testbin) {
      auto graph2 = deserialize_bin(output_path.c_str());
      serialize_txt(std::format("{}_test", output_path).c_str(), graph2);
    }
  }

}
