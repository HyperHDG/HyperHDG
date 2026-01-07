#include <vector>
#include <cstdint>
#include <utility> // for pair
#include <fstream>
#include <cstdio>
#include <cstring>
#include <algorithm>

#include <nanoflann.hpp>
#include <petsc.h>
#include <petscviewerhdf5.h>

using Real = double;
using ID = unsigned int;

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
    fprintf(stderr, "error: couldn't open '%s'\n", path);
    return {};
  }
  for (std::string line; std::getline(nodes_file, line); ) {
    ID id; Real x, y, z;
    if (4 == sscanf(line.c_str(), "%d,%lf,%lf,%lf", &id, &x, &y, &z)) {
      nodes.push_back({x,y,z});
    } else {
      if (strcmp("Id,x,y,z", line.c_str()) == 0)
        continue;
      else {
        fprintf(stderr, "error: %s: couldn't parse line '%s'\n", path, line.c_str());
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
    fprintf(stderr, "error: couldn't open '%s'\n", path);
    return {};
  }
  for (std::string line; std::getline(fibers_file, line); ) {
    ID f; ID u, v;
    if (3 == sscanf(line.c_str(), "%d,%d,%d", &f, &u, &v)) {
      fibers.push_back({u,v});
    } else {
      if (strcmp("Id,node1,node2", line.c_str()) == 0)
        continue;
      else {
        fprintf(stderr, "error: %s: couldn't parse line '%s'\n", path, line.c_str());
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
    fprintf(stderr, "error: couldn't open '%s'\n", path);
    return {};
  }
  for (std::string line; std::getline(connections_file, line); ) {
    ID c,f1,f2; Real a1,a2;
    if (5 == sscanf(line.c_str(), "%d,%d,%d,%lf,%lf", &c, &f1, &f2, &a1, &a2)) {
      connections.push_back({f1,f2,a1,a2});
    } else {
      if (strcmp("Id,fiber1,fiber2,a1,a2", line.c_str()) == 0)
        continue;
      else {
        fprintf(stderr, "error: %s: couldn't parse line '%s'\n", path, line.c_str());
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
    fprintf(stderr, "error: couldn't open '%s'\n", path);
    return {};
  }
  for (std::string line; std::getline(fiber_props_file, line); ) {
    ID f;
    Prop prop; // 6 structural + 2*3 normal
    int total_chars_read = 0, chars_read = 0;
    const char* cline = line.c_str();

    if (1 != sscanf(cline, "%d%n,", &f, &chars_read)) {
      if (0 == strcmp(cline, "Id,EA,kG_1A,kG_2A,G_xI_x,E_1I_1,E_2I_2,n_11,n_12,n_13,n_21,n_22,n_23"))
        continue;
      fprintf(stderr, "error: %s: sscanf format '%%d' invalid for '%s'\n", path, cline);
      return {};
    }

    total_chars_read += chars_read+1;

    for (ID i = 0; i < 12; i++) {
      if (1 != sscanf(cline+total_chars_read, "%lf%n", &prop[i], &chars_read)) {
        fprintf(stderr, "error: %s: sscanf format '%%lf' invalid for '%s'\n", path, cline+total_chars_read);
        return {};
      }
      total_chars_read += chars_read+1; // skip ','
    }

    props.push_back(prop);
  }

  return props;
}

ID compute_vertex_type(const Point& vertex, const Point& max_p, const Point& min_p) {
    if (vertex[0] - min_p[0] < 1e-6 * (max_p[0] - min_p[0]) || max_p[0] - vertex[0] < 1e-6 * (max_p[0] - min_p[0]) ||
        vertex[1] - min_p[1] < 1e-6 * (max_p[1] - min_p[1]) || max_p[1] - vertex[1] < 1e-6 * (max_p[1] - min_p[1]))
      return 1;
    else
      return 0;
}

void compute_types(const std::vector<Point>& vertices, const std::vector<Edge>& edges, std::vector<ID>& node_types, std::vector<Edge>& types) {
  node_types.resize(vertices.size());
  types.resize(edges.size()*2);

  // Calculate the bounding box (min/max x, y, z) for all vertices
  Point min_p = {1e10, 1e10, 1e10};
  Point max_p = {1e-10, 1e-10, 1e-10};
  for (const Point& vertex : vertices) {
    for (size_t i = 0; i < 3; i++) {
      min_p[i] = std::min(min_p[i], vertex[i]);
      max_p[i] = std::max(max_p[i], vertex[i]);
    }
  }
  for (size_t n = 0; n < vertices.size(); n++)
    node_types[n] = compute_vertex_type(vertices[n], max_p, min_p);

  for (size_t m = 0; m < edges.size(); m++) {
    auto edge = edges[m];
    types[m] = {node_types[edge.first], node_types[edge.second]};
  }
}

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

int main(int argc, char** argv) {
  PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

  PetscBool is_set, help;
  char input_folder[PATH_MAX] = ".";
  char output_path[PATH_MAX] = "network.h5";
  char buf[PATH_MAX];
  const size_t bufsz = sizeof(buf);

  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "make_geo options", NULL);
  PetscCall(PetscOptionsString("-i", "input folder", NULL, input_folder, input_folder, PATH_MAX, &is_set));
  PetscCall(PetscOptionsString("-o", "output path", NULL, output_path, output_path, PATH_MAX, &is_set));
  PetscOptionsEnd();

  PetscCall(PetscOptionsGetBool(NULL, NULL, "-help", &help, &is_set));
  if (help) {
    PetscOptionsView(NULL, PETSC_VIEWER_STDOUT_WORLD);
    PetscFinalize();
    return 0;
  }

  printf("input:\n");
  printf("  folder: %s\n", input_folder);

  snprintf(buf, bufsz, "%s/nodes.csv", input_folder);
  std::vector<Point> nodes = read_nodes(buf);
  printf("  nodes: %zu\n", nodes.size());

  snprintf(buf, bufsz, "%s/fibers.csv", input_folder);
  std::vector<Edge> fibers = read_fibers(buf);
  printf("  fibers: %zu\n", fibers.size());

  snprintf(buf, bufsz, "%s/connections.csv", input_folder);
  std::vector<Connection> connections = read_connections(buf);
  printf("  connections: %zu\n", connections.size());

  snprintf(buf, bufsz, "%s/fibersProps.csv", input_folder);
  std::vector<Prop> fiber_props = read_props(buf);
  printf("  fibersProps: %zu\n", fiber_props.size());

  snprintf(buf, bufsz, "%s/connectionsProp.csv", input_folder);
  std::vector<Prop> connection_props = read_props(buf);
  printf("  connectionProps: %zu\n", connection_props.size());

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
  printf("# building kdtree\n");
  kdtree.buildIndex();

  printf("# generating vertex/edge list\n");

  std::vector<Point> vertices;
  std::vector<ID> is_in_vertices(pcloud.pts.size(), 0);
  std::vector<ID> vertices_idx(pcloud.pts.size(), (ID)-1);

  std::vector<Edge> edges;

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

  printf("#   from connections\n");

  // if len(vertices) == 0: vertices = np.vstack((point_a, point_b))
  if (pcloud.pts.size() >= 2) {
    vertices.push_back(pcloud.pts[0]);
    vertices.push_back(pcloud.pts[1]);
    vertices_idx[0] = 0;
    vertices_idx[1] = 1;
    is_in_vertices[0] = 1;
    is_in_vertices[1] = 1;
    edges.push_back({0,1});
  }

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

  printf("#   from fibers\n");

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
        fprintf(stderr, "ERROR: expected connection point to be found");
        return 1;
      }
      ID index_b = find_merged_id(point_ba, 1e-10);
      if (index_b == (ID)-1) {
        fprintf(stderr, "ERROR: expected connection point to be found");
        return 1;
      }

      if (index_a != index_b) {
        edges.push_back({index_a, index_b});
        edge_props.push_back(fiber_props[fid]);
      }
    }
  }

  snprintf(buf, bufsz, "%s.geo.h5", output_path);

  printf("output:\n");
  printf("  vertices_size: %zu\n", vertices.size());
  printf("  edges_size: %zu\n", edges.size());
  printf("  path: %s\n", buf);

  std::vector<ID> node_types;
  std::vector<Edge> types;
  compute_types(vertices, edges, node_types, types);

  IS is_edges, types_faces, types_points;
  Vec points, properties;
  PetscViewer viewer;

  PetscCall(ISCreateGeneral(PETSC_COMM_SELF, edges.size()*2, (PetscInt*)edges.data(), PETSC_USE_POINTER, &is_edges));
  PetscCall(ISSetBlockSize(is_edges, 2));
  PetscCall(PetscObjectSetName((PetscObject)is_edges, "edges"));

  PetscCall(ISCreateGeneral(PETSC_COMM_SELF, node_types.size(), (PetscInt*)node_types.data(), PETSC_USE_POINTER, &types_points));
  PetscCall(PetscObjectSetName((PetscObject)types_points, "types_points"));

  PetscCall(ISCreateGeneral(PETSC_COMM_SELF, types.size(), (PetscInt*)types.data(), PETSC_USE_POINTER, &types_faces));
  PetscCall(ISSetBlockSize(types_faces, 2));
  PetscCall(PetscObjectSetName((PetscObject)types_faces, "types_faces"));

  PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 3, vertices.size()*3, (PetscReal*)vertices.data(), &points));
  PetscCall(PetscObjectSetName((PetscObject)points, "points"));
  PetscCall(VecCreateSeqWithArray(PETSC_COMM_SELF, 12, edge_props.size()*12, (PetscReal*)edge_props.data(), &properties));
  PetscCall(PetscObjectSetName((PetscObject)properties, "properties"));

  PetscCall(PetscViewerHDF5Open(PETSC_COMM_SELF, buf, FILE_MODE_WRITE, &viewer));
  PetscCall(PetscViewerHDF5SetCompress(viewer, PETSC_TRUE));

  PetscCall(PetscViewerHDF5PushGroup(viewer, "/domain"));
  PetscCall(ISView(is_edges, viewer));
  PetscCall(ISView(types_points, viewer));
  PetscCall(ISView(types_faces, viewer));
  PetscCall(VecView(points, viewer));
  PetscCall(VecView(properties, viewer));

  PetscCall(VecDestroy(&points));
  PetscCall(VecDestroy(&properties));
  PetscCall(ISDestroy(&is_edges));
  PetscCall(ISDestroy(&types_points));
  PetscCall(ISDestroy(&types_faces));
  PetscCall(PetscViewerDestroy(&viewer));

  PetscCall(PetscFinalize());
}
