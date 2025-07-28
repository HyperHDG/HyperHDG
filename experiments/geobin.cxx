#include "geobin.hxx"

#include <bxzstr.hpp>
#include <fstream>
#include <fmtlog/fmtlog.h>

namespace geobin {

std::vector<Point> read_nodes(const char* path) {
  std::vector<Point> nodes;
  std::fstream nodes_file(path);
  if (!nodes_file) {
    std::println("error: couldn't open {}", path);
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
    if (3 == sscanf(line.c_str(), "%d,%d,%d", &f, &u, &v)) {
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
    if (5 == sscanf(line.c_str(), "%d,%d,%d,%lf,%lf", &c, &f1, &f2, &a1, &a2)) {
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

    if (1 != sscanf(cline, "%d%n,", &f, &chars_read)) {
      if (0 == strcmp(cline, "Id,EA,kG_1A,kG_2A,G_xI_x,E_1I_1,E_2I_2,n_11,n_12,n_13,n_21,n_22,n_23"))
        continue;
      std::println(stderr, "error: {}: sscanf format '%d' invalid for '{}'", path, cline);
      return {};
    }

    total_chars_read += chars_read+1;

    for (ID i = 0; i < 12; i++) {
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

void GraphEdgeList::to_access(graph_access& graph_acc) {
  // is directed -> 2*edges
  graph_acc.start_construction(vertices.size(), 2*edges.size());

  std::vector<std::vector<NodeID>> adjacency(vertices.size());
  for (const geobin::Edge& edge : edges) {
    adjacency[edge.first].push_back(edge.second);
    adjacency[edge.second].push_back(edge.first);
  }

  for (NodeID n = 0; n < vertices.size(); n++) {
    NodeID nn = graph_acc.new_node();
    graph_acc.setNodeWeight(nn, 1);
    for (const NodeID neighbor : adjacency[n]) {
      EdgeID e = graph_acc.new_edge(nn, neighbor);
      graph_acc.setEdgeWeight(e, 1);
    }
  }

  graph_acc.finish_construction();
}

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
    u32 left = 0, right = 0;
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
    for (ID i = 1; i < 12; i++)
      std::print(gfile, "  {:.18e}", prop[i]);
    std::print(gfile, "\n");
  }

  std::ofstream pfile(std::format("{}.pts", output_path));
  for (const Point& vertex : vertices)
    std::print(pfile, "{:.18e} {:.18e} {:.18e}\n", vertex[0], vertex[1], vertex[2]);
}

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
    .size = m * 2 * sizeof(ID),
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

  bxz::ofstream file(std::format("{}.geo.bin.zstd", output_path), bxz::zstd);
  file.exceptions(std::ofstream::badbit | std::ofstream::failbit);
  file.write((char*)&header, sizeof(header));
  file.write((char*)&graph.vertices[0], points.size);
  file.write((char*)&graph.edges[0], hypernodes_of_hyperedges.size);
  file.write((char*)&graph.types[0], types_of_hyperfaces.size);
  file.write((char*)&graph.edges[0], points_of_hyperedges.size);
  file.write((char*)&graph.edge_props[0], hyperedge_properties.size);
}

GraphEdgeList deserialize_bin(const char* input_path) {
  GraphEdgeList graph;

  bxz::ifstream file(input_path, std::ios::binary);
  file.exceptions(std::ifstream::badbit | std::ifstream::failbit);
  std::array<char, sizeof(GeoBinHeader)> header_buf;
  file.read(header_buf.data(), sizeof(GeoBinHeader));
  GeoBinHeader* header = (GeoBinHeader*)header_buf.data();
  assert(0 == std::strcmp(header->magic, "GEOBIN1")); // extra trailing \0
  assert(3 == header->space_dim);
  assert(1 == header->hyperedge_dim);

  graph.vertices.resize(header->n_points);
  graph.edges.resize(header->n_hyperedges);
  graph.edge_props.resize(header->n_hyperedges);
  graph.types.resize(header->n_hyperedges);

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
