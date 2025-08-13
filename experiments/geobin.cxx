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

ID compute_vertex_type(const Point& vertex, const Point& max_p, const Point& min_p) {
    if (vertex[0] - min_p[0] < 1e-6 * (max_p[0] - min_p[0]) || max_p[0] - vertex[0] < 1e-6 * (max_p[0] - min_p[0]) ||
        vertex[1] - min_p[1] < 1e-6 * (max_p[1] - min_p[1]) || max_p[1] - vertex[1] < 1e-6 * (max_p[1] - min_p[1]))
      return 1;
    else
      return 0;
}

Graph& compute_types(Graph& graph) {
  graph.node_types.resize(graph.vertices.size());
  graph.types.resize(graph.edges.size());

  // Calculate the bounding box (min/max x, y, z) for all vertices
  Point min_p = {1e10};
  Point max_p = {1e-10};
  for (const Point& vertex : graph.vertices) {
    for (geobin::u64 i = 0; i < 3; i++) {
      min_p[i] = std::min(min_p[i], vertex[i]);
      max_p[i] = std::max(max_p[i], vertex[i]);
    }
  }
  for (geobin::ID n = 0; n < graph.vertices.size(); n++)
    graph.node_types[n] = compute_vertex_type(graph.vertices[n], max_p, min_p);

  for (const Edge& edge : graph.edges)
    graph.types.push_back({graph.node_types[edge.first], graph.node_types[edge.second]});

  return graph;
}

void serialize_txt(const char* output_path, const Graph& graph) {
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

void serialize_bin(const char* output_path, const Graph& graph) {
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

Graph deserialize_bin(const char* input_path) {
  Graph graph;

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

  geobin::ID nverts = graph.vertices.size();
  geobin::ID nedges = graph.edges.size();
  std::vector<std::vector<geobin::ID>> adjacency(nverts+1);
  for (const geobin::Edge& edge : graph.edges) {
    adjacency[edge.first].push_back(edge.second);
    adjacency[edge.second].push_back(edge.first);
  }
  graph.xadj.resize(nverts+1, 0);
  graph.adjncy.resize(nedges*2, 0); // *2 for directed repr
  geobin::ID edges_so_far = 0;
  for (geobin::ID n = 0; n < nverts+1; n++) {
    graph.xadj[n] = edges_so_far;
    for (geobin::ID nid = 0; nid < adjacency[n].size(); nid++) {
      assert(nid + edges_so_far < nedges*2);
      graph.adjncy[nid + edges_so_far] = adjacency[n][nid];
    }
    edges_so_far += adjacency[n].size();
  }
  assert(edges_so_far == nedges*2);

  compute_types(graph);

  return graph;
}

void serialize_domains(const char* path, const std::vector<std::vector<geobin::ID>>& domains) {
  bxz::ofstream file(std::format("{}.dom.zstd", path), bxz::zstd);

  DataTable ioffset_table = {
    .name = "IOFFSET",
    .offset = sizeof(DomainsHeader),
    .size = (1+domains.size())*sizeof(geobin::ID),
  };

  DataTable domains_table = {
    .name = "DOMAINS",
    .offset = ioffset_table.offset + ioffset_table.size,
    .size = 0,
  };
  for (u64 p = 0; p < domains.size(); p++)
    domains_table.size += domains[p].size() * sizeof(geobin::ID);

  DomainsHeader hdr = {
    .magic = "DOMAIN1",
    .idsize = sizeof(geobin::ID),
    .n_domains = domains.size(),
    .tables = {
      ioffset_table,
      domains_table,
    }
  };
  file.write((char*)&hdr, sizeof(hdr));

  geobin::ID ioff = 0;
  std::vector<geobin::ID> ioffsets(1+domains.size());
  for (u64 p = 0; p < domains.size(); p++) {
    ioffsets[p] = ioff;
    ioff += domains[p].size();
  }
  ioffsets[domains.size()] = ioff;
  file.write((char*)&ioffsets[0], ioffsets.size()*sizeof(geobin::ID));

  for (geobin::ID p = 0; p < domains.size(); p++)
    file.write((char*)&domains[p][0], domains[p].size()*sizeof(geobin::ID));
}

void serialize_graph_partition_vtu(
  const geobin::Graph& graph,
  const std::vector<geobin::ID>& partition,
  const char* file_path
) {
  std::ofstream file(std::format("{}.vtu", file_path));
  if (!file.is_open()) {
    loge("could not opt file {}.vtu", file_path);
    return;
  }

  std::print(file, "<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  std::print(file, "  <UnstructuredGrid>\n");
  std::print(file, "    <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n", graph.vertices.size(), graph.edges.size());
  std::print(file, "      <Points>\n");
  std::print(file, "        <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for (const auto& point : graph.vertices)
      std::print(file, "          {} {} {}\n", point[0], point[1], point[2]);
  std::print(file, "        </DataArray>\n");
  std::print(file, "      </Points>\n");
  std::print(file, "      <Cells>\n");
  std::print(file, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n");
  std::print(file, "          ");
  for (const auto& edge : graph.edges)
      std::print(file, "{} {} ", edge.first, edge.second);
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  // offsets: the cumulative sum of the number of points in each cell.
  // for VTK_LINE (2 points per cell), this will be 2, 4, 6, ...
  std::print(file, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n");
  std::print(file, "          ");
  for (size_t i = 0; i < graph.edges.size(); ++i)
      std::print(file, "{} ", (i + 1) * 2);
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  std::print(file, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  std::print(file, "          ");
  for (size_t i = 0; i < graph.edges.size(); ++i)
      std::print(file, "3 "); // VTK_LINE type
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  std::print(file, "      </Cells>\n");
  std::print(file, "      <PointData>\n");
  std::print(file, "        <DataArray type=\"Int32\" Name=\"PartitionID\" format=\"ascii\">\n");
  std::print(file, "          ");
  for (geobin::ID id : partition)
      std::print(file, "{} ", id);
  std::print(file, "\n");
  std::print(file, "        </DataArray>\n");
  std::print(file, "      </PointData>\n");
  std::print(file, "    </Piece>\n");
  std::print(file, "  </UnstructuredGrid>\n");
  std::print(file, "</VTKFile>\n");
}

}
