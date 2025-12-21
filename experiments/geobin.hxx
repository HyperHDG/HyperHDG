#ifndef GEOBIN_HXX
#define GEOBIN_HXX

#include <cstdint>
#include <vector>
#include <array>
#include <cstring>

namespace geobin {

using u64 = uint64_t;
using u32 = uint32_t;
using u8 = uint8_t;
using Real = double;
using ID = u32;

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

struct Graph {
  std::vector<Edge> edges;
  std::vector<Point> vertices;
  std::vector<Prop> edge_props;
  std::vector<Edge> types;
  std::vector<ID> node_types;
  std::vector<ID> xadj;
  std::vector<ID> adjncy;
};

struct DataTable {
  char name[8];
  u64 offset;     // from file start
  u64 size;       // in bytes
};

struct GeoBinHeader {
  char magic[8]; // 'GEOBIN1\0'
  u64 space_dim;
  u64 hyperedge_dim;
  u64 n_points;
  u64 n_hypernodes;
  u64 n_hyperedges;
  DataTable tables[5];
};

struct DomainsHeader {
  char magic[8]; // 'DOMAIN1\0'
  u64 idsize; // == sizeof(ID)
  u64 n_domains;
  DataTable tables[2];
};

std::vector<Point> read_nodes(const char* path);

std::vector<Edge> read_fibers(const char* path);

std::vector<Connection> read_connections(const char* path);

std::vector<Prop> read_props(const char* path);

ID compute_vertex_type(const Point& vertex, const Point& max_p, const Point& min_p);
Graph& compute_types(Graph& graph);

void serialize_txt(const char* output_path, const Graph& graph);

void serialize_bin(const char* output_path, const Graph& graph);

Graph deserialize_bin(const char* input_path);

void serialize_domains(const char* path, const std::vector<std::vector<ID>>& domains);

void serialize_graph_partition_vtu(
  const geobin::Graph& graph,
  const std::vector<geobin::ID>& partition,
  const char* file_path
);

}

#endif // GEOBIN_HXX
