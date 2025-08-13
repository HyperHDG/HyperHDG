#ifndef GEOBIN_HXX
#define GEOBIN_HXX

#include <cstdint>
#include <vector>
#include <array>
#include <cstring>
#include <print>
#include <format>

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

struct Graph {
  std::vector<Edge> edges;
  std::vector<Point> vertices;
  std::vector<Prop> edge_props;
  std::vector<Edge> types;
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

Graph& compute_types(Graph& graph);

void serialize_txt(const char* output_path, const Graph& graph);

void serialize_bin(const char* output_path, const Graph& graph);

Graph deserialize_bin(const char* input_path);

void serialize_domains(const char* path, const std::vector<std::vector<ID>>& domains);

}

#endif // GEOBIN_HXX
