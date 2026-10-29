#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/hy_assert.hxx>

#include <hdf5.h>

#include <algorithm>
#include <array>
#include <cstdint>
#include <string>
#include <type_traits>
#include <vector>

// Forward declaration: the struct is defined in read_domain.hxx, which includes this header.
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          template <typename...>
          typename vectorT,
          typename pointT,
          typename hyEdge_index_t,
          typename hyNode_index_t,
          typename pt_index_t>
struct DomainInfo;

/*!*************************************************************************************************
 * \brief   The .geo.h5 domain file format and the HDF5 calls to read and write it.
 *
 * A .geo.h5 file holds the same information as a .geo file (cf. read_domain_geo) in the group
 * \c /domain, as plain two-dimensional datasets that NumPy/h5py read directly:
 *
 *     /domain/points        (n_points,     space_dim)      float   coordinates
 *     /domain/edges         (n_hyEdges,    2^hyEdge_dim)   int     point indices per hyperedge
 *     /domain/hypernodes    (n_hyEdges,    2 hyEdge_dim)   int     hypernode indices per hyperedge
 *     /domain/types_faces   (n_hyEdges,    2 hyEdge_dim)   int     face types (boundary conditions)
 *     /domain/types_points  (n_points,     1)              int     derived (hyEdge_dim 1): max
 *                                                                  type of the faces at a point
 *     /domain/properties    (n_hyEdges,    n_properties)   float   optional per-hyperedge data
 *
 * with the attributes \c space_dim, \c hyperedge_dim and \c n_hypernodes on \c /domain (integers;
 * further attributes are ignored). \c hypernodes is optional: without it, hypernodes are numbered
 * like points (one hypernode per point), which is the case for networks (hyEdge_dim 1) and the
 * convention of the Python network generators. Files written by \c write_domain_hdf5 for
 * hyEdge_dim 1 contain no \c hypernodes dataset if the numbering agrees.
 *
 * \c tools/geo2h5 converts a .geo file into this format. The group \c /VTKHDF (optional, written by
 *default) is a VTKHDF unstructured grid view of the same data for ParaView: \c Points and the
 *point/cell data are links into \c /domain, and the \c Connectivity of a network is a virtual
 *dataset over \c /domain/edges, so the view costs nothing but a few small arrays. ParaView chooses
 *its reader by file extension, so open a .geo.h5 file with File > Open and pick "VTKHDF Reader"
 *when asked.
 *
 * Reading is one \c H5Dread per dataset into the (contiguous) containers of \c DomainInfo;
 * types are converted by HDF5, so int32 files load into unsigned int indices and float32 points
 * into double. The layout is row-major per hyperedge, which is what a later per-rank hyperslab
 * read needs.
 *
 * Requires HyperHDG to be configured with HDF5 (CMake option \c HYPERHDG_HDF5, which defines
 * \c HYPERHDG_HAVE_HDF5). The functions throw \c std::runtime_error on malformed files.
 *
 * \authors   Joseph Holten, Karlsruhe Institute of Technology, 2026.
 **************************************************************************************************/
namespace HDF5Domain
{
/*!*************************************************************************************************
 * \brief   HDF5 identifier that is closed when it goes out of scope.
 **************************************************************************************************/
class Id
{
 public:
  Id(const hid_t id, herr_t (*close)(hid_t), const std::string& what) : id_(id), close_(close)
  {
    hy_check(id_ >= 0, "HDF5: " << what << " failed.");
  }
  ~Id()
  {
    if (id_ >= 0)
      close_(id_);
  }
  Id(const Id&) = delete;
  Id& operator=(const Id&) = delete;
  operator hid_t() const { return id_; }

 private:
  hid_t id_;
  herr_t (*close_)(hid_t);
};

/*!*************************************************************************************************
 * \brief   The HDF5 memory type of a C++ arithmetic type.
 **************************************************************************************************/
template <typename T>
hid_t native_type()
{
  if constexpr (std::is_same_v<T, double>)
    return H5T_NATIVE_DOUBLE;
  else if constexpr (std::is_same_v<T, float>)
    return H5T_NATIVE_FLOAT;
  else if constexpr (std::is_same_v<T, long double>)
    return H5T_NATIVE_LDOUBLE;
  else if constexpr (std::is_same_v<T, int>)
    return H5T_NATIVE_INT;
  else if constexpr (std::is_same_v<T, unsigned int>)
    return H5T_NATIVE_UINT;
  else if constexpr (std::is_same_v<T, long>)
    return H5T_NATIVE_LONG;
  else if constexpr (std::is_same_v<T, unsigned long>)
    return H5T_NATIVE_ULONG;
  else if constexpr (std::is_same_v<T, long long>)
    return H5T_NATIVE_LLONG;
  else if constexpr (std::is_same_v<T, unsigned long long>)
    return H5T_NATIVE_ULLONG;
  else if constexpr (std::is_same_v<T, short>)
    return H5T_NATIVE_SHORT;
  else if constexpr (std::is_same_v<T, unsigned short>)
    return H5T_NATIVE_USHORT;
  else if constexpr (std::is_same_v<T, signed char>)
    return H5T_NATIVE_SCHAR;
  else if constexpr (std::is_same_v<T, unsigned char>)
    return H5T_NATIVE_UCHAR;
  else
    static_assert(!std::is_same_v<T, T>, "No HDF5 type for this C++ type.");
}

/*!*************************************************************************************************
 * \brief   True if \c name exists below \c loc.
 **************************************************************************************************/
inline bool exists(const hid_t loc, const char* name)
{
  return H5Lexists(loc, name, H5P_DEFAULT) > 0;
}

/*!*************************************************************************************************
 * \brief   Rows and columns of a one- or two-dimensional dataset (a 1d dataset has one column).
 **************************************************************************************************/
inline std::array<hsize_t, 2> shape(const hid_t loc, const char* name)
{
  Id dset(H5Dopen2(loc, name, H5P_DEFAULT), H5Dclose, std::string("open dataset ") + name);
  Id space(H5Dget_space(dset), H5Sclose, "get dataspace");
  const int rank = H5Sget_simple_extent_ndims(space);
  hy_check(rank == 1 || rank == 2,
           "HDF5: dataset " << name << " must be 1d or 2d, but has rank " << rank << ".");
  hsize_t dims[2] = {0, 1};
  H5Sget_simple_extent_dims(space, dims, nullptr);
  return {dims[0], rank == 2 ? dims[1] : 1};
}

/*!*************************************************************************************************
 * \brief   Read a whole dataset of \c rows x \c cols entries into \c buffer (row-major).
 **************************************************************************************************/
template <typename T>
void read(const hid_t loc, const char* name, const hsize_t rows, const hsize_t cols, T* buffer)
{
  const std::array<hsize_t, 2> dims = shape(loc, name);
  hy_check(dims[0] == rows && dims[1] == cols,
           "HDF5: dataset " << name << " has shape (" << dims[0] << ", " << dims[1]
                            << "), expected (" << rows << ", " << cols << ").");
  if (rows * cols == 0)
    return;
  Id dset(H5Dopen2(loc, name, H5P_DEFAULT), H5Dclose, std::string("open dataset ") + name);
  hy_check(H5Dread(dset, native_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer) >= 0,
           "HDF5: reading dataset " << name << " failed.");
}

/*!*************************************************************************************************
 * \brief   Read an integer attribute of \c loc if it exists; return whether it does.
 **************************************************************************************************/
inline bool read_attribute(const hid_t loc, const char* name, long long& value)
{
  if (H5Aexists(loc, name) <= 0)
    return false;
  Id attr(H5Aopen(loc, name, H5P_DEFAULT), H5Aclose, std::string("open attribute ") + name);
  hy_check(H5Aread(attr, H5T_NATIVE_LLONG, &value) >= 0,
           "HDF5: reading attribute " << name << " failed.");
  return true;
}

/*!*************************************************************************************************
 * \brief   Write an integer attribute to \c loc.
 **************************************************************************************************/
inline void write_attribute(const hid_t loc, const char* name, const long long value)
{
  Id space(H5Screate(H5S_SCALAR), H5Sclose, "create scalar dataspace");
  Id attr(H5Acreate2(loc, name, H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose,
          std::string("create attribute ") + name);
  hy_check(H5Awrite(attr, H5T_NATIVE_LLONG, &value) >= 0,
           "HDF5: writing attribute " << name << " failed.");
}

/*!*************************************************************************************************
 * \brief   Write a fixed-length string attribute to \c loc (as NumPy's bytes_ does).
 **************************************************************************************************/
inline void write_attribute(const hid_t loc, const char* name, const std::string& value)
{
  Id type(H5Tcopy(H5T_C_S1), H5Tclose, "copy string type");
  H5Tset_size(type, value.size());
  H5Tset_strpad(type, H5T_STR_NULLPAD);
  Id space(H5Screate(H5S_SCALAR), H5Sclose, "create scalar dataspace");
  Id attr(H5Acreate2(loc, name, type, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose,
          std::string("create attribute ") + name);
  hy_check(H5Awrite(attr, type, value.data()) >= 0,
           "HDF5: writing attribute " << name << " failed.");
}

/*!*************************************************************************************************
 * \brief   Write a dataset of \c rows x \c cols entries from \c data (row-major).
 *
 * \c cols = 0 writes a one-dimensional dataset of \c rows entries. Datasets with at least one
 * entry are chunked and gzip compressed (level 1, the level at which network files compress
 * well without slowing down the write).
 **************************************************************************************************/
template <typename T>
void write(const hid_t loc,
           const char* name,
           const hsize_t rows,
           const hsize_t cols,
           const T* data,
           const bool compress = true)
{
  const int rank = cols == 0 ? 1 : 2;
  const hsize_t dims[2] = {rows, cols};
  Id space(H5Screate_simple(rank, dims, nullptr), H5Sclose, "create dataspace");
  Id dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, "create dataset property list");
  if (compress && rows > 0)
  {
    const hsize_t chunk[2] = {std::min<hsize_t>(rows, hsize_t(1) << 16), cols};
    H5Pset_chunk(dcpl, rank, chunk);
    H5Pset_deflate(dcpl, 1);
  }
  Id dset(H5Dcreate2(loc, name, native_type<T>(), space, H5P_DEFAULT, dcpl, H5P_DEFAULT), H5Dclose,
          std::string("create dataset ") + name);
  if (rows * (cols == 0 ? 1 : cols) > 0)
    hy_check(H5Dwrite(dset, native_type<T>(), H5S_ALL, H5S_ALL, H5P_DEFAULT, data) >= 0,
             "HDF5: writing dataset " << name << " failed.");
}

/*!*************************************************************************************************
 * \brief   Write the VTKHDF unstructured-grid view of \c /domain into \c file.
 *
 * Hyperedges become VTK lines, quads or hexahedra (hyEdge_dim 1, 2, 3). Points are linked for
 * three space dimensions and zero-padded otherwise (VTK wants three coordinates). For networks
 * the connectivity is a virtual dataset over \c /domain/edges; for quads and hexahedra the
 * lexicographic vertex order of a hypercube is permuted into VTK's order.
 **************************************************************************************************/
template <unsigned int hyEdge_dim, unsigned int space_dim, typename PointsT, typename EdgesT>
void write_vtkhdf_view(const hid_t file,
                       const PointsT& points,
                       const EdgesT& edges,
                       const bool has_types_points,
                       const bool has_properties)
{
  static_assert(hyEdge_dim >= 1 && hyEdge_dim <= 3 && space_dim <= 3,
                "The VTKHDF view exists for lines, quads and hexahedra in up to 3 dimensions.");
  constexpr unsigned int n_verts = 1 << hyEdge_dim;
  constexpr unsigned char vtk_type = hyEdge_dim == 1 ? 3 : hyEdge_dim == 2 ? 9 : 12;
  const hsize_t n_points = points.size(), n_cells = edges.size(), n_conn = n_cells * n_verts;

  Id root(H5Gcreate2(file, "VTKHDF", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
          "create group VTKHDF");
  {
    const long long version[2] = {2, 0};
    const hsize_t two = 2;
    Id space(H5Screate_simple(1, &two, nullptr), H5Sclose, "create dataspace");
    Id attr(H5Acreate2(root, "Version", H5T_STD_I64LE, space, H5P_DEFAULT, H5P_DEFAULT), H5Aclose,
            "create attribute Version");
    hy_check(H5Awrite(attr, H5T_NATIVE_LLONG, version) >= 0, "HDF5: writing Version failed.");
  }
  write_attribute(root, "Type", std::string("UnstructuredGrid"));

  if constexpr (space_dim == 3)
    hy_check(H5Lcreate_soft("/domain/points", root, "Points", H5P_DEFAULT, H5P_DEFAULT) >= 0,
             "HDF5: linking Points failed.");
  else
  {
    using coord_t = typename PointsT::value_type::value_type;
    std::vector<coord_t> padded(3 * n_points, coord_t(0));
    for (hsize_t pt = 0; pt < n_points; ++pt)
      for (unsigned int dim = 0; dim < space_dim; ++dim)
        padded[3 * pt + dim] = points[pt][dim];
    write(root, "Points", n_points, 3, padded.data());
  }

  if constexpr (hyEdge_dim == 1)
  {
    using index_t = typename EdgesT::value_type::value_type;
    Id vspace(H5Screate_simple(1, &n_conn, nullptr), H5Sclose, "create virtual dataspace");
    const hsize_t src_dims[2] = {n_cells, n_verts};
    Id sspace(H5Screate_simple(2, src_dims, nullptr), H5Sclose, "create source dataspace");
    Id dcpl(H5Pcreate(H5P_DATASET_CREATE), H5Pclose, "create dataset property list");
    hy_check(H5Pset_virtual(dcpl, vspace, ".", "/domain/edges", sspace) >= 0,
             "HDF5: virtual mapping of Connectivity failed.");
    Id dset(H5Dcreate2(root, "Connectivity", native_type<index_t>(), vspace, H5P_DEFAULT, dcpl,
                       H5P_DEFAULT),
            H5Dclose, "create dataset Connectivity");
  }
  else
  {
    // lexicographic (x fastest) hypercube vertices -> VTK_QUAD / VTK_HEXAHEDRON order
    constexpr std::array<unsigned int, 8> order = {0, 1, 3, 2, 4, 5, 7, 6};
    std::vector<long long> conn(n_conn);
    for (hsize_t cell = 0; cell < n_cells; ++cell)
      for (unsigned int v = 0; v < n_verts; ++v)
        conn[cell * n_verts + v] = edges[cell][order[v]];
    write(root, "Connectivity", n_conn, 0, conn.data());
  }

  std::vector<long long> offsets(n_cells + 1);
  for (hsize_t cell = 0; cell <= n_cells; ++cell)
    offsets[cell] = cell * n_verts;
  write(root, "Offsets", n_cells + 1, 0, offsets.data());
  const std::vector<unsigned char> types(n_cells, vtk_type);
  write(root, "Types", n_cells, 0, types.data());
  const long long counts[3] = {(long long)n_points, (long long)n_cells, (long long)n_conn};
  write(root, "NumberOfPoints", 1, 0, counts, false);
  write(root, "NumberOfCells", 1, 0, counts + 1, false);
  write(root, "NumberOfConnectivityIds", 1, 0, counts + 2, false);

  if (has_types_points)
  {
    Id pd(H5Gcreate2(root, "PointData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
          "create group PointData");
    H5Lcreate_soft("/domain/types_points", pd, "types_points", H5P_DEFAULT, H5P_DEFAULT);
  }
  if (has_properties)
  {
    Id cd(H5Gcreate2(root, "CellData", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
          "create group CellData");
    H5Lcreate_soft("/domain/properties", cd, "properties", H5P_DEFAULT, H5P_DEFAULT);
  }
}

}  // end of namespace HDF5Domain

/*!*************************************************************************************************
 * \brief   Whether \c filename names a .geo.h5 domain file (an HDF5 file with a group /domain).
 **************************************************************************************************/
inline bool is_domain_hdf5(const std::string& filename)
{
  if (H5Fis_accessible(filename.c_str(), H5P_DEFAULT) <= 0)
    return false;
  const hid_t file = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  if (file < 0)
    return false;
  const bool has_domain = H5Lexists(file, "/domain", H5P_DEFAULT) > 0;
  H5Fclose(file);
  return has_domain;
}

/*!*************************************************************************************************
 * \brief   Read a .geo.h5 domain file (cf. namespace HDF5Domain for the format).
 *
 * \tparam  hyEdge_dim      The local dimension of a hyperedge.
 * \tparam  space_dim       The dimension of the surrounding space.
 * \tparam  vectorT         The typename of a large vector holding e.g. all points.
 * \tparam  pointT          The typename of a point.
 * \tparam  hyEdge_index_t  The index type for hyperedges.
 * \tparam  hyNode_index_t  The index type for hypernodes.
 * \tparam  pt_index_t      The index type for points.
 * \param   filename        Name of the .geo.h5 file.
 * \retval  domain_info     Topological and geometrical information of the hypergraph.
 **************************************************************************************************/
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          template <typename...>
          typename vectorT,
          typename pointT,
          typename hyEdge_index_t,
          typename hyNode_index_t,
          typename pt_index_t>
DomainInfo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
read_domain_hdf5(const std::string& filename)
{
  using namespace HDF5Domain;
  using coord_t = typename pointT::value_type;
  static_assert(sizeof(pointT) == space_dim * sizeof(coord_t),
                "Points must be contiguous arrays of their coordinates.");
  constexpr unsigned int n_verts = 1 << hyEdge_dim, n_faces = 2 * hyEdge_dim;

  Id file(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose,
          "open file " + filename);
  Id domain(H5Gopen2(file, "/domain", H5P_DEFAULT), H5Gclose, "open group /domain of " + filename);

  long long attr;
  if (read_attribute(domain, "space_dim", attr))
    hy_check(attr == space_dim, "space_dim in " << filename << " is " << attr << ", but should be "
                                                << space_dim << "!");
  if (read_attribute(domain, "hyperedge_dim", attr))
    hy_check(attr == hyEdge_dim, "hyperedge_dim in " << filename << " is " << attr
                                                     << ", but should be " << hyEdge_dim << "!");

  const std::array<hsize_t, 2> points_shape = shape(domain, "points");
  const std::array<hsize_t, 2> edges_shape = shape(domain, "edges");
  hy_check(points_shape[1] == space_dim, "points in " << filename << " have " << points_shape[1]
                                                      << " coordinates, but should have "
                                                      << space_dim << "!");
  hy_check(edges_shape[1] == n_verts, "edges in " << filename << " have " << edges_shape[1]
                                                  << " points, but should have " << n_verts << "!");
  const pt_index_t n_points = points_shape[0];
  const hyEdge_index_t n_hyEdges = edges_shape[0];
  const bool has_hypernodes = exists(domain, "hypernodes");
  hy_check(has_hypernodes || hyEdge_dim == 1,
           filename << " has no dataset hypernodes, which is only optional for hyperedge_dim 1.");

  hyNode_index_t n_hyNodes = n_points;
  if (has_hypernodes && !read_attribute(domain, "n_hypernodes", attr))
    attr = -1;
  DomainInfo<hyEdge_dim, space_dim, vectorT, pointT, hyEdge_index_t, hyNode_index_t, pt_index_t>
    domain_info(n_points, n_hyEdges, has_hypernodes && attr >= 0 ? attr : n_hyNodes, n_points);

  read(domain, "points", n_points, space_dim,
       reinterpret_cast<coord_t*>(domain_info.points.data()));
  read(domain, "edges", n_hyEdges, n_verts, domain_info.points_hyEdge.data()->data());
  read(domain, "types_faces", n_hyEdges, n_faces, domain_info.hyFaces_hyEdge.data()->data());
  if (has_hypernodes)
  {
    read(domain, "hypernodes", n_hyEdges, n_faces, domain_info.hyNodes_hyEdge.data()->data());
    if (attr < 0)  // no n_hypernodes attribute: count them
    {
      for (const auto& nodes : domain_info.hyNodes_hyEdge)
        for (const hyNode_index_t node : nodes)
          n_hyNodes = std::max<hyNode_index_t>(n_hyNodes, node + 1);
      domain_info.n_hyNodes = n_hyNodes;
    }
  }
  else
    for (hyEdge_index_t edge = 0; edge < n_hyEdges; ++edge)
      for (unsigned int i = 0; i < n_faces; ++i)
        domain_info.hyNodes_hyEdge[edge][i] = domain_info.points_hyEdge[edge][i];

  if (exists(domain, "properties"))
  {
    const std::array<hsize_t, 2> props_shape = shape(domain, "properties");
    hy_check(props_shape[0] == n_hyEdges,
             "properties in " << filename << " has " << props_shape[0] << " rows, but there "
                              << "are " << n_hyEdges << " hyperedges!");
    domain_info.n_properties = props_shape[1];
    std::vector<coord_t> flat(n_hyEdges * domain_info.n_properties);
    read(domain, "properties", n_hyEdges, domain_info.n_properties, flat.data());
    domain_info.hyEdge_properties.resize(n_hyEdges);
    for (hyEdge_index_t edge = 0; edge < n_hyEdges; ++edge)
      domain_info.hyEdge_properties[edge].assign(
        flat.begin() + edge * domain_info.n_properties,
        flat.begin() + (edge + 1) * domain_info.n_properties);
  }

  return domain_info;
}  // end of read_domain_hdf5

/*!*************************************************************************************************
 * \brief   Write a domain as a .geo.h5 file (cf. namespace HDF5Domain for the format).
 *
 * An existing file is overwritten. With \c vtkhdf (default), the file also contains the VTKHDF
 * view for ParaView, which exists for hyperedge dimensions 1 to 3 in up to 3 space dimensions
 * and is skipped otherwise.
 *
 * \param   filename        Name of the file to write.
 * \param   domain_info     Topological and geometrical information of the hypergraph.
 * \param   vtkhdf          Also write the VTKHDF view.
 **************************************************************************************************/
template <unsigned int hyEdge_dim,
          unsigned int space_dim,
          template <typename...>
          typename vectorT,
          typename pointT,
          typename hyEdge_index_t,
          typename hyNode_index_t,
          typename pt_index_t>
void write_domain_hdf5(const std::string& filename,
                       const DomainInfo<hyEdge_dim,
                                        space_dim,
                                        vectorT,
                                        pointT,
                                        hyEdge_index_t,
                                        hyNode_index_t,
                                        pt_index_t>& domain_info,
                       const bool vtkhdf = true)
{
  using namespace HDF5Domain;
  using coord_t = typename pointT::value_type;
  static_assert(sizeof(pointT) == space_dim * sizeof(coord_t),
                "Points must be contiguous arrays of their coordinates.");
  constexpr unsigned int n_verts = 1 << hyEdge_dim, n_faces = 2 * hyEdge_dim;
  const hsize_t n_points = domain_info.n_points, n_hyEdges = domain_info.n_hyEdges;

  Id file(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT), H5Fclose,
          "create file " + filename);
  Id domain(H5Gcreate2(file, "domain", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT), H5Gclose,
            "create group domain");
  write_attribute(domain, "space_dim", space_dim);
  write_attribute(domain, "hyperedge_dim", hyEdge_dim);
  write_attribute(domain, "n_hypernodes", (long long)domain_info.n_hyNodes);

  write(domain, "points", n_points, space_dim,
        reinterpret_cast<const coord_t*>(domain_info.points.data()));
  write(domain, "edges", n_hyEdges, n_verts, domain_info.points_hyEdge.data()->data());
  write(domain, "types_faces", n_hyEdges, n_faces, domain_info.hyFaces_hyEdge.data()->data());

  // hypernodes: only if they are not numbered like the points (never the case for a network read
  // from a Python generator, always for hyEdge_dim > 1, where the counts differ)
  bool nodes_like_points = hyEdge_dim == 1 && domain_info.n_hyNodes == domain_info.n_points;
  for (hyEdge_index_t edge = 0; nodes_like_points && edge < n_hyEdges; ++edge)
    for (unsigned int i = 0; i < n_faces; ++i)
      nodes_like_points &=
        domain_info.hyNodes_hyEdge[edge][i] == (hyNode_index_t)domain_info.points_hyEdge[edge][i];
  if (!nodes_like_points)
    write(domain, "hypernodes", n_hyEdges, n_faces, domain_info.hyNodes_hyEdge.data()->data());

  // types_points: the strongest (largest) type among the faces meeting at a point
  if (nodes_like_points)
  {
    std::vector<unsigned int> types_points(n_points, 0);
    for (hyEdge_index_t edge = 0; edge < n_hyEdges; ++edge)
      for (unsigned int i = 0; i < n_faces; ++i)
        types_points[domain_info.points_hyEdge[edge][i]] = std::max<unsigned int>(
          types_points[domain_info.points_hyEdge[edge][i]], domain_info.hyFaces_hyEdge[edge][i]);
    write(domain, "types_points", n_points, 1, types_points.data());
  }

  const bool has_properties = domain_info.n_properties > 0;
  if (has_properties)
  {
    std::vector<coord_t> flat;
    flat.reserve(n_hyEdges * domain_info.n_properties);
    for (const auto& props : domain_info.hyEdge_properties)
    {
      hy_check(
        props.size() == domain_info.n_properties,
        "Every hyperedge must have n_properties = " << domain_info.n_properties << " properties.");
      flat.insert(flat.end(), props.begin(), props.end());
    }
    write(domain, "properties", n_hyEdges, domain_info.n_properties, flat.data());
  }

  if constexpr (hyEdge_dim <= 3 && space_dim <= 3)
    if (vtkhdf)
      write_vtkhdf_view<hyEdge_dim, space_dim>(file, domain_info.points, domain_info.points_hyEdge,
                                               nodes_like_points, has_properties);
}  // end of write_domain_hdf5
