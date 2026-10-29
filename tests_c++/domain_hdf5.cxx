#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/read_domain.hxx>

#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

// .geo.h5 domain files (include/HyperHDG/read_domain_hdf5.hxx): a domain read from a .geo file,
// written with write_domain_hdf5 and read back through read_domain is the same domain, for a
// network in 3d (with properties), a network in 2d (padded VTKHDF points), hyperedges of
// dimension 2 (a hypernodes dataset) and a network in 5d (no VTKHDF view). The file
// domains/cross.geo.h5, written with h5py in the layout of the Python network generators,
// equals domains/cross.geo.

template <typename InfoT>
void expect_equal(const InfoT& a, const InfoT& b, const std::string& what)
{
  hy_check(a.n_hyEdges == b.n_hyEdges && a.n_hyNodes == b.n_hyNodes && a.n_points == b.n_points,
           what << ": sizes differ (" << a.n_hyEdges << "/" << b.n_hyEdges << " hyperedges, "
                << a.n_hyNodes << "/" << b.n_hyNodes << " hypernodes, " << a.n_points << "/"
                << b.n_points << " points).");
  for (unsigned int pt = 0; pt < a.n_points; ++pt)
    for (unsigned int dim = 0; dim < a.points[pt].size(); ++dim)
      hy_check(a.points[pt][dim] == b.points[pt][dim], what << ": point " << pt << " differs.");
  for (unsigned int edge = 0; edge < a.n_hyEdges; ++edge)
  {
    hy_check(a.hyNodes_hyEdge[edge] == b.hyNodes_hyEdge[edge],
             what << ": hypernodes of hyperedge " << edge << " differ.");
    hy_check(a.hyFaces_hyEdge[edge] == b.hyFaces_hyEdge[edge],
             what << ": face types of hyperedge " << edge << " differ.");
    hy_check(a.points_hyEdge[edge] == b.points_hyEdge[edge],
             what << ": points of hyperedge " << edge << " differ.");
  }
  hy_check(a.n_properties == b.n_properties && a.hyEdge_properties == b.hyEdge_properties,
           what << ": properties differ.");
}

bool has_group(const std::string& filename, const char* group)
{
  HDF5Domain::Id file(H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose,
                      "open " + filename);
  return HDF5Domain::exists(file, group);
}

template <unsigned int hyEdge_dim, unsigned int space_dim>
void roundtrip(const std::string& geo, const bool with_properties, const bool expect_vtkhdf)
{
  auto info = read_domain<hyEdge_dim, space_dim>(geo);
  if (with_properties)
  {
    info.n_properties = 2;
    info.hyEdge_properties.resize(info.n_hyEdges);
    for (unsigned int edge = 0; edge < info.n_hyEdges; ++edge)
      info.hyEdge_properties[edge] = {double(edge), 10. * edge + 0.5};
  }

  const std::string h5 = (std::filesystem::temp_directory_path() /
                          (std::filesystem::path(geo).stem().string() + "_domain_hdf5_test.geo.h5"))
                           .string();
  write_domain_hdf5(h5, info);
  hy_check(is_domain_hdf5(h5), h5 << " is not recognized as a .geo.h5 file.");
  hy_check(has_group(h5, "VTKHDF") == expect_vtkhdf,
           h5 << (expect_vtkhdf ? " lacks" : " unexpectedly has") << " the VTKHDF view.");
  hy_check(has_group(h5, "/domain/hypernodes") == (hyEdge_dim > 1),
           h5 << ": hypernodes dataset expected exactly for hyEdge_dim > 1.");

  const auto back = read_domain<hyEdge_dim, space_dim>(h5);
  expect_equal(info, back, geo + " -> " + h5);
  std::filesystem::remove(h5);
  std::cout << "  " << geo << ": roundtrip ok" << std::endl;
}

int main()
{
  roundtrip<1, 3>("domains/cross.geo", true, true);
  roundtrip<1, 2>("domains/injection_test.geo", false, true);
  roundtrip<2, 3>("domains/leVeque_hg.geo", false, true);
  roundtrip<1, 5>("domains/simplex_1_5.geo", false, false);

  hy_check(!is_domain_hdf5("domains/cross.geo"), "A .geo file must not pass as HDF5.");
  expect_equal(read_domain<1, 3>("domains/cross.geo"), read_domain<1, 3>("domains/cross.geo.h5"),
               "domains/cross.geo vs domains/cross.geo.h5");
  std::cout << "  domains/cross.geo.h5 (h5py) equals domains/cross.geo" << std::endl;

  std::cout << "domain_hdf5: all checks passed" << std::endl;
  return 0;
}
