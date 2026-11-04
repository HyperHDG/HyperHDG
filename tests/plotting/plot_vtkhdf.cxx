#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/read_domain_hdf5.hxx>

#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// End-to-end check of the VTKHDF plot (include/HyperHDG/plot.hxx, PlotOptions::vtkhdf): every
// hypernode is Dirichlet with the value time * x, so the exact LDG-H solution of degree one is
// u = time * x with flux q = -time * e_x, and the plotted values must be these, step by step.
// Three steps are plotted into one file, which is then read back with HDF5 and checked
// against the coordinates it contains: steps and times, the values block of every step, the
// mesh arrays (offsets, types, connectivity) and the Steps bookkeeping ParaView reads.
// tests_c++/plot_vtkhdf_paraview.py opens the same files through ParaView's reader.

template <unsigned int space_dimT, typename param_float_t = double>
struct TestParameters
{
  static constexpr std::array<unsigned int, 27U> dirichlet_nodes{
    0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
    14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 1.;
  }
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>& point,
                                       const param_float_t time = 0.)
  {
    return time * point[0];
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

template <typename T>
std::vector<T> read_1d(const hid_t loc, const char* name)
{
  const auto dims = HDF5Domain::shape(loc, name);
  std::vector<T> data(dims[0] * dims[1]);
  HDF5Domain::read(loc, name, dims[0], dims[1], data.data());
  return data;
}

template <unsigned int hyEdge_dim, unsigned int space_dim>
std::string do_test()
{
  const std::vector<double> times = {0., 0.5, 1.};
  const std::vector<unsigned int> num_elements(space_dim, 2U);
  GlobalLoop::Elliptic<
    Topology::Cubic<hyEdge_dim, space_dim>, Geometry::UnitCube<hyEdge_dim, space_dim, double>,
    NodeDescriptor::Cubic<hyEdge_dim, space_dim>,
    LocalSolver::Diffusion<hyEdge_dim, 1, 2, TestParameters, double>, std::vector<double> >
    problem(num_elements, 1.);
  const std::string name =
    "plot_vtkhdf-" + std::to_string(hyEdge_dim) + "-" + std::to_string(space_dim);
  problem.plot_option("fileName", name);
  problem.plot_option("fileEnding", "vtkhdf");
  const std::vector<double> lambda = problem.zero_vector();
  for (const double time : times)
    problem.plot_solution(lambda, time);  // fileNumber increments: mesh, then steps 0, 1, 2

  constexpr unsigned int n_verts = 1 << hyEdge_dim, n_comp = hyEdge_dim + 1;
  constexpr unsigned char vtk_type = hyEdge_dim == 1 ? 3 : hyEdge_dim == 2 ? 8 : 11;
  const std::string file = "output/" + name + ".vtkhdf";
  using HDF5Domain::Id;
  Id h5(H5Fopen(file.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT), H5Fclose, "open " + file);
  Id root(H5Gopen2(h5, "VTKHDF", H5P_DEFAULT), H5Gclose, "open VTKHDF");
  Id steps(H5Gopen2(root, "Steps", H5P_DEFAULT), H5Gclose, "open Steps");
  Id pdata(H5Gopen2(root, "PointData", H5P_DEFAULT), H5Gclose, "open PointData");
  Id pdo(H5Gopen2(steps, "PointDataOffsets", H5P_DEFAULT), H5Gclose, "open PointDataOffsets");

  long long n_steps = -1;
  hy_check(
    HDF5Domain::read_attribute(steps, "NSteps", n_steps) && n_steps == (long long)times.size(),
    file << ": NSteps is " << n_steps << ", expected " << times.size() << ".");
  const auto step_values = read_1d<double>(steps, "Values");
  for (unsigned int s = 0; s < times.size(); ++s)
    hy_check(step_values[s] == times[s], file << ": step " << s << " has time " << step_values[s]);

  const auto n_points = read_1d<long long>(root, "NumberOfPoints")[0];
  const auto n_cells = read_1d<long long>(root, "NumberOfCells")[0];
  const auto n_conn = read_1d<long long>(root, "NumberOfConnectivityIds")[0];
  hy_check(n_points == n_verts * n_cells && n_conn == n_verts * n_cells && n_cells > 0,
           file << ": " << n_points << " points, " << n_cells << " cells, " << n_conn << " ids.");

  const auto points = read_1d<float>(root, "Points");
  const auto connectivity = read_1d<long long>(root, "Connectivity");
  const auto offsets = read_1d<long long>(root, "Offsets");
  const auto types = read_1d<unsigned char>(root, "Types");
  hy_check((long long)offsets.size() == n_cells + 1 && (long long)types.size() == n_cells,
           file << ": sizes of Offsets / Types.");
  for (long long cell = 0; cell < n_cells; ++cell)
  {
    hy_check(offsets[cell] == cell * n_verts && offsets[cell + 1] == (cell + 1) * n_verts &&
               types[cell] == vtk_type,
             file << ": cell " << cell << " has wrong offset or type.");
    for (unsigned int v = 0; v < n_verts; ++v)
    {
      const long long id = connectivity[cell * n_verts + v];
      hy_check(id >= 0 && id < n_points, file << ": connectivity id " << id << " out of range.");
      for (unsigned int w = 0; w < v; ++w)
        hy_check(id != connectivity[cell * n_verts + w],
                 file << ": cell " << cell << " repeats a vertex.");
    }
  }

  const auto values_shape = HDF5Domain::shape(pdata, "values");
  hy_check(values_shape[0] == (hsize_t)n_points * times.size() && values_shape[1] == n_comp,
           file << ": values has shape (" << values_shape[0] << ", " << values_shape[1] << ").");
  const auto values = read_1d<float>(pdata, "values");
  const auto values_offsets = read_1d<long long>(pdo, "values");
  for (unsigned int s = 0; s < times.size(); ++s)
  {
    hy_check(values_offsets[s] == (long long)s * n_points, file << ": PointDataOffsets/values.");
    for (long long p = 0; p < n_points; ++p)
    {
      const float* row = values.data() + ((hsize_t)s * n_points + p) * n_comp;
      // u is the last component; the flux is edge-local, so it equals -t e_x only if the
      // hyperedges span the space (otherwise it vanishes on edges orthogonal to x)
      const float x = points[3 * p], u = row[n_comp - 1], q0 = row[0];
      hy_check(std::abs(u - times[s] * x) < 1e-5, file << ": step " << s << ", point " << p
                                                       << " (x = " << x << "): u = " << u
                                                       << ", expected " << times[s] * x);
      if constexpr (hyEdge_dim == space_dim)
        hy_check(std::abs(q0 + times[s]) < 1e-5, file << ": step " << s << ", point " << p
                                                      << ": q_x = " << q0 << ", expected "
                                                      << -times[s]);
    }
  }
  std::cout << "  " << file << ": " << n_cells << " cells, " << times.size() << " steps ok"
            << std::endl;
  return file;
}

int main()
{
  do_test<1, 1>();
  do_test<1, 2>();
  do_test<1, 3>();
  do_test<2, 2>();
  do_test<2, 3>();
  do_test<3, 3>();
  std::cout << "plot_vtkhdf: all checks passed" << std::endl;
  return 0;
}
