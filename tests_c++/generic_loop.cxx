// Design-validation harness for GlobalLoop::Generic: the consolidated loop must reproduce the
// established loops elementwise on the same problems (no driver is migrated yet).
//
// Part A: diffusion on a cubic grid -- Generic vs Elliptic (jacobian/residual/jacobian_mat and
//         the configure-time zero_indices vs the read_dirichlet_indices pattern).
//
// hy_check (not hy_assert) so the release build fails loudly too.

#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/local_solver/diffusion_uniform_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>

#include <cmath>
#include <map>
#include <numeric>
#include <vector>

using namespace std;

// deterministic pseudo-random trace values in [-1, 1] (no <random> to keep runs reproducible
// across standard libraries)
static double lcg_value(unsigned int& state)
{
  state = 1664525u * state + 1013904223u;
  return static_cast<double>(state) / 2147483648. - 1.;
}

// accumulate a COO triplet set into (row, col) -> summed value (matches the duplicate-summing
// semantics of PETSc MatSetValuesCOO / scipy.sparse)
template <typename MatT>
static map<pair<unsigned int, unsigned int>, double> coo_to_map(const MatT& mat)
{
  map<pair<unsigned int, unsigned int>, double> result;
  for (size_t i = 0; i < mat.value_vec.size(); ++i)
    result[{mat.row_vec[i], mat.col_vec[i]}] += mat.value_vec[i];
  return result;
}

static void check_equal(const vector<double>& a, const vector<double>& b, const char* what)
{
  hy_check(a.size() == b.size(), what << ": size mismatch " << a.size() << " vs " << b.size());
  for (size_t i = 0; i < a.size(); ++i)
    hy_check(abs(a[i] - b[i]) < 1e-12,
             what << ": entry " << i << " differs, " << a[i] << " vs " << b[i]);
}

static int test_part_a()
{
  const vector<unsigned int> num_elements = {4, 2, 2};

  using Elliptic = GlobalLoop::Elliptic<Topology::Cubic<1, 3>, Geometry::UnitCube<1, 3>,
                                        NodeDescriptor::Cubic<1, 3>,
                                        LocalSolver::DiffusionUniform<1, 1, 2, double>,
                                        vector<double> >;
  using Generic = GlobalLoop::Generic<Topology::Cubic<1, 3>, Geometry::UnitCube<1, 3>,
                                      NodeDescriptor::Cubic<1, 3>,
                                      LocalSolver::DiffusionUniform<1, 1, 2, double>,
                                      vector<double> >;

  Elliptic elliptic(num_elements, num_elements, 1.);
  Generic generic(num_elements, num_elements, 1.);

  const unsigned int n = elliptic.size_of_system();
  hy_check(generic.size_of_system() == n, "system sizes must agree");

  vector<double> x(n), out(n, 0.);
  unsigned int seed = 42u;
  for (auto& v : x)
    v = lcg_value(seed);

  // --- 1. empty zero set: operators must agree entrywise --------------------------------------
  generic.jacobian(x, out);
  check_equal(out, elliptic.trace_to_flux(x), "A1 jacobian vs trace_to_flux");

  generic.residual(x, out);
  check_equal(out, elliptic.residual_flux(x), "A1 residual vs residual_flux");

  sparse_mat<vector<double> > mat;
  generic.jacobian_mat(mat);
  {
    const auto map_g = coo_to_map(mat);
    const auto map_e = coo_to_map(elliptic.trace_to_flux_mat());
    for (const auto& [key, val] : map_e)
    {
      const auto it = map_g.find(key);
      const double val_g = (it == map_g.end()) ? 0. : it->second;
      hy_check(abs(val - val_g) < 1e-12, "A1 jacobian_mat: entry (" << key.first << ", "
               << key.second << ") differs, " << val_g << " vs " << val);
    }
    for (const auto& [key, val] : map_g)
      hy_check(map_e.count(key) || abs(val) < 1e-12,
               "A1 jacobian_mat: spurious entry (" << key.first << ", " << key.second << ")");
  }

  // --- 2. zero set = all cubic boundary types; cross-validate against the established
  //        read_dirichlet_indices pattern ------------------------------------------------------
  vector<unsigned int> boundary_types(26);
  iota(boundary_types.begin(), boundary_types.end(), 1u);
  generic.configure(boundary_types);

  const auto& zero_indices = generic.zero_indices_global();
  hy_check(!zero_indices.empty(), "A2: the cubic grid has boundary nodes");
  hy_check(zero_indices.size() < n, "A2: not every node is a boundary node");
  elliptic.read_dirichlet_indices(
    vector<unsigned int>(zero_indices.begin(), zero_indices.end()));

  generic.jacobian(x, out);
  check_equal(out, elliptic.trace_to_flux(x), "A2 jacobian vs trace_to_flux (zeroed)");

  generic.residual(x, out);
  check_equal(out, elliptic.residual_flux(x), "A2 residual vs residual_flux (zeroed)");

  // --- 3. jacobian_mat with zero set: rows/cols blanked, single diagonal 1 ---------------------
  generic.jacobian_mat(mat);  // buffer reuse across assemblies is part of the contract
  {
    const auto map_g = coo_to_map(mat);
    const auto map_e = coo_to_map(elliptic.trace_to_flux_mat());
    const auto is_zero = [&](unsigned int idx)
    { return binary_search(zero_indices.begin(), zero_indices.end(), idx); };

    for (const auto& [key, val] : map_g)
      if (is_zero(key.first) || is_zero(key.second))
        hy_check(abs(val - (key.first == key.second ? 1. : 0.)) < 1e-12,
                 "A3: zeroed row/col entry (" << key.first << ", " << key.second
                 << ") must be " << (key.first == key.second ? 1. : 0.) << ", got " << val);
    for (const auto idx : zero_indices)
    {
      const auto it = map_g.find({idx, idx});
      hy_check(it != map_g.end() && abs(it->second - 1.) < 1e-12,
               "A3: zero index " << idx << " needs a diagonal 1");
    }
    for (const auto& [key, val] : map_e)
      if (!is_zero(key.first) && !is_zero(key.second))
      {
        const auto it = map_g.find(key);
        const double val_g = (it == map_g.end()) ? 0. : it->second;
        hy_check(abs(val - val_g) < 1e-12, "A3: interior entry (" << key.first << ", "
                 << key.second << ") differs, " << val_g << " vs " << val);
      }
  }

  return 0;
}

int main()
{
  return test_part_a();
}
