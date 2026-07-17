// Design-validation harness for GlobalLoop::Generic: the consolidated loop must reproduce the
// established loops elementwise on the same problems (no driver is migrated yet).
//
// Part A: diffusion on a cubic grid -- Generic vs Elliptic (jacobian/residual/jacobian_mat and
//         the configure-time zero_indices vs the read_dirichlet_indices pattern).
// Part B: TimoshenkoWave (wave4) on a file domain -- Generic vs Hyperbolic, including the Gauss
//         step protocol: postprocess with stage >= 0 vs set_data, stage == -1 vs finalize_step;
//         states compared via residuals and errors (no linear solve needed).
//
// hy_check (not hy_assert) so the release build fails loudly too.

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/global_loop/generic.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>
#include <HyperHDG/local_solver/diffusion_uniform_ldgh.hxx>
#include <HyperHDG/local_solver/timowave.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>

#include "../experiments/timowave4.hxx"

// the file-domain dispatcher (read_domain.hxx) checks PETSC_COMM_WORLD under the PETSc build,
// which is valid only after PetscInitialize
#if defined(HYPERHDG_PETSC)
#include <petscsys.h>
#elif defined(HYPERHDG_MPI)
#include <mpi.h>
#endif

#include <cmath>
#include <complex>
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
static auto coo_to_map(const MatT& mat)
{
  using value_t = typename std::decay_t<decltype(mat.value_vec)>::value_type;
  map<pair<unsigned int, unsigned int>, value_t> result;
  for (size_t i = 0; i < mat.value_vec.size(); ++i)
    result[{mat.row_vec[i], mat.col_vec[i]}] += mat.value_vec[i];
  return result;
}

// compare two COO maps entrywise (missing entries count as zero)
template <typename MapT>
static void check_mat_equal(const MapT& map_a, const MapT& map_b, const char* what,
                            const double tol = 1e-12)
{
  const auto value_or_zero = [](const MapT& m, const auto& key)
  {
    const auto it = m.find(key);
    return it == m.end() ? typename MapT::mapped_type(0.) : it->second;
  };
  for (const auto& [key, val] : map_a)
    hy_check(abs(val - value_or_zero(map_b, key)) < tol,
             what << ": entry (" << key.first << ", " << key.second << ") differs");
  for (const auto& [key, val] : map_b)
    hy_check(abs(val - value_or_zero(map_a, key)) < tol,
             what << ": entry (" << key.first << ", " << key.second << ") differs");
}

template <typename Scalar>
static void check_equal(const vector<Scalar>& a, const vector<Scalar>& b, const char* what,
                        const double tol = 1e-12)
{
  hy_check(a.size() == b.size(), what << ": size mismatch " << a.size() << " vs " << b.size());
  for (size_t i = 0; i < a.size(); ++i)
    hy_check(abs(a[i] - b[i]) < tol,
             what << ": entry " << i << " differs, " << a[i] << " vs " << b[i]);
}

template <typename Scalar>
static Scalar random_scalar(unsigned int& seed)
{
  if constexpr (requires { Scalar{0., 0.}.imag(); })
  {
    const double re = lcg_value(seed);
    return Scalar{re, lcg_value(seed)};
  }
  else
    return lcg_value(seed);
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
  check_mat_equal(coo_to_map(mat), coo_to_map(elliptic.trace_to_flux_mat()), "A1 jacobian_mat");

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

// Part B: Generic vs Hyperbolic on TimoshenkoWave/wave4. The step protocol is solve-free: both
// loops receive the SAME synthetic stage trace, so equal advanced states <=> equal residuals
// (the residual reads the recombined per-edge state) and equal errors (which read it too).
template <unsigned int poly_deg, unsigned int stages, typename Scalar>
static int test_part_b(const string& domain)
{
  using LSol =
    LocalSolver::TimoshenkoWave<1, 3, poly_deg, 2 * poly_deg, TestTimoWave4, double, stages>;
  using Hyperbolic = GlobalLoop::Hyperbolic<Topology::File<1, 3>, Geometry::File<1, 3>,
                                            NodeDescriptor::File<1, 3>, LSol, vector<Scalar> >;
  using Generic = GlobalLoop::Generic<Topology::File<1, 3>, Geometry::File<1, 3>,
                                      NodeDescriptor::File<1, 3>, LSol, vector<Scalar> >;

  const vector<double> vals = {1., 0.03125};  // {tau, dt}
  const double dt = vals[1];
  Hyperbolic hyp(domain, vals);
  Generic gen(domain, vals);

  const unsigned int n = hyp.size_of_system();
  hy_check(gen.size_of_system() == n, "B: system sizes must agree");

  // initial state at t = 0 (also sets up the per-edge state on both loops)
  const vector<Scalar> init_h = hyp.make_initial(hyp.zero_vector());
  vector<Scalar> init_g = gen.zero_vector();
  gen.initialize(init_g, 0.);
  check_equal(init_g, init_h, "B initialize vs make_initial", 1e-10);

  unsigned int seed = 7u;
  vector<Scalar> x(n), out_g(n, 0.), out_h(n, 0.), out_j(n, 0.), zero(n, 0.);
  for (auto& v : x)
    v = random_scalar<Scalar>(seed);

  sparse_mat<vector<Scalar> > mat;
  for (int rep = 0; rep < static_cast<int>(Hyperbolic::n_gauss_reps()); ++rep)
  {
    const Gauss::StageTime st{dt, rep};

    gen.residual(x, out_g, st);
    hyp.residual_flux2(x, out_h, st);
    check_equal(out_g, out_h, "B residual vs residual_flux2", 1e-9);

    // jacobian: A x = residual(x) - residual(0) (the residual is affine in the trace and both
    // calls share the zeroing), plus the direct apply where Hyperbolic offers it (s = 1)
    gen.jacobian(x, out_j, st);
    gen.residual(zero, out_g, st);
    for (unsigned int i = 0; i < n; ++i)
      out_h[i] -= out_g[i];
    check_equal(out_j, out_h, "B jacobian vs residual difference", 1e-9);
    if constexpr (stages == 1)
      check_equal(out_j, hyp.trace_to_flux(x, dt), "B jacobian vs trace_to_flux", 1e-9);

    gen.jacobian_mat(mat, st);  // reused buffer across reps: part of the contract
    check_mat_equal(coo_to_map(mat),
                    coo_to_map(hyp.template trace_to_flux_mat<unsigned int, vector<Scalar> >(st)),
                    "B jacobian_mat vs trace_to_flux_mat", 1e-9);
  }

  // Gauss step protocol: identical synthetic stage traces into both loops, then the old
  // finalize_step path vs the new postprocess(stage = -1) sentinel.
  vector<Scalar> zeta(n);
  for (auto& v : zeta)
    v = random_scalar<Scalar>(seed);
  for (int rep = 0; rep < static_cast<int>(Hyperbolic::n_gauss_reps()); ++rep)
  {
    hyp.set_data(zeta, Gauss::StageTime{dt, rep});
    gen.postprocess(zeta, Gauss::StageTime{dt, rep});
  }
  hyp.finalize_step();
  gen.postprocess(zeta, Gauss::StageTime{dt, -1});

  gen.residual(zero, out_g, Gauss::StageTime{2 * dt, 0});
  hyp.residual_flux2(zero, out_h, Gauss::StageTime{2 * dt, 0});
  check_equal(out_g, out_h, "B advanced state via residual", 1e-9);
  check_equal(gen.errors(zero, 2 * dt), hyp.errors(zero, 2 * dt), "B advanced state via errors",
              1e-10);

  return 0;
}

int main()
{
#if defined(HYPERHDG_PETSC)
  PetscInitialize(nullptr, nullptr, nullptr, nullptr);
#elif defined(HYPERHDG_MPI)
  MPI_Init(nullptr, nullptr);
#endif
  int result = test_part_a();
  result += test_part_b<1, 1, double>("domains/single1.geo");
  result += test_part_b<1, 1, double>("domains/cross2.geo");
  result += test_part_b<2, 2, complex<double> >("domains/single1.geo");
#if defined(HYPERHDG_PETSC)
  PetscFinalize();
#elif defined(HYPERHDG_MPI)
  MPI_Finalize();
#endif
  return result;
}
