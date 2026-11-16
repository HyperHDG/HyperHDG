// Timoshenko beam network waves with GlobalLoop::Hyperbolic on the cross of domains/cross2.geo (run
// in the source tree) with the manufactured solution wave4 of python/parameters/timowave.hxx: HDG
// in space, s = (degree + 1) / 2 Gauss stages in time, solve up to T = 1 on 4, 8, 16 subintervals
// per beam and check that the relative L2 error of the state converges with order degree + 1 on the
// last refinement. The time steps are small enough (temporal order 2s = degree + 1) to keep the
// temporal error below the spatial one. The spatial phase px = -pi/2 makes the Dirichlet data
// vanish, so that Gauss collocation keeps its full order (no stiff order reduction). Every Gauss
// stage representative's condensed system (complex for s >= 2) is assembled and factorized once
// (dense LU) and solved every time step. One executable per DEGREE, set by CMake
// (tests/CMakeLists.txt). The records are checked against tests/expected/timowave_p<DEGREE>.yaml
// (tests/records.hxx).
#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/global_loop/hyperbolic.hxx>
#include <HyperHDG/local_solver/timowave.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/file.hxx>
#include <HyperHDG/wrapper/lapack.hxx>

#include <python/parameters/timowave.hxx>
#include <tests/records.hxx>

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <string>
#include <type_traits>
#include <vector>

// Set per executable by tests/CMakeLists.txt; the default only serves editors and linters that
// compile this file without the CMake flags.
#ifndef DEGREE
#define DEGREE 1
#endif

constexpr unsigned int degree = DEGREE;
constexpr unsigned int stages = (degree + 1) / 2;
constexpr unsigned int n_steps = 256;
constexpr double end_time = 1.;

// The stage systems of s >= 2 Gauss stages are complex.
using scalar = std::conditional_t<stages == 1, double, std::complex<double>>;
using Loop = GlobalLoop::Hyperbolic<
  Topology::File<1, 3>,
  Geometry::File<1, 3>,
  NodeDescriptor::File<1, 3>,
  LocalSolver::TimoshenkoWave<1, 3, degree, 2 * degree, TestTimoWave4, double, stages>,
  std::vector<scalar>>;

// Dense LU factors of the condensed system of one Gauss stage representative.
struct Factored
{
  std::vector<scalar> lu;
  std::vector<int> ipiv;
};

Factored factorize(Loop& loop, const int rep)
{
  const unsigned int size = loop.size_of_system();
  const auto mat = loop.trace_to_flux_mat(Gauss::StageTime{0., rep});
  Factored result{std::vector<scalar>((std::size_t)size * size, 0.), std::vector<int>(size)};
  for (std::size_t k = 0; k < mat.value_vec.size(); ++k)  // column-major, duplicates summed
    result.lu[(std::size_t)mat.col_vec[k] * size + mat.row_vec[k]] += mat.value_vec[k];
  Wrapper::lapack_factorize(size, result.lu.data(), result.ipiv.data());
  return result;
}

int main()
{
  TestTimoWave4<3, double>::px = -M_PI / 2;
  const double dt = end_time / n_steps;

  int failures = 0;
  double previous = 0.;
  for (unsigned int n = 4; n <= 16; n *= 2)
  {
    Loop loop(std::string("domains/cross2.geo"), std::vector<double>{1., dt});
    loop.set_refinement(n);

    const unsigned int size = loop.size_of_system();
    const std::vector<scalar> zero = loop.zero_vector();
    std::vector<Factored> factored;
    for (unsigned int rep = 0; rep < Loop::n_gauss_reps(); ++rep)
      factored.push_back(factorize(loop, rep));

    // Maxima over all time steps of the errors and of the norms of the analytic solution; the
    // components are the state (L2), the trace (length-weighted skeleton norm) and the dual pair
    // (n, m) (L2).
    std::vector<double> error(3, 0.), norm(3, 0.);
    std::vector<scalar> trace = zero;
    auto record_maxima = [&](const double time)
    {
      const auto e = loop.errors(trace, time), a = loop.norms(trace, time);
      for (unsigned int i = 0; i < 3; ++i)
      {
        error[i] = std::max(error[i], (double)e[i]);
        norm[i] = std::max(norm[i], (double)a[i]);
      }
    };

    loop.make_initial(zero, 0.);
    record_maxima(0.);
    for (unsigned int step = 1; step <= n_steps; ++step)
    {
      const double time = step * dt;
      for (unsigned int rep = 0; rep < Loop::n_gauss_reps(); ++rep)
      {
        const Gauss::StageTime stage_time{time, (int)rep};
        trace = loop.residual_flux(zero, stage_time);
        for (scalar& value : trace)
          value = -value;
        Wrapper::lapack_solve_factored(size, 1, factored[rep].lu.data(), factored[rep].ipiv.data(),
                                       trace.data());
        loop.set_data(trace, stage_time);
      }
      loop.finalize_step();
      record_maxima(time);
    }

    const double e_rel = error[0] / norm[0], e_trace = error[1] / norm[1],
                 e_dual = error[2] / norm[2];
    // At degree 5, e_trace (down to 3e-11) is at the round-off of the dense LU solves and depends
    // on the BLAS (MKL and OpenBLAS differ by 1e-2 relative at n = 16): printed, but not recorded.
    if constexpr (degree < 5)
      tests::record(
        "- {degree: %u, stages: %u, n: %u, nt: %u, dofs: %u, e_rel: %.6e, e_trace: %.6e, "
        "e_dual: %.6e}\n",
        degree, stages, n, n_steps, size, e_rel, e_trace, e_dual);
    else
    {
      tests::record(
        "- {degree: %u, stages: %u, n: %u, nt: %u, dofs: %u, e_rel: %.6e, e_dual: %.6e}\n", degree,
        stages, n, n_steps, size, e_rel, e_dual);
      std::printf("  e_trace: %.6e\n", e_trace);
    }
    if (!std::isfinite(e_rel) || !std::isfinite(e_trace) || !std::isfinite(e_dual))
    {
      std::printf("  FAIL: an error is not finite\n");
      ++failures;
    }

    // Only the last rate is checked: the coarse meshes are preasymptotic.
    if (n == 16)
    {
      const double rate = std::log2(previous / e_rel);
      std::printf("  rate: %.3f (expected at least %.1f)\n", rate, degree + 1 - 0.3);
      if (!(rate >= degree + 1 - 0.3))
        ++failures;
    }
    previous = e_rel;
  }
  // At degree 5, the round-off also reaches e_dual: Accelerate (macOS) differs from MKL by 1e-4
  // relative at n = 16.
  return tests::check(failures, degree < 5 ? 1e-5 : 1e-3);
}
