// SparseLA::cg on small dense operators: every stop reason once, the initial guess, the monitor,
// and the legacy SparseLA::conjugate_gradient on top of it.
#include <HyperHDG/sparse_la.hxx>

#include <cmath>
#include <cstdio>
#include <limits>
#include <vector>

using Vec = std::vector<double>;

int failures = 0;
void check(const bool ok, const char* what)
{
  std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what);
  failures += !ok;
}

// The 1D Laplacian tridiag(-1, 2, -1) times scale: symmetric, definite of the sign of scale.
auto laplacian(const double scale)
{
  return [scale](const Vec& in, Vec& out)
  {
    const std::size_t n = in.size();
    for (std::size_t i = 0; i < n; ++i)
      out[i] = scale * (2. * in[i] - (i > 0 ? in[i - 1] : 0.) - (i + 1 < n ? in[i + 1] : 0.));
  };
}
double residual_norm(const Vec& b, const Vec& x, const double scale)
{
  Vec ax(x.size());
  laplacian(scale)(x, ax);
  double sum = 0.;
  for (std::size_t i = 0; i < x.size(); ++i)
    sum += (b[i] - ax[i]) * (b[i] - ax[i]);
  return std::sqrt(sum);
}

int main()
{
  const unsigned int n = 50;
  // A generic right-hand side: for b = (1, ..., 1), CG on this operator reaches ||r|| = 0 exactly.
  Vec b(n);
  double b_norm = 0.;
  for (unsigned int i = 0; i < n; ++i)
  {
    b[i] = 1. + std::sin(1. + i);
    b_norm += b[i] * b[i];
  }
  b_norm = std::sqrt(b_norm);
  SparseLA::SolverOptions options;
  options.rtol = 1e-10;

  {
    Vec x(n, 0.);
    unsigned int calls = 0;
    SparseLA::SolverOptions monitored = options;
    monitored.monitor = [&](unsigned int, double) { ++calls; };
    const auto result = SparseLA::cg(laplacian(1.), b, x, monitored);
    check(result.reason == SparseLA::StopReason::converged_rtol, "SPD: converged_rtol");
    check(result.iterations <= n, "SPD: at most n iterations");
    check(residual_norm(b, x, 1.) <= 1e-9 * b_norm, "SPD: the true residual is small");
    check(calls == result.iterations + 1, "the monitor sees every iteration and the end");
    check(result.rhs_norm == b_norm && result.initial_residual_norm == b_norm,
          "zero guess: r_0 = b");
  }
  {
    Vec x(n, 0.);
    const auto result = SparseLA::cg(laplacian(-1.), b, x, options);
    check(result.converged() && residual_norm(b, x, -1.) <= 1e-9 * b_norm,
          "negative definite: converged");
  }
  {
    Vec x(n, 0.);
    const auto first = SparseLA::cg(laplacian(1.), b, x, options);
    const auto again = SparseLA::cg(laplacian(1.), b, x, options);
    check(again.converged() && again.iterations <= 1 && again.initial_residual_norm < 1e-9,
          "the solution as initial guess: (almost) no iterations");
    check(first.iterations > again.iterations, "the initial guess is used");
  }
  {
    Vec x(n, 0.);
    SparseLA::SolverOptions few = options;
    few.max_iterations = 3;
    const auto result = SparseLA::cg(laplacian(1.), b, x, few);
    check(result.reason == SparseLA::StopReason::diverged_max_iterations && result.iterations == 3,
          "diverged_max_iterations after max_iterations");
  }
  {
    Vec x(n, 0.);
    const Vec zero(n, 0.);
    const auto result = SparseLA::cg(laplacian(1.), zero, x, options);
    check(result.reason == SparseLA::StopReason::converged_atol && result.iterations == 0,
          "zero right-hand side: converged_atol without iterations");
  }
  {
    Vec x(n, 0.);
    SparseLA::SolverOptions absolute = options;
    absolute.rtol = 0.;
    absolute.atol = 1e-6;
    const auto result = SparseLA::cg(laplacian(1.), b, x, absolute);
    check(result.reason == SparseLA::StopReason::converged_atol && result.residual_norm <= 1e-6,
          "converged_atol with rtol = 0");
  }
  {
    Vec x(n, 0.);
    SparseLA::SolverOptions tight = options;
    tight.dtol = 0.5;  // diverged as soon as ||r|| > ||r_0|| / 2 without having converged
    const auto result = SparseLA::cg(laplacian(1.), b, x, tight);
    check(result.reason == SparseLA::StopReason::diverged_dtol, "diverged_dtol");
  }
  {
    Vec x(n, 0.);
    const auto result =
      SparseLA::cg([](const Vec&, Vec& out) { out.assign(out.size(), 0.); }, b, x, options);
    check(result.reason == SparseLA::StopReason::diverged_breakdown, "zero operator: breakdown");
  }
  {
    Vec x(n, 0.);
    const auto result = SparseLA::cg(
      [](const Vec&, Vec& out)
      { out.assign(out.size(), std::numeric_limits<double>::quiet_NaN()); }, b, x, options);
    check(result.reason == SparseLA::StopReason::diverged_nan, "NaN operator: diverged_nan");
  }
  {
    struct Problem
    {
      Vec trace_to_flux(const Vec& in)
      {
        Vec out(in.size());
        laplacian(1.)(in, out);
        return out;
      }
    } problem;
    const Vec x = SparseLA::conjugate_gradient(b, problem, 0, 1e-10);
    check(residual_norm(b, x, 1.) <= 1e-9, "legacy conjugate_gradient: solves");
    bool thrown = false;
    try
    {
      SparseLA::conjugate_gradient(b, problem, 3, 1e-10);
    }
    catch (const SparseLA::SolveException&)
    {
      thrown = true;
    }
    check(thrown, "legacy conjugate_gradient: throws SolveException without convergence");
  }
  return failures;
}
