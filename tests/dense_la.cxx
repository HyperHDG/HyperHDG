// Point / SmallVec / SmallMat of dense_la.hxx in dimensions 1 to 10, double and float:
// construction, comparison, arithmetic, norms, the QR decomposition, and the handedness of the QR
// frame of an edge direction in 3D. The checks do not rely on hy_assert, so they also run with
// NDEBUG.
#include <HyperHDG/dense_la.hxx>

#include <array>
#include <cmath>
#include <cstdio>
#include <random>
#include <utility>

int failures = 0;
void check(const bool ok, const char* what, const unsigned int dim, const char* type)
{
  if (!ok)
    std::printf("FAIL  %s (dimension %u, %s)\n", what, dim, type);
  failures += !ok;
}
template <typename T>
bool almost_equal(const T a, const T b)
{
  return std::abs(a - b) < (T)1e-5;
}
// Every component of a equals f(component of b, component of c).
template <unsigned int dim, typename T, class F>
bool componentwise(const Point<dim, T>& a, const Point<dim, T>& b, const Point<dim, T>& c, F f)
{
  bool ok = true;
  for (unsigned int i = 0; i < dim; ++i)
    ok = ok && a[i] == f(b[i], c[i]);
  return ok;
}

template <unsigned int dim, typename T>
void test_point(const char* type)
{
  auto expect = [type](const bool ok, const char* what) { check(ok, what, dim, type); };
  std::mt19937 generator(dim);
  std::uniform_real_distribution<T> uniform(-10., 10.);
  auto nonzero = [&]()
  {
    T value;
    do
      value = uniform(generator);
    while (value == (T)0.);
    return value;
  };

  Point<dim, T> zero_point;
  SmallVec<dim, T> zero_vec;
  bool zeros = true;
  for (unsigned int i = 0; i < dim; ++i)
    zeros = zeros && zero_point[i] == (T)0. && zero_vec[i] == (T)0.;
  expect(zeros, "default construction gives zeros");

  Point<dim, T> a, b;
  for (unsigned int i = 0; i < dim; ++i)
    a[i] = nonzero(), b[i] = nonzero();
  const T s = nonzero();
  const Point<dim, T> ac(a), bc(b);

  // ---- construction and comparison
  expect(a == ac && !(a != ac) && !(a < ac) && !(ac < a), "a copy is equal, not unequal, not less");
  Point<dim, T> moved(std::move(Point<dim, T>(ac)));
  Point<dim, T> from_array = ac.data(), from_moved_array(std::move(Point<dim, T>(ac).data()));
  expect(moved == ac && from_array == ac && from_moved_array == ac,
         "move and array constructors give the same point");
  Point<dim, T> changed = ac;
  changed[dim - 1] -= (T)1.;
  expect(changed != ac && changed < ac && !(ac < changed), "a smaller last component: less");
  expect(moved == ac, "copies are independent");

  // ---- arithmetic
  auto plus = [](T x, T y) { return x + y; };
  auto minus = [](T x, T y) { return x - y; };
  auto times = [](T x, T y) { return x * y; };
  auto divided = [](T x, T y) { return x / y; };
  Point<dim, T> scalar_point;
  for (unsigned int i = 0; i < dim; ++i)
    scalar_point[i] = s;
  Point<dim, T> c = ac;
  expect(componentwise(c += s, ac, scalar_point, plus), "+= scalar");
  c = ac;
  expect(componentwise(c -= s, ac, scalar_point, minus), "-= scalar");
  c = ac;
  expect(componentwise(c *= s, ac, scalar_point, times), "*= scalar");
  c = ac;
  expect(componentwise(c /= s, ac, scalar_point, divided), "/= scalar");
  c = ac;
  expect(componentwise(c += bc, ac, bc, plus), "+= point");
  c = ac;
  expect(componentwise(c -= bc, ac, bc, minus), "-= point");
  c = ac;
  expect(componentwise(c *= bc, ac, bc, times), "*= point");
  c = ac;
  expect(componentwise(c /= bc, ac, bc, divided), "/= point");
  expect(componentwise(Point<dim, T>(ac + bc), ac, bc, plus), "point + point");
  expect(componentwise(Point<dim, T>(ac - bc), ac, bc, minus), "point - point");
  expect(componentwise(Point<dim, T>(hada_prod(ac, bc)), ac, bc, times), "hada_prod");
  expect(componentwise(Point<dim, T>(hada_divi(ac, bc)), ac, bc, divided), "hada_divi");
  expect(componentwise(Point<dim, T>(s + ac), scalar_point, ac, plus) &&
           componentwise(Point<dim, T>(ac + s), scalar_point, ac, plus),
         "scalar + point, point + scalar");
  expect(componentwise(Point<dim, T>(s - ac), scalar_point, ac, minus) &&
           componentwise(Point<dim, T>(ac - s), ac, scalar_point, minus),
         "scalar - point, point - scalar");
  expect(componentwise(Point<dim, T>(s * ac), scalar_point, ac, times) &&
           componentwise(Point<dim, T>(ac * s), scalar_point, ac, times),
         "scalar * point, point * scalar");
  expect(componentwise(Point<dim, T>(s / ac), scalar_point, ac, divided) &&
           componentwise(Point<dim, T>(ac / s), ac, scalar_point, divided),
         "scalar / point, point / scalar");
  T product = 0.;
  for (unsigned int i = 0; i < dim; ++i)
    product += ac[i] * bc[i];
  expect(scalar_product(ac, bc) == product, "scalar_product");

  // ---- norms
  expect(almost_equal(norm_1(ac), norm_p(ac, (T)1.)) && almost_equal(norm_2(ac), norm_p(ac, (T)2.)),
         "norm_1 and norm_2 equal norm_p with p = 1, 2");
  bool bounded = true;
  for (unsigned int i = 0; i < dim; ++i)
    bounded = bounded && std::abs(ac[i]) <= norm_infty(ac);
  expect(bounded, "components bounded by norm_infty");

  // ---- QR: the first column of Q is the normalized first column, the columns are orthonormal
  SmallMat<dim, 2, T> matrix;
  matrix.set_column(0, ac);
  matrix.set_column(1, bc);
  const SmallMat<dim, dim, T> q = qr_decomp_q(matrix);
  expect(norm_2(ac / norm_2(ac) - q.get_column(0)) < (T)1e-5 ||
           norm_2(ac / norm_2(ac) + q.get_column(0)) < (T)1e-5,
         "QR: first column of Q is the normalized first column");
  bool orthonormal = true;
  for (unsigned int i = 0; i < dim; ++i)
  {
    orthonormal = orthonormal && almost_equal(norm_2(q.get_column(i)), (T)1.);
    for (unsigned int j = 0; j < i; ++j)
      orthonormal =
        orthonormal && almost_equal(scalar_product(q.get_column(i), q.get_column(j)), (T)0.);
  }
  expect(orthonormal, "QR: the columns of Q are orthonormal");
}

// The edge-local frame (inner_normal(0), outer_normal(0), outer_normal(1)) is the Q of the QR
// decomposition of the edge direction. It must be right-handed, det(Q) = +1, for every direction,
// in particular for axis-aligned ones, where LAPACK emits a trivial Householder reflector.
void test_qr_handedness()
{
  const std::array<std::array<double, 3>, 8> directions = {{{{1., 0., 0.}},
                                                            {{-1., 0., 0.}},
                                                            {{0., 1., 0.}},
                                                            {{0., -1., 0.}},
                                                            {{0., 0., 1.}},
                                                            {{0., 0., -1.}},
                                                            {{1., 1., 1.}},
                                                            {{2., -3., 5.}}}};
  for (const auto& dir : directions)
  {
    SmallMat<3, 1, double> matrix;
    for (unsigned int i = 0; i < 3; ++i)
      matrix(i, 0) = dir[i];
    SmallSquareMat<3, double> q;
    SmallSquareMat<1, double> r;
    qr_decomp(matrix, q, r);
    const double det = q(0, 0) * (q(1, 1) * q(2, 2) - q(1, 2) * q(2, 1)) -
                       q(0, 1) * (q(1, 0) * q(2, 2) - q(1, 2) * q(2, 0)) +
                       q(0, 2) * (q(1, 0) * q(2, 1) - q(1, 1) * q(2, 0));
    if (!(det > 0.5))
      std::printf("FAIL  QR frame not right-handed: det = %g for (%g, %g, %g)\n", det, dir[0],
                  dir[1], dir[2]);
    failures += !(det > 0.5);
  }
}

int main()
{
  [&]<unsigned int... d>(std::integer_sequence<unsigned int, d...>)
  {
    (test_point<d + 1, double>("double"), ...);
    (test_point<d + 1, float>("float"), ...);
  }(std::make_integer_sequence<unsigned int, 10>());
  test_qr_handedness();
  std::printf("%s: %d failed checks\n", failures ? "FAILED" : "ok", failures);
  return failures;
}
