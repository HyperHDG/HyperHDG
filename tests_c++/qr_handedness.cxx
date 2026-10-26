#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/hy_assert.hxx>

#include <array>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>

// The edge-local frame (inner_normal(0), outer_normal(0), outer_normal(1)) is the QR matrix
// Q of the edge direction.  It must be right-handed, i.e. det(Q) = +1, for EVERY direction --
// in particular for axis-aligned ones, where LAPACK emits a trivial Householder reflector.
int main()
{
  const std::array<std::array<double, 3>, 8> directions = {{
    {{1., 0., 0.}}, {{-1., 0., 0.}}, {{0., 1., 0.}}, {{0., -1., 0.}},
    {{0., 0., 1.}}, {{0., 0., -1.}}, {{1., 1., 1.}}, {{2., -3., 5.}}}};

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

    hy_check(det > 0.5, "QR frame must be right-handed (det(Q) = +1), but det = "
                          << det << " for direction (" << dir[0] << ", " << dir[1] << ", "
                          << dir[2] << ").");
  }

  return 0;
}
