/*!*************************************************************************************************
 * \file    stiffness_matrix.cxx
 * \brief   Check the assembly of the condensed (stiffness) matrix.
 *
 * The matrix assembled by trace_to_flux_mat (with the local matrix being computed per hyperedge
 * and with the local matrix being computed only once) is compared to the matrix whose columns are
 * obtained by applying the matrix-free operator trace_to_flux to all unit vectors. This is done
 * - for the local solver Diffusion, which provides trace_to_flux_mat (one factorization per
 *   hyperedge), on a uniform mesh with mixed boundary conditions, on a uniform mesh with a space
 *   dependent diffusion coefficient, and on a non-uniform hypergraph, and
 * - for the local solver DiffusionUniform, which does not provide trace_to_flux_mat (fallback).
 *
 * If the executable is called with the argument "misuse", factorize_once is used for a space
 * dependent diffusion coefficient. In debug mode, this has to be detected (and the program aborts).
 **************************************************************************************************/

#include <HyperHDG/geometry/file.hxx>
#include <HyperHDG/geometry/unit_cube.hxx>
#include <HyperHDG/global_loop/elliptic.hxx>
#include <HyperHDG/local_solver/diffusion_ldgh.hxx>
#include <HyperHDG/local_solver/diffusion_uniform_ldgh.hxx>
#include <HyperHDG/node_descriptor/cubic.hxx>
#include <HyperHDG/node_descriptor/file.hxx>
#include <HyperHDG/topology/cubic.hxx>
#include <HyperHDG/topology/file.hxx>

#include <cmath>
#include <csignal>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

/*!*************************************************************************************************
 * \brief   Parameters with constant coefficients and mixed boundary conditions.
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct ConstantParameters
{
  static constexpr std::array<unsigned int, 2U> dirichlet_nodes{1, 2};
  static constexpr std::array<unsigned int, 24U> neumann_nodes{
    3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>&,
                                               const param_float_t = 0.)
  {
    return 2.;
  }
  static param_float_t analytic_result(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  static param_float_t right_hand_side(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 1.;
  }
  static param_float_t dirichlet_value(const Point<space_dimT, param_float_t>&,
                                       const param_float_t = 0.)
  {
    return 0.;
  }
  static param_float_t neumann_value(const Point<space_dimT, param_float_t>&,
                                     const param_float_t = 0.)
  {
    return 0.;
  }
};

/*!*************************************************************************************************
 * \brief   Parameters with space dependent diffusion coefficient (all boundaries Dirichlet).
 **************************************************************************************************/
template <unsigned int space_dimT, typename param_float_t = double>
struct VariableParameters : public ConstantParameters<space_dimT, param_float_t>
{
  static constexpr std::array<unsigned int, 26U> dirichlet_nodes{
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26};
  static constexpr std::array<unsigned int, 0U> neumann_nodes{};
  static param_float_t inverse_diffusion_coeff(const Point<space_dimT, param_float_t>& point,
                                               const param_float_t = 0.)
  {
    return 1. + point[0] + 2. * point[0] * point[0];
  }
};

/*!*************************************************************************************************
 * \brief   Dense matrix from triplets (summing up duplicates), stored column-major.
 **************************************************************************************************/
template <typename sparse_matT>
vector<double> dense_from_triplets(sparse_matT& mat, const unsigned int n)
{
  vector<double> dense(n * n, 0.);
  for (unsigned int k = 0; k < mat.get_values().size(); ++k)
    dense[mat.get_cols()[k] * n + mat.get_rows()[k]] += mat.get_values()[k];
  return dense;
}

/*!*************************************************************************************************
 * \brief   Compare assembled matrix to matrix obtained by applying trace_to_flux to unit vectors.
 **************************************************************************************************/
template <typename ProblemT>
bool check(ProblemT& problem, const string& name, const bool check_factorize_once)
{
  const unsigned int n = problem.size_of_system();
  vector<double> reference(n * n, 0.);
  for (unsigned int col = 0; col < n; ++col)
  {
    vector<double> unit_vector = problem.zero_vector();
    unit_vector[col] = 1.;
    vector<double> column = problem.trace_to_flux(unit_vector);
    for (unsigned int row = 0; row < n; ++row)
      reference[col * n + row] = column[row];
  }

  double max_entry = 0.;
  for (double entry : reference)
    max_entry = max(max_entry, abs(entry));

  bool successful = true;
  for (const bool factorize_once : {false, true})
  {
    if (factorize_once && !check_factorize_once)
      continue;
    auto sparse = problem.trace_to_flux_mat(0., factorize_once);
    successful &= (sparse.get_values().size() == sparse.get_rows().size() &&
                   sparse.get_values().size() == sparse.get_cols().size());
    vector<double> dense = dense_from_triplets(sparse, n);
    double max_diff = 0.;
    for (unsigned int k = 0; k < n * n; ++k)
      max_diff = max(max_diff, abs(dense[k] - reference[k]));
    const bool ok = max_entry > 0. && max_diff <= 1e-12 * max_entry;
    cout << name << " (factorize_once = " << factorize_once << "): size " << n
         << ", max. deviation " << max_diff << ", max. entry " << max_entry
         << (ok ? "  OK" : "  FAILED") << endl;
    successful &= ok;
  }
  return successful;
}

int main(int argc, char* argv[])
{
  if (argc > 1 && string(argv[1]) == "misuse")
  {
    // hy_assert aborts. Turn this into a regular failure (exit code 1) for CTest's WILL_FAIL.
    signal(SIGABRT, [](int) { _Exit(1); });
    GlobalLoop::Elliptic<Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>,
                         NodeDescriptor::Cubic<2, 2>,
                         LocalSolver::Diffusion<2, 1, 2, VariableParameters, double> >
      problem(vector<unsigned int>{4, 4});
    problem.trace_to_flux_mat(0., true);  // Must abort in debug mode.
#ifdef NDEBUG
    cout << "Compiled with NDEBUG, the misuse of factorize_once cannot be detected." << endl;
    return 1;  // Test is expected to fail, cf. CMakeLists.txt.
#endif
    return 0;
  }

  bool successful = true;

  {
    GlobalLoop::Elliptic<Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>,
                         NodeDescriptor::Cubic<2, 2>,
                         LocalSolver::Diffusion<2, 2, 4, ConstantParameters, double> >
      problem(vector<unsigned int>{3, 2});
    successful &= check(problem, "Diffusion 2D, uniform, mixed boundary", true);
  }
  {
    GlobalLoop::Elliptic<Topology::Cubic<3, 3>, Geometry::UnitCube<3, 3, double>,
                         NodeDescriptor::Cubic<3, 3>,
                         LocalSolver::Diffusion<3, 1, 2, ConstantParameters, double> >
      problem(vector<unsigned int>{2, 2, 1});
    successful &= check(problem, "Diffusion 3D, uniform, mixed boundary", true);
  }
  {
    GlobalLoop::Elliptic<Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>,
                         NodeDescriptor::Cubic<2, 2>,
                         LocalSolver::Diffusion<2, 1, 4, VariableParameters, double> >
      problem(vector<unsigned int>{3, 3});
    successful &= check(problem, "Diffusion 2D, space dependent coefficient", false);
  }
  {
    GlobalLoop::Elliptic<Topology::File<1, 2>, Geometry::File<1, 2>, NodeDescriptor::File<1, 2>,
                         LocalSolver::Diffusion<1, 2, 4, VariableParameters, double> >
      problem(string("domains/simplex_1_2.geo"));
    successful &= check(problem, "Diffusion on non-uniform hypergraph", false);
  }
  {
    GlobalLoop::Elliptic<Topology::File<1, 3>, Geometry::File<1, 3>, NodeDescriptor::File<1, 3>,
                         LocalSolver::DiffusionUniform<1, 1, 2, double> >
      problem(string("domains/simplex_1_3.geo"), 1.);
    successful &= check(problem, "DiffusionUniform (fallback without trace_to_flux_mat)", true);
  }

  return successful ? 0 : 1;
}
