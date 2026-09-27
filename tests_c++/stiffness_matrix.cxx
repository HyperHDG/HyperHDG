/*!*************************************************************************************************
 * \file    stiffness_matrix.cxx
 * \brief   Check the assembly of the condensed (stiffness) matrix and the operator cache.
 *
 * The matrix assembled by trace_to_flux_mat (with the local matrix being computed per hyperedge
 * and with the local matrix being computed only once) and the cached operator (cache_operator) are
 * compared to the matrix whose columns are obtained by applying the matrix-free operator
 * trace_to_flux to all unit vectors. This is done
 * - for the local solver Diffusion, which provides trace_to_flux_mat (one factorization per
 *   hyperedge), on a uniform mesh with mixed boundary conditions, on a uniform mesh with a space
 *   dependent diffusion coefficient, on a non-uniform hypergraph, and on a refined hypergraph with
 *   hypernodes shared by more than two hyperedges (also after changing the refinement), and
 * - for the local solver DiffusionUniform, which does not provide trace_to_flux_mat (fallback).
 *
 * Moreover, it is checked that
 * - using factorize_once for a space dependent diffusion coefficient throws std::invalid_argument
 *   (in debug and in release mode), and
 * - exceptions of the local solver are forwarded to the caller (also if the assembly is
 *   parallelized with OpenMP, where they must not leave the parallel region).
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
#include <iostream>
#include <stdexcept>
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
 * \brief   Local solver that throws for all interior hyperedges (i.e. without boundary faces).
 **************************************************************************************************/
struct ThrowingDiffusion : public LocalSolver::Diffusion<2, 1, 2, ConstantParameters, double>
{
  using LocalSolver::Diffusion<2, 1, 2, ConstantParameters, double>::Diffusion;
  template <typename hyEdgeT>
  SmallSquareMat<2 * 2 * n_glob_dofs_per_node(), double> trace_to_flux_mat(
    hyEdgeT& hyper_edge,
    const double time = 0.) const
  {
    bool interior = true;
    for (unsigned int face = 0; face < 4; ++face)
      interior &= hyper_edge.node_descriptor[face] == 0;
    if (interior)
      throw std::runtime_error("ThrowingDiffusion");
    return LocalSolver::Diffusion<2, 1, 2, ConstantParameters, double>::trace_to_flux_mat(
      hyper_edge, time);
  }
};

/*!*************************************************************************************************
 * \brief   Check that function throws an exception of type ExceptionT.
 **************************************************************************************************/
template <typename ExceptionT, typename FunctionT>
bool throws(const string& name, FunctionT&& function)
{
  bool thrown = false;
  try
  {
    function();
  }
  catch (const ExceptionT& exc)
  {
    thrown = true;
  }
  cout << name << (thrown ? ": exception thrown  OK" : ": no exception  FAILED") << endl;
  return thrown;
}

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
 * \brief   Matrix whose columns are trace_to_flux applied to all unit vectors (column-major).
 **************************************************************************************************/
template <typename ProblemT>
vector<double> operator_columns(ProblemT& problem)
{
  const unsigned int n = problem.size_of_system();
  vector<double> matrix(n * n, 0.);
  for (unsigned int col = 0; col < n; ++col)
  {
    vector<double> unit_vector = problem.zero_vector();
    unit_vector[col] = 1.;
    vector<double> column = problem.trace_to_flux(unit_vector);
    for (unsigned int row = 0; row < n; ++row)
      matrix[col * n + row] = column[row];
  }
  return matrix;
}

/*!*************************************************************************************************
 * \brief   Compare the (cached) trace_to_flux to a reference matrix.
 **************************************************************************************************/
template <typename ProblemT>
bool check_cached(ProblemT& problem,
                  const vector<double>& reference,
                  const double max_entry,
                  const string& name)
{
  const vector<double> cached = operator_columns(problem);
  double max_diff = 0.;
  for (unsigned int k = 0; k < reference.size(); ++k)
    max_diff = max(max_diff, abs(cached[k] - reference[k]));
  const bool ok = max_entry > 0. && max_diff <= 1e-12 * max_entry;
  cout << name << ": max. deviation " << max_diff << (ok ? "  OK" : "  FAILED") << endl;
  return ok;
}

/*!*************************************************************************************************
 * \brief   Compare assembled matrix to matrix obtained by applying trace_to_flux to unit vectors.
 **************************************************************************************************/
template <typename ProblemT>
bool check(ProblemT& problem, const string& name, const bool check_factorize_once)
{
  const unsigned int n = problem.size_of_system();
  const vector<double> reference = operator_columns(problem);

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
  for (const bool factorize_once : {false, true})
  {
    if (factorize_once && !check_factorize_once)
      continue;
    problem.cache_operator(true, factorize_once);
    successful &=
      check_cached(problem, reference, max_entry,
                   name + " (cache_operator, factorize_once = " + to_string(factorize_once) + ")");
    problem.cache_operator(false);
  }
  return successful;
}

int main()
{
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
  {
    // Hypernodes of the refined cross are shared by up to four hyperedges.
    GlobalLoop::Elliptic<Topology::File<1, 2>, Geometry::File<1, 2>, NodeDescriptor::File<1, 2>,
                         LocalSolver::Diffusion<1, 2, 4, ConstantParameters, double> >
      problem(string("domains/cross_2d.geo"));
    problem.set_refinement(2);
    successful &= check(problem, "Diffusion on refined cross", false);
    // The cached operator has to be reassembled after the refinement has changed.
    problem.set_refinement(3);
    const vector<double> reference = operator_columns(problem);
    double max_entry = 0.;
    for (double entry : reference)
      max_entry = max(max_entry, abs(entry));
    problem.set_refinement(2);
    problem.cache_operator(true);
    problem.trace_to_flux(problem.zero_vector());  // Assemble the cache for refinement level 2.
    problem.set_refinement(3);
    successful &=
      check_cached(problem, reference, max_entry, "Diffusion on cross after refinement");
  }

  // factorize_once must be detected to be wrong for space dependent coefficients.
  {
    GlobalLoop::Elliptic<Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>,
                         NodeDescriptor::Cubic<2, 2>,
                         LocalSolver::Diffusion<2, 1, 2, VariableParameters, double> >
      problem(vector<unsigned int>{4, 4});
    successful &= throws<invalid_argument>("Misuse of factorize_once in trace_to_flux_mat",
                                           [&]() { problem.trace_to_flux_mat(0., true); });
    problem.cache_operator(true, true);
    successful &= throws<invalid_argument>("Misuse of factorize_once in cache_operator",
                                           [&]() { problem.trace_to_flux(problem.zero_vector()); });
  }

  // Exceptions of the local solver must be forwarded (also from parallel loops).
  {
    GlobalLoop::Elliptic<Topology::Cubic<2, 2>, Geometry::UnitCube<2, 2, double>,
                         NodeDescriptor::Cubic<2, 2>, ThrowingDiffusion>
      problem(vector<unsigned int>{8, 8});
    for (const bool factorize_once : {false, true})
    {
      const string suffix = " (factorize_once = " + to_string(factorize_once) + ")";
      successful &= throws<runtime_error>("Exception in trace_to_flux_mat" + suffix,
                                          [&]() { problem.trace_to_flux_mat(0., factorize_once); });
      problem.cache_operator(true, factorize_once);
      successful &= throws<runtime_error>("Exception in cached trace_to_flux" + suffix,
                                          [&]() { problem.trace_to_flux(problem.zero_vector()); });
      problem.cache_operator(false);
    }
  }

  return successful ? 0 : 1;
}
