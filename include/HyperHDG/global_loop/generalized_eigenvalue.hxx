#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/dense_la.hxx>
#include <HyperHDG/global_loop/prototype.hxx>
#include <HyperHDG/hdg_hypergraph.hxx>
#include <HyperHDG/hy_assert.hxx>
#include <HyperHDG/hypercube.hxx>
#include <HyperHDG/plot.hxx>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Generalized eigenvalue problem in the bulk and skeletal unknowns of an HDG method.
 *
 * This global loop assembles the stiffness matrix \f$A\f$ and the mass matrix \f$M\f$ of the
 * generalized eigenvalue problem \f$A x = \lambda M x\f$, where \f$x\f$ contains the skeletal
 * unknowns and the bulk unknowns. Its typical local solver is LocalSolver::DiffusionEigsWeighted.
 * In contrast to the condensed eigenvalue loops (NonlinearEigenvalue, ShiftedInverseEigenvalue),
 * the bulk unknowns are kept, which makes the problem linear in the eigenvalue and allows for
 * mass on the skeleton (e.g. on hypernodes). Both matrices are symmetric, \f$A\f$ is positive
 * definite (after removal of the Dirichlet unknowns), and \f$M\f$ is positive semi-definite. Thus,
 * the problem can be solved by sparse eigenvalue solvers in shift-and-invert mode, e.g. by
 * \c scipy.sparse.linalg.eigsh.
 *
 * The unknowns are ordered as follows: The first n_skeleton_dofs() unknowns are the skeletal
 * unknowns (in the order of the hypernode factory), followed by n_loc_bulk_dofs() unknowns for
 * every hyperedge (in the order of the hyperedges). The matrices are returned in coordinate format
 * (entries with the same row and column index have to be summed up), and the Dirichlet unknowns
 * (cf. dirichlet_indices) have to be removed by the user.
 *
 * Moreover, the global loop provides the values of a vector of unknowns at tensorial points of the
 * hyperedges and of their faces together with the physical coordinates of these points. Together
 * with Gaussian quadrature, this allows to evaluate integrals of discrete functions (e.g. of
 * eigenfunctions) times arbitrary functions, e.g., to compare discrete and exact eigenfunctions.
 *
 * \tparam  TopologyT       Class type containing topological information.
 * \tparam  GeometryT       Class type containing geometrical information.
 * \tparam  NodeDescriptorT Class type containing the information of nodes of hyperedges.
 * \tparam  LocalSolverT    Class type of the local solver.
 * \tparam  LargeVecT       Class type of large, global vector.
 * \tparam  dof_index_t     Index type of hyperedges. Default is \c unsigned \c int.
 *
 * \authors   Andreas Rupp, Saarland University, 2026.
 **************************************************************************************************/
template <class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename LargeVecT = std::vector<double>,
          typename dof_index_t = unsigned int>
class GeneralizedEigenvalue
{
  /*!***********************************************************************************************
   * \brief   Floating type is determined by floating type of large vector's entries.
   ************************************************************************************************/
  using dof_value_t = typename LargeVecT::value_type;
  /*!***********************************************************************************************
   * \brief   Dimension of a hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int hyEdge_dim = TopologyT::hyEdge_dim();
  /*!***********************************************************************************************
   * \brief   Dimension of the surrounding space.
   ************************************************************************************************/
  static constexpr unsigned int space_dim = TopologyT::space_dim();
  /*!***********************************************************************************************
   * \brief   Number of skeletal unknowns per hypernode.
   ************************************************************************************************/
  static constexpr unsigned int n_dofs_per_node = LocalSolverT::n_glob_dofs_per_node();
  /*!***********************************************************************************************
   * \brief   Number of bulk unknowns per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_bulk_dofs = LocalSolverT::n_loc_bulk_dofs();
  /*!***********************************************************************************************
   * \brief   Number of local unknowns per hyperedge.
   ************************************************************************************************/
  static constexpr unsigned int n_loc_dofs = LocalSolverT::n_loc_dofs();

 private:
  /*!***********************************************************************************************
   * \brief   Instantiation of a hypergraph.
   ************************************************************************************************/
  HDGHyperGraph<n_dofs_per_node,
                TopologyT,
                GeometryT,
                NodeDescriptorT,
                typename LocalSolverT::data_type>
    hyper_graph_;
  /*!***********************************************************************************************
   * \brief   Instantiation of a local solver.
   ************************************************************************************************/
  const LocalSolverT local_solver_;
  /*!***********************************************************************************************
   * \brief   Struct encoding the options for plotting.
   ************************************************************************************************/
  PlotOptions plot_options;
  /*!***********************************************************************************************
   * \brief   Global indices of the local unknowns of a hyperedge.
   ************************************************************************************************/
  template <typename hyEdgeT>
  std::array<dof_index_t, n_loc_dofs> local_to_global(hyEdgeT& hyper_edge,
                                                      const dof_index_t hyEdge_index) const
  {
    std::array<dof_index_t, n_loc_dofs> indices;
    const dof_index_t first_bulk = n_skeleton_dofs() + hyEdge_index * n_bulk_dofs;
    for (unsigned int i = 0; i < n_bulk_dofs; ++i)
      indices[i] = first_bulk + i;
    const auto hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
    std::array<dof_index_t, n_dofs_per_node> dof_indices;
    for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
    {
      hyper_graph_.hyNode_factory().get_dof_indices(hyEdge_hyNodes[face], dof_indices);
      for (unsigned int l = 0; l < n_dofs_per_node; ++l)
        indices[n_bulk_dofs + face * n_dofs_per_node + l] = dof_indices[l];
    }
    return indices;
  }
  /*!***********************************************************************************************
   * \brief   Assemble a global matrix from the local matrices returned by a member function.
   ************************************************************************************************/
  template <typename LocalMatrixFunT>
  sparse_mat<LargeVecT> assemble(LocalMatrixFunT local_matrix)
  {
    sparse_mat<LargeVecT> matrix;
    matrix.row_vec.reserve(hyper_graph_.n_hyEdges() * n_loc_dofs * n_loc_dofs);
    matrix.col_vec.reserve(hyper_graph_.n_hyEdges() * n_loc_dofs * n_loc_dofs);
    matrix.value_vec.reserve(hyper_graph_.n_hyEdges() * n_loc_dofs * n_loc_dofs);
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      auto hyper_edge = hyper_graph_[index];
      const auto loc_mat = local_matrix(hyper_edge);
      const auto indices = local_to_global(hyper_edge, index);
      for (unsigned int j = 0; j < n_loc_dofs; ++j)
        for (unsigned int i = 0; i < n_loc_dofs; ++i)
          if (loc_mat(i, j) != 0.)
          {
            matrix.row_vec.push_back(indices[i]);
            matrix.col_vec.push_back(indices[j]);
            matrix.value_vec.push_back(loc_mat(i, j));
          }
    }
    return matrix;
  }
  /*!***********************************************************************************************
   * \brief   Reference point of a face's tensorial point.
   ************************************************************************************************/
  static Point<hyEdge_dim, dof_value_t> face_reference_point(const unsigned int face,
                                                             const unsigned int point,
                                                             const std::vector<double>& abscissas)
  {
    Point<hyEdge_dim, dof_value_t> ref_point;
    unsigned int index = point;
    for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
      if (dim == face / 2)
        ref_point[dim] = (face % 2 == 0) ? 0. : 1.;
      else
      {
        ref_point[dim] = abscissas[index % abscissas.size()];
        index /= abscissas.size();
      }
    return ref_point;
  }
  /*!***********************************************************************************************
   * \brief   Reference point of a hyperedge's tensorial point.
   ************************************************************************************************/
  static Point<hyEdge_dim, dof_value_t> bulk_reference_point(const unsigned int point,
                                                             const std::vector<double>& abscissas)
  {
    Point<hyEdge_dim, dof_value_t> ref_point;
    unsigned int index = point;
    for (unsigned int dim = 0; dim < hyEdge_dim; ++dim)
    {
      ref_point[dim] = abscissas[index % abscissas.size()];
      index /= abscissas.size();
    }
    return ref_point;
  }
  /*!***********************************************************************************************
   * \brief   Integer power.
   ************************************************************************************************/
  static unsigned int int_pow(const unsigned int base, const unsigned int exponent)
  {
    unsigned int result = 1;
    for (unsigned int i = 0; i < exponent; ++i)
      result *= base;
    return result;
  }

 public:
  /*!***********************************************************************************************
   * \brief   Constructor for class containing a HyperGraph and a local solver.
   *
   * \param   construct_topo    Information to construct a topology.
   * \param   construct_geom    Information to construct a geometry.
   * \param   construct_loc_sol Information to construct a local solver.
   ************************************************************************************************/
  GeneralizedEigenvalue(const typename TopologyT::constructor_value_type& construct_topo,
                        const typename GeometryT::constructor_value_type& construct_geom,
                        const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo, construct_geom), local_solver_(construct_loc_sol)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Constructor for class containing a HyperGraph and a local solver.
   *
   * \param   construct_topo    Information to construct a topology.
   * \param   construct_loc_sol Information to construct a local solver.
   ************************************************************************************************/
  GeneralizedEigenvalue(const typename TopologyT::constructor_value_type& construct_topo,
                        const typename LocalSolverT::constructor_value_type& construct_loc_sol)
  : hyper_graph_(construct_topo), local_solver_(construct_loc_sol)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Constructor for class containing a HyperGraph and a local solver.
   *
   * \param   construct_topo    Information to construct a topology.
   ************************************************************************************************/
  GeneralizedEigenvalue(const typename TopologyT::constructor_value_type& construct_topo)
  : hyper_graph_(construct_topo)
  {
    static_assert(TopologyT::hyEdge_dim() == GeometryT::hyEdge_dim(),
                  "Hyperedge dimension of topology and geometry must be equal!");
    static_assert(TopologyT::space_dim() == GeometryT::space_dim(),
                  "Space dimension of topology and geometry must be equal!");
    static_assert(TopologyT::hyEdge_dim() == LocalSolverT::hyEdge_dim(),
                  "Hyperedge dimension of hypergraph and local solver must be equal!");
  }
  /*!***********************************************************************************************
   * \brief   Number of skeletal unknowns (including the Dirichlet unknowns).
   ************************************************************************************************/
  dof_index_t n_skeleton_dofs() const { return hyper_graph_.n_global_dofs(); }
  /*!***********************************************************************************************
   * \brief   Number of hyperedges.
   ************************************************************************************************/
  dof_index_t n_hyEdges() const { return hyper_graph_.n_hyEdges(); }
  /*!***********************************************************************************************
   * \brief   Number of bulk unknowns per hyperedge.
   ************************************************************************************************/
  dof_index_t n_loc_bulk_dofs() const { return n_bulk_dofs; }
  /*!***********************************************************************************************
   * \brief   Number of all unknowns, i.e., size of the matrices (including Dirichlet unknowns).
   ************************************************************************************************/
  dof_index_t size_of_system() const { return n_skeleton_dofs() + n_hyEdges() * n_bulk_dofs; }
  /*!***********************************************************************************************
   * \brief   Return indices of the (skeletal) Dirichlet unknowns.
   ************************************************************************************************/
  std::vector<unsigned int> dirichlet_indices()
  {
    std::vector<unsigned int> dirichlet_indices;
    std::array<dof_index_t, n_dofs_per_node> dof_indices;
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      auto hyper_edge = hyper_graph_[index];
      const auto hyNode_types = local_solver_.node_types(hyper_edge);
      const auto hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
      for (unsigned int face = 0; face < hyNode_types.size(); ++face)
        if (hyNode_types[face] == 1)
        {
          hyper_graph_.hyNode_factory().get_dof_indices(hyEdge_hyNodes[face], dof_indices);
          for (unsigned int l = 0; l < dof_indices.size(); ++l)
            dirichlet_indices.push_back(dof_indices[l]);
        }
    }
    std::sort(dirichlet_indices.begin(), dirichlet_indices.end());
    dirichlet_indices.erase(std::unique(dirichlet_indices.begin(), dirichlet_indices.end()),
                            dirichlet_indices.end());
    return dirichlet_indices;
  }
  /*!***********************************************************************************************
   * \brief   Stiffness matrix in coordinate format.
   ************************************************************************************************/
  sparse_mat<LargeVecT> stiffness_mat()
  {
    return assemble([&](auto& hyper_edge) { return local_solver_.stiffness_matrix(hyper_edge); });
  }
  /*!***********************************************************************************************
   * \brief   Mass matrix in coordinate format.
   ************************************************************************************************/
  sparse_mat<LargeVecT> mass_mat()
  {
    return assemble([&](auto& hyper_edge) { return local_solver_.mass_matrix(hyper_edge); });
  }
  /*!***********************************************************************************************
   * \brief   Coefficients of all hyperedges (n_hyEdges() times the number of coefficients).
   ************************************************************************************************/
  LargeVecT coefficients()
    requires requires { LocalSolverT::n_coefficients(); }
  {
    LargeVecT result;
    result.reserve(hyper_graph_.n_hyEdges() * LocalSolverT::n_coefficients());
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      auto hyper_edge = hyper_graph_[index];
      const auto coeffs = LocalSolverT::coefficients(hyper_edge);
      result.insert(result.end(), coeffs.begin(), coeffs.end());
    }
    return result;
  }
  /*!***********************************************************************************************
   * \brief   Values of the bulk unknowns at tensorial points of all hyperedges.
   *
   * \param   x_vec         Vector of all unknowns (of length size_of_system()).
   * \param   abscissas     One-dimensional abscissas in [0,1] defining the tensorial points.
   * \retval  values        Value at point p of hyperedge e at p + e * (number of points).
   ************************************************************************************************/
  LargeVecT bulk_values(const LargeVecT& x_vec, const std::vector<double>& abscissas)
  {
    hy_assert(x_vec.size() == size_of_system(), "The vector must contain all unknowns!");
    const unsigned int n_points = int_pow(abscissas.size(), hyEdge_dim);
    LargeVecT values(hyper_graph_.n_hyEdges() * n_points);
    SmallVec<n_bulk_dofs, dof_value_t> coeffs;
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      for (unsigned int i = 0; i < n_bulk_dofs; ++i)
        coeffs[i] = x_vec[n_skeleton_dofs() + index * n_bulk_dofs + i];
      for (unsigned int pt = 0; pt < n_points; ++pt)
        values[index * n_points + pt] =
          LocalSolverT::bulk_value(coeffs, bulk_reference_point(pt, abscissas));
    }
    return values;
  }
  /*!***********************************************************************************************
   * \brief   Physical coordinates of the tensorial points used in bulk_values.
   *
   * \retval  points        Coordinate k of point p of hyperedge e at k + space_dim * (p + e *
   *(number of points)).
   ************************************************************************************************/
  LargeVecT bulk_points(const std::vector<double>& abscissas)
  {
    const unsigned int n_points = int_pow(abscissas.size(), hyEdge_dim);
    LargeVecT points(hyper_graph_.n_hyEdges() * n_points * space_dim);
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      auto hyper_edge = hyper_graph_[index];
      for (unsigned int pt = 0; pt < n_points; ++pt)
      {
        const auto point = hyper_edge.geometry.map_ref_to_phys(bulk_reference_point(pt, abscissas));
        for (unsigned int k = 0; k < space_dim; ++k)
          points[(index * n_points + pt) * space_dim + k] = point[k];
      }
    }
    return points;
  }
  /*!***********************************************************************************************
   * \brief   Measures (lengths, areas, volumes) of all hyperedges.
   ************************************************************************************************/
  LargeVecT hyEdge_measures()
  {
    LargeVecT measures(hyper_graph_.n_hyEdges());
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
      measures[index] = hyper_graph_[index].geometry.area();
    return measures;
  }
  /*!***********************************************************************************************
   * \brief   Values of the skeletal unknowns at tensorial points of all faces of all hyperedges.
   *
   * Faces that are shared by several hyperedges appear several times (once per hyperedge).
   *
   * \param   x_vec         Vector of all unknowns (or of the skeletal unknowns only).
   * \param   abscissas     One-dimensional abscissas in [0,1] defining the tensorial points.
   * \retval  values        Value at point p of face f of hyperedge e at p + (number of points) *
   *                        (f + 2 * hyEdge_dim * e).
   ************************************************************************************************/
  LargeVecT face_values(const LargeVecT& x_vec, const std::vector<double>& abscissas)
  {
    hy_assert(x_vec.size() >= n_skeleton_dofs(), "The vector must contain the skeletal unknowns!");
    const unsigned int n_points = int_pow(abscissas.size(), hyEdge_dim - 1);
    LargeVecT values(hyper_graph_.n_hyEdges() * 2 * hyEdge_dim * n_points);
    std::array<dof_value_t, n_dofs_per_node> dof_values;
    SmallVec<n_dofs_per_node, dof_value_t> coeffs;
    Point<std::max(hyEdge_dim - 1, 1U), dof_value_t> face_point;
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      auto hyper_edge = hyper_graph_[index];
      const auto hyEdge_hyNodes = hyper_edge.topology.get_hyNode_indices();
      for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
      {
        hyper_graph_.hyNode_factory().get_dof_values(hyEdge_hyNodes[face], x_vec, dof_values);
        for (unsigned int l = 0; l < n_dofs_per_node; ++l)
          coeffs[l] = dof_values[l];
        for (unsigned int pt = 0; pt < n_points; ++pt)
        {
          unsigned int pt_index = pt;
          for (unsigned int dim = 0; dim < hyEdge_dim - 1; ++dim)
          {
            face_point[dim] = abscissas[pt_index % abscissas.size()];
            pt_index /= abscissas.size();
          }
          values[(index * 2 * hyEdge_dim + face) * n_points + pt] =
            LocalSolverT::face_value(coeffs, face_point);
        }
      }
    }
    return values;
  }
  /*!***********************************************************************************************
   * \brief   Physical coordinates of the tensorial points used in face_values.
   ************************************************************************************************/
  LargeVecT face_points(const std::vector<double>& abscissas)
  {
    const unsigned int n_points = int_pow(abscissas.size(), hyEdge_dim - 1);
    LargeVecT points(hyper_graph_.n_hyEdges() * 2 * hyEdge_dim * n_points * space_dim);
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      auto hyper_edge = hyper_graph_[index];
      for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
        for (unsigned int pt = 0; pt < n_points; ++pt)
        {
          const auto point =
            hyper_edge.geometry.map_ref_to_phys(face_reference_point(face, pt, abscissas));
          for (unsigned int k = 0; k < space_dim; ++k)
            points[((index * 2 * hyEdge_dim + face) * n_points + pt) * space_dim + k] = point[k];
        }
    }
    return points;
  }
  /*!***********************************************************************************************
   * \brief   Node types (as given in the domain) of all faces of all hyperedges.
   ************************************************************************************************/
  std::vector<unsigned int> face_node_types()
  {
    std::vector<unsigned int> types(hyper_graph_.n_hyEdges() * 2 * hyEdge_dim);
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      auto hyper_edge = hyper_graph_[index];
      for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
        types[index * 2 * hyEdge_dim + face] = hyper_edge.node_descriptor[face];
    }
    return types;
  }
  /*!***********************************************************************************************
   * \brief   Hypernode indices of all faces of all hyperedges.
   ************************************************************************************************/
  std::vector<unsigned int> face_node_indices()
  {
    std::vector<unsigned int> nodes(hyper_graph_.n_hyEdges() * 2 * hyEdge_dim);
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      const auto hyEdge_hyNodes = hyper_graph_.hyEdge_topology(index).get_hyNode_indices();
      for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
        nodes[index * 2 * hyEdge_dim + face] = hyEdge_hyNodes[face];
    }
    return nodes;
  }
  /*!***********************************************************************************************
   * \brief   Measures (lengths, areas) of all faces of all hyperedges.
   ************************************************************************************************/
  LargeVecT face_measures()
  {
    LargeVecT measures(hyper_graph_.n_hyEdges() * 2 * hyEdge_dim);
    for (dof_index_t index = 0; index < hyper_graph_.n_hyEdges(); ++index)
    {
      auto hyper_edge = hyper_graph_[index];
      for (unsigned int face = 0; face < 2 * hyEdge_dim; ++face)
        measures[index * 2 * hyEdge_dim + face] = hyper_edge.geometry.face_area(face);
    }
    return measures;
  }
  /*!***********************************************************************************************
   * \brief   Set plot option and return old plot option.
   ************************************************************************************************/
  std::string plot_option(const std::string& option, std::string value = "")
  {
    return set_plot_option(plot_options, option, value);
  }
  /*!***********************************************************************************************
   * \brief   Plot an eigenfunction in vtu format.
   *
   * The bulk values are reconstructed from the skeletal values and the eigenvalue by the local
   * solver, cf. LocalSolver::DiffusionEigsWeighted::bulk_values.
   *
   * \param   x_vec         A vector containing (at least) the skeletal unknowns.
   * \param   eigenvalue    The eigenvalue.
   ************************************************************************************************/
  void plot_solution(const std::vector<dof_value_t>& x_vec, const dof_value_t eigenvalue = 0.)
  {
    std::vector<dof_value_t> skeleton(x_vec.begin(), x_vec.begin() + n_skeleton_dofs());
    plot(hyper_graph_, local_solver_, skeleton, plot_options, eigenvalue);
  }
  /*!***********************************************************************************************
   * \brief   Return refinement level.
   ************************************************************************************************/
  unsigned int get_refinement() { return hyper_graph_.get_refinement(); }
  /*!***********************************************************************************************
   * \brief   Set refinement level.
   ************************************************************************************************/
  void set_refinement(unsigned int level) { hyper_graph_.set_refinement(level); }
};  // end of class GeneralizedEigenvalue

}  // end of namespace GlobalLoop
