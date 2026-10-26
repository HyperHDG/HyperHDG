#pragma once  // Ensure that file is included only once in a single compilation.

#include <HyperHDG/global_loop/mat_generate.hxx>
#include <HyperHDG/global_loop/prototype.hxx>

namespace GlobalLoop
{
/*!*************************************************************************************************
 * \brief   Assembly of the condensed matrix of a global loop.
 *
 * Base class (CRTP) of the linear global loops. The global loop \c GlobalLoopT, which is passed as
 * a template together with its template arguments, has to declare this class a friend, and to
 * provide the members \c hyper_graph_ and \c local_solver_. Passing the template (instead of the
 * class) allows to derive via
 * \c CondensedOperator<GlobalLoop, TopologyT, ...> without spelling out the derived class.
 *
 * \tparam  GlobalLoopT     The template of the global loop.
 * \tparam  TopologyT       Class type containing topological information.
 * \tparam  GeometryT       Class type containing geometrical information.
 * \tparam  NodeDescriptorT Class type containing the information of nodes of hyperedges.
 * \tparam  LocalSolverT    Class type of the local solver.
 * \tparam  LargeVecT       Type of large (global) vectors.
 * \tparam  dof_index_t     Index type of hyperedges.
 **************************************************************************************************/
template <template <class, class, class, class, typename, typename> class GlobalLoopT,
          class TopologyT,
          class GeometryT,
          class NodeDescriptorT,
          class LocalSolverT,
          typename LargeVecT,
          typename dof_index_t>
class CondensedOperator
{
  /*!***********************************************************************************************
   * \brief   The derived global loop.
   ************************************************************************************************/
  using GlobalLoop =
    GlobalLoopT<TopologyT, GeometryT, NodeDescriptorT, LocalSolverT, LargeVecT, dof_index_t>;
  /*!***********************************************************************************************
   * \brief   Floating type is determined by floating type of large vector's entries.
   ************************************************************************************************/
  using dof_value_t = typename LargeVecT::value_type;

 public:
  /*!***********************************************************************************************
   * \brief   Assemble the condensed matrix.
   *
   * Function that assembles the condensed matrix \f$A\f$ of the HDG method, i.e., the matrix of
   * the linear map evaluated by \c trace_to_flux, in coordinate (triplet) format. Entries with the
   * same row and column index have to be summed up (as done by SciPy's \c csr_matrix).
   *
   * If the local solver provides \c trace_to_flux_mat, its local matrix is factorized once per
   * hyperedge. Otherwise, the local problem is solved for every local unit vector.
   *
   * \tparam  hyNode_index_t  Typename of the hypernode index. Defaults to \c unsigned \c int.
   * \param   param           Parameter of trace_to_flux, i.e., the time at which given analytic
   *                          functions will be evaluated (or the shift for eigenvalue problems).
   * \param   factorize_once  If true, the local matrix is only computed (and factorized) once for
   *                          every combination of face types (boundary conditions) and reused for
   *                          all hyperedges. This is only correct if all hyperedges have the same
   *                          local matrix, e.g., for uniform meshes with coefficients that are
   *                          constant in space. It is wrong for space dependent diffusion tensors
   *                          or non-uniform meshes, which is detected for a sample of hyperedges
   *                          (for all hyperedges in debug mode) and results in an exception.
   *                          Defaults to false (local matrix per hyperedge).
   * \retval  sparse_mat      Condensed matrix in coordinate format.
   ************************************************************************************************/
  template <typename hyNode_index_t = dof_index_t>
  sparse_mat<LargeVecT> trace_to_flux_mat(const dof_value_t param = 0.,
                                          const bool factorize_once = false)
  {
    GlobalLoop& loop = static_cast<GlobalLoop&>(*this);
    return MatGenerate::assemble_triplets<hyNode_index_t, LargeVecT>(
      loop.hyper_graph_, loop.local_solver_, param, factorize_once);
  }
};

}  // end of namespace GlobalLoop
