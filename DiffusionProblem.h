/*!*************************************************************************************************************
 * @brief	This is an example problem.
 * 
 * This class contains functions to define and solve Poisson's equation, i.e.,
 * @f[
 *  - \Delta u = 0 \quad \text{ in } \Omega, \qquad u = u_\text D \quad \text{ on } \partial \Omega
 * @f]
 * in a spatial domain @f$\Omega \subset \mathsf R^d@f$. Here, @f$d@f$ is the spatial dimension @c space_dim,
 * @f$\Omega@f$ is a regular graph (@c hyperedge_dim = 1) or hypergraph whose hyperedges are surfaces (@c
 * hyperedge_dim = 2) or volumes (@c hyperedge_dim = 3).
 * 
 * \todo	The loop in matrix_vector_multiply() only combines properties of HyperGraph with local solvers,
 * 			right? Dirichlet boundary conditions? Post filtering!
 * 
 * @tparam	hyperedge_dim		Dimension of a hyperedge, i.e., 1 is for PDEs defined on graphs, 2 is for PDEs
 * 								defined on surfaces, and 3 is for PDEs defined on volumes.
 * @tparam	space_dim			The dimension of the space, the object is located in. This number should be
 * 								larger than or equal to hyperedge_dim.
 * @tparam	polynomial_degree	The polynomial degree of test and trial functions.
 * 
 * @authors	Guido Kanschat, University of Heidelberg, 2019.
 * @authors	Andreas Rupp, University of Heidelberg, 2019.
 **************************************************************************************************************/

#ifndef DIFFUSIONPROBLEM_H
#define DIFFUSIONPROBLEM_H

#include "HyperGraph.h"
#include "DiffusionSolver.h"

template <unsigned int hyperedge_dim, unsigned int space_dim, unsigned int polynomial_degree>
class DiffusionProblemRegular
{
  private:
    HyperGraph < local_dof_amount_node(hyperedge_dim, polynomial_degree),
                 HyperGraph_Cubic< hyperedge_dim, space_dim >,
                 HyperGraph_Cubic_UnitCube< hyperedge_dim, space_dim > >
               hyper_graph_;
    std::vector<unsigned int> dirichlet_indices_;
    DiffusionSolver_RegularQuad<hyperedge_dim, polynomial_degree, 2 * polynomial_degree> local_solver_;
  public:
	/*!*********************************************************************************************************
     * @brief	Example problem constructor.
	 * 
	 * Constructor for class containing a HyperGraph and a DiffusionSolver that solve a diffusion problem on a
	 * regular mesh.
	 * 
	 * @param 	num_elements A std::vector containing the amount of mesh elements per spatial dimension.
	 **********************************************************************************************************/
    DiffusionProblemRegular(std::vector<int> num_elements);
    /*!*********************************************************************************************************
     * @brief	Read indices of Dirichlet type hypernodes/faces.
	 * 
	 * Read the indices ot the hypernodes/faces that are of Dirichlet type and therefore do not contain degrees
	 * of freedom that are allowed to change within iterations of the iterative solver and other processes. In
	 * contrast, these degrees of freedom are set by the user.
	 * 
	 * The user creates a vector that contains the coefficients of the corresponding degrees of freedom (read by
	 * this function) and defines the Dirichlet values by this choice. The remaining elements of the global
	 * vector of unknowns (which is @b not the vector @c indices are supposed to be zero).
	 * 
	 * @param 	indices	A std::vector containing the (global) indices of Dirichlet type hypernodes/faces.
	 **********************************************************************************************************/
    void read_dirichlet_indices(std::vector<int> indices);
    /*!*********************************************************************************************************
     * @brief	Returns vector of appropriate size for the predefined problem.
	 * 
	 * Returns a vector containing only the value zero, but of the size @f$n@f$ which is also the number which
	 * is returned if \c size_of_system() is evaluated.
	 * 
	 * @retval	zero	A @c std::vector of the correct size for the unknowns of the given problem.
	 **********************************************************************************************************/
    std::vector<double> return_zero_vector();
    /*!*********************************************************************************************************
     * @brief	Evaluate condensed matrix-vector product.
	 * 
	 * Function that evaluates the condensed, matrix-free version of the matrix-vector product @f$A x = y@f$,
	 * where @f$A@f$ is the condensed matrix of the LDG-H method, @f$x@f$ is the vector of parameters to define
	 * the skeletal variable @f$\lambda@f$, and @f$y@f$ is the resulting vector, which has the same size as the
	 * input vector @f$x@f$. 
	 * 
	 * @param	x_vec	A @c std::vector containing the input vector @f$x@f$.
	 * @retval	y_vec	A @c std::vector containing the product @f$y = Ax@f$.
	 **********************************************************************************************************/
    std::vector<double> matrix_vector_multiply(std::vector<double> x_vec);
    /*!*********************************************************************************************************
     * @brief	Determine size of condensed system for the skeletal unknowns.
	 * 
	 * Function that returns the size @f$n@f$ of the @f$n \times n@f$ linear, sparse system @f$Ax = b@f$ that is
	 * solved by the program in a matrix-free fashion.
	 * This function is needed to define a @c LinearOperator from Python's @c scipy.sparse.linalg package which
	 * can be used to define iterative solvers for sparse systems.
	 * 
	 * @retval	n		An @c int which Python needs and actually is a parsed @c unsigned @c int.
	 **********************************************************************************************************/
    int size_of_system();
    /*!*********************************************************************************************************
     * @brief	Plot solution in vtu format.
	 * 
	 * Function that plots the solution of the problem to a predefined file.
	 * 
	 * @param	lambda	A vector of unknowns containing the 
	 * @retval	file	A file in the output directory.
	 **********************************************************************************************************/
    void plot_solution(std::vector<double> lambda);
};

#endif
